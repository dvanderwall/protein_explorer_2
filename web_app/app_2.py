"""
Flask web application for the KinoPlex protein explorer application.

This version uses database queries instead of file loading, but maintains
the same functionality and output as the original application.
"""



from flask import Flask, render_template, request, jsonify, redirect, url_for, session
import os
import sys
import logging
import networkx as nx
import requests
import numpy as np
from Bio.PDB import PDBParser, Selection, NeighborSearch
import io
import re
import json
from typing import Dict, List, Optional
import shutil
import math
import time








# Add the parent directory to the path to import protein_explorer
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import protein_explorer as pe

from protein_explorer.data.scaffold import (
    get_uniprot_id_from_gene,
    get_protein_by_id, 
    get_alphafold_structure
)


# Import database module through the package for consistency
from protein_explorer.db import (
    init_db, get_phosphosite_data, get_phosphosites_batch,
    find_structural_matches, find_structural_matches_batch,
    find_sequence_matches, find_sequence_matches_batch,
    get_kinase_scores, get_kinase_scores_batch
)

# Import phosphosite analysis functions
from protein_explorer.analysis.phospho import analyze_phosphosites
from protein_explorer.analysis.phospho_analyzer_2 import (
    get_phosphosite_data,
    enhance_phosphosite, 
    enhance_structural_matches,
    enhance_structural_matches_group,
    analyze_phosphosite_context,
    enhance_site_visualization, 
    create_comparative_motif_visualization,
    analyze_residue_distributions
)

# Import visualization functions
from protein_explorer.visualization.network import create_phosphosite_network

from protein_explorer.analysis.sequence_analyzer_2 import (
    find_sequence_matches,
    find_sequence_matches_with_connections,
    analyze_motif_conservation,
    create_sequence_network_data,
    get_motif_enrichment,
    create_sequence_motif_visualization
)

from protein_explorer.analysis.kinase_predictor_2 import (
    load_kinase_scores, 
    get_site_kinase_scores,
    predict_kinases,
    get_heatmap_data,
    get_kinase_comparison_data,
    get_known_kinase_info,
    categorize_kinases_by_family
)

# Import the new network-based functions
from protein_explorer.analysis.network_kinase_predictor_2 import (
    get_similar_sites, compute_aggregated_kinase_scores, predict_kinases_network,
    get_network_heatmap_data, get_network_kinase_comparison, 
    get_kinase_family_distribution_network
)
from protein_explorer.visualization.protein_phosphosite_network import create_phosphosite_network_visualization
from protein_explorer.visualization.protein_sequence_phosphosite_network import create_sequence_network_visualization

from protein_explorer.analysis.Cantley_Kinase_Scorer import (
    process_matches_for_kinase_scoring,
    create_heatmap_visualization,
    create_top_kinases_visualization,
    create_protein_kinase_report,
    get_html_protein_kinase_report,
    create_site_selector_ui,
    process_site_specific_matches,
    generate_site_analysis_visualizations
)

MATCH_CACHE = {}

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize database connection
init_db()

# Ensure cache directory exists
def ensure_cache_directory():
    """Ensure the cache directory exists at application startup."""
    # Paths for cache
    cache_dir = os.path.expanduser("~/.protein_explorer/cache")
    
    # Create cache directory if it doesn't exist
    if not os.path.exists(cache_dir):
        try:
            print(f"Creating cache directory: {cache_dir}")
            os.makedirs(cache_dir, exist_ok=True)
            print(f"Cache directory created successfully: {cache_dir}")
        except Exception as e:
            # If we can't create in home directory, use a temporary directory
            import tempfile
            alt_cache_dir = os.path.join(tempfile.gettempdir(), "protein_explorer_cache")
            print(f"Failed to create cache in home directory: {e}")
            print(f"Using alternative cache directory: {alt_cache_dir}")
            
            try:
                os.makedirs(alt_cache_dir, exist_ok=True)
            except Exception as e3:
                print(f"Failed to create alternative cache directory: {e3}")
                print("Application may have issues with caching")

# Run cache directory initialization
cache_dir = os.path.expanduser("~/.protein_explorer/cache")
if os.path.exists(cache_dir):
    shutil.rmtree(cache_dir)
    os.makedirs(cache_dir, exist_ok=True)
    print("Cache cleared.")
else:
    print("Cache directory does not exist.")
ensure_cache_directory()

app = Flask(__name__)
# Add this after creating your Flask app (where you have "app = Flask(__name__)")
app.secret_key = os.environ.get('SECRET_KEY', 'dev_key_change_in_production')  # More secure to get from environment




@app.route('/')
def index():
    """Render the home page."""
    return render_template('index.html')

@app.route('/search', methods=['GET', 'POST'])
def search():
    """Search for a protein by UniProt ID or gene symbol."""
    if request.method == 'POST':
        # Get search parameters
        identifier = request.form.get('identifier', '')
        id_type = request.form.get('id_type', 'uniprot')
        
        if not identifier:
            return render_template('search.html', error="Please enter an identifier")
        
        try:
            # Redirect to protein page
            return redirect(url_for('protein', identifier=identifier, id_type=id_type))
        except Exception as e:
            logger.error(f"Error in search: {e}")
            return render_template('search.html', error=str(e))
    
    # GET request
    return render_template('search.html')

@app.route('/protein/<identifier>')
def protein(identifier):
    """Display protein information and structure with full supplementary structural data."""
    try:
        print(f"DEBUG: Processing protein view for identifier: {identifier}")
        id_type = request.args.get('id_type', 'uniprot')
        
        # Get protein data
        if id_type.lower() == 'uniprot':
            print(f"DEBUG: Getting protein by UniProt ID")
            protein_data = pe.data.get_protein_by_id(uniprot_id=identifier)
        else:
            print(f"DEBUG: Getting protein by gene symbol")
            pre_filt_id = get_uniprot_id_from_gene(identifier)
            protein_data = get_protein_by_id(pre_filt_id)
        
        print(f"DEBUG: Protein data retrieved, has_structure: {protein_data.get('has_structure', False)}")
        
        # TEMPORARY FIX: For specific proteins, manually set has_structure=True
        if identifier in ['P04637', 'P53_HUMAN']:
            protein_data['has_structure'] = True
            print(f"DEBUG: Manually set has_structure=True for {identifier}")
        
        # Create cache directory if it doesn't exist
        import os
        cache_dir = os.path.expanduser("~/.protein_explorer/cache")
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir, exist_ok=True)
            print(f"DEBUG: Created cache directory: {cache_dir}")
        
        # Check if protein has a structure
        structure_html = None
        structure = None
        phosphosite_html = None
        phosphosites = None
        
        if protein_data.get('has_structure', False):
            print(f"DEBUG: Protein has structure, retrieving...")
            
            # Get structure 
            try:
                # Try direct download first for known proteins
                if identifier in ['P04637', 'P53_HUMAN']:
                    import requests
                    url = f"https://alphafold.ebi.ac.uk/files/AF-{protein_data['uniprot_id']}-F1-model_v4.pdb"
                    response = requests.get(url)
                    if response.status_code == 200:
                        structure = response.text
                        print(f"DEBUG: Got structure directly, length: {len(structure)}")
                    else:
                        # Fall back to normal method
                        structure = pe.data.get_alphafold_structure(protein_data['uniprot_id'])
                else:
                    # Use normal method
                    structure = pe.data.get_alphafold_structure(protein_data['uniprot_id'])
                
                # Get the protein sequence from metadata if available
                sequence = None
                if protein_data.get('metadata', {}).get('sequence', {}).get('value'):
                    sequence = protein_data['metadata']['sequence']['value']
                
                # Create structure visualization
                if structure:
                    print(f"DEBUG: Creating structure visualization")
                    structure_html = pe.visualization.visualize_structure(
                        structure,
                        sequence=sequence
                    )
                    
                    # Analyze phosphorylation sites
                    if sequence:
                        print(f"DEBUG: Analyzing phosphorylation sites")
                        try:
                            # Use phospho.py for site detection
                            phosphosites = pe.analysis.phospho.analyze_phosphosites(
                                sequence, structure, uniprot_id=protein_data.get('uniprot_id')
                            )
                            print(f"DEBUG: Found {len(phosphosites)} potential phosphorylation sites")
                            
                            # Get site IDs for batch query to database
                            site_ids = [f"{protein_data['uniprot_id']}_{site['resno']}" for site in phosphosites]
                            
                            # Batch retrieve supplementary data
                            supp_data_dict = get_phosphosites_batch(site_ids)
                            
                            # Enhance phosphosites with additional metrics from database
                            for site in phosphosites:
                                site_id = f"{protein_data['uniprot_id']}_{site['resno']}"
                                supp_data = supp_data_dict.get(site_id)
                                
                                if supp_data:
                                    # Add database fields to site data
                                    for key in ['surface_accessibility', 'site_plddt', 
                                                'polar_aa_percent', 'nonpolar_aa_percent', 
                                                'acidic_aa_percent', 'basic_aa_percent',
                                                'SITE_+/-7_AA']:
                                        if key in supp_data and supp_data[key] is not None:
                                            # Map database field names to our field names
                                            if key == 'SITE_+/-7_AA':
                                                site['motif'] = supp_data[key]
                                            elif key == 'site_plddt':
                                                site['mean_plddt'] = supp_data[key]
                                            else:
                                                # Convert to snake_case
                                                site_key = key.lower()
                                                site[site_key] = supp_data[key]
                            
                            from protein_explorer.analysis.enhanced_table import enhance_phosphosite_table
                            phosphosite_html = enhance_phosphosite_table(phosphosites, protein_data['uniprot_id'])
                            
                        except Exception as e:
                            print(f"DEBUG: Error analyzing phosphosites: {e}")
                            import traceback
                            print(traceback.format_exc())
                            
                            # Fall back to basic analysis without supplementary data
                            phosphosites = pe.analysis.phospho.analyze_phosphosites(
                                sequence, structure, uniprot_id=protein_data.get('uniprot_id')
                            )
                            
                            # Use the enhanced table function even with basic data
                            try:
                                from protein_explorer.analysis.enhanced_table import enhance_phosphosite_table
                                phosphosite_html = enhance_phosphosite_table(phosphosites, protein_data['uniprot_id'])
                            except Exception as e2:
                                print(f"DEBUG: Error using enhanced table: {e2}")
                                # Fall back to original HTML generation if enhanced_table fails
                                
                                phosphosite_html = f"""
                                <div class="card mt-4">
                                    <div class="card-header">
                                        <h5 class="mb-0">Phosphorylation Site Analysis</h5>
                                    </div>
                                    <div class="card-body p-0">
                                        <div class="table-responsive">
                                            <table class="table table-striped table-hover phosphosite-table">
                                                <thead class="thead-light">
                                                    <tr>
                                                        <th>Site</th>
                                                        <th>Motif (-7 to +7)</th>
                                                        <th>Mean pLDDT</th>
                                                        <th>Nearby Residues (10Å)</th>
                                                        <th>Known in PhosphositePlus</th>
                                                    </tr>
                                                </thead>
                                                <tbody id="phosphosite-table">
                                """
                                
                                for site in phosphosites:
                                    # Add data attributes to the row for better visualization
                                    data_attrs = f'''
                                        data-site="{site['site']}" 
                                        data-resno="{site['resno']}" 
                                        data-type="{site['site'][0] if 'site' in site else ''}"
                                        data-plddt="{site.get('mean_plddt', 0)}" 
                                        data-nearby="{site.get('nearby_count', 0)}"
                                        data-known="{str(site.get('is_known', False)).lower()}"
                                    '''
                                    
                                    phosphosite_html += f"""
                                    <tr {data_attrs}>
                                        <td><a href="/site/{protein_data['uniprot_id']}/{site['site']}" class="text-decoration-none">
                                            <strong id="site-{site['resno']}">{site['site']}</strong>
                                        </a></td>
                                        <td><code class="motif-sequence">{site['motif']}</code></td>
                                        <td>{site['mean_plddt']}</td>
                                        <td>{site['nearby_count']}</td>
                                        <td>{"Yes" if site.get('is_known', False) else "No"}</td>
                                    </tr>
                                    """
                                
                                phosphosite_html += """
                                                </tbody>
                                            </table>
                                        </div>
                                    </div>
                                </div>
                                
                                <script>
                                    document.addEventListener('DOMContentLoaded', function() {
                                        const siteLinks = document.querySelectorAll('.site-link');
                                        siteLinks.forEach(link => {
                                            link.addEventListener('click', function(e) {
                                                e.preventDefault();
                                                const resno = this.getAttribute('data-resno');
                                                const sequenceSpans = document.querySelectorAll('.sequence-viewer span');
                                                if (sequenceSpans.length > 0) {
                                                    const index = parseInt(resno) - 1;
                                                    if (index >= 0 && index < sequenceSpans.length) {
                                                        sequenceSpans[index].click();
                                                    }
                                                }
                                            });
                                        });
                                    });
                                </script>
                                """
            except Exception as e:
                print(f"DEBUG: Error getting structure: {e}")
                structure_html = f'<div class="alert alert-danger">Error loading structure: {str(e)}</div>'
        
        print(f"DEBUG: Rendering template with phosphosite_html: {'Present' if phosphosite_html else 'None'}")
        
        # Make sure we include necessary scripts for visualization
        extra_scripts = """
        <script src="https://cdn.jsdelivr.net/npm/chart.js@3.9.1/dist/chart.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-annotation@2.0.1/dist/chartjs-plugin-annotation.min.js"></script>
        <script>
        // Check if Chart.js loaded successfully
        document.addEventListener('DOMContentLoaded', function() {
            console.log('DOM loaded, Chart.js status:', typeof Chart !== 'undefined' ? 'Loaded' : 'Not loaded');
            
            // If not loaded, try to load it again
            if (typeof Chart === 'undefined') {
                console.warn('Chart.js not loaded, attempting to load it again');
                const script = document.createElement('script');
                script.src = 'https://cdn.jsdelivr.net/npm/chart.js@3.9.1/dist/chart.min.js';
                script.onload = function() {
                    console.log('Chart.js loaded successfully on retry');
                    // Load visualization script
                    const vizScript = document.createElement('script');
                    vizScript.src = "{{ url_for('static', filename='js/phosphosite-visualization.js') }}";
                    document.head.appendChild(vizScript);
                };
                document.head.appendChild(script);
            }
        });
        </script>
        """
        
        return render_template(
            'protein.html',
            protein=protein_data,
            structure_html=structure_html,
            phosphosites=phosphosites,  # Pass the phosphosites data
            phosphosite_html=phosphosite_html,  # Pass the HTML for rendering
            extra_scripts=extra_scripts  # Add extra scripts for visualization
        )
    except Exception as e:
        print(f"DEBUG: Exception in protein view: {e}")
        logger.error(f"Error in protein view: {e}")
        return render_template('error.html', error=str(e))

@app.route('/analyze', methods=['GET', 'POST'])
def analyze():
    """Analyze multiple proteins."""
    if request.method == 'POST':
        # Get proteins from form
        proteins_text = request.form.get('proteins', '')
        analysis_type = request.form.get('analysis_type', 'network')
        
        # Parse protein list
        protein_list = [p.strip() for p in proteins_text.split(',') if p.strip()]
        
        if not protein_list:
            return render_template('analyze.html', error="Please enter at least one protein")
        
        try:
            # Determine ID type (assume all are the same type)
            id_type = 'gene' if '_' not in protein_list[0] else 'uniprot'
            
            # Convert to UniProt IDs if needed
            if id_type == 'gene':
                uniprot_ids = []
                for gene in protein_list:
                    try:
                        protein = pe.data.get_protein_by_id(gene_symbol=gene)
                        uniprot_ids.append(protein['uniprot_id'])
                    except Exception as e:
                        logger.warning(f"Could not convert gene {gene} to UniProt ID: {e}")
                protein_list = uniprot_ids
            
            # Perform the analysis
            if analysis_type == 'network':
                # Build network
                network = pe.navigation.build_interaction_network(protein_list, max_depth=1)
                
                # Create network visualization
                network_html = pe.visualization.visualize_network(
                    network,
                    highlight_nodes=protein_list,
                    title=f"Protein Interaction Network"
                )
                
                # Find common interactors
                common_interactors = pe.navigation.find_common_interactors(network, protein_list)
                
                return render_template(
                    'analyze_results.html',
                    analysis_type=analysis_type,
                    proteins=protein_list,
                    network_html=network_html,
                    common_interactors=common_interactors,
                    node_count=network.number_of_nodes(),
                    edge_count=network.number_of_edges()
                )
            elif analysis_type == 'structure':
                # Get structures and compare
                structures = {}
                for uniprot_id in protein_list:
                    try:
                        structure = pe.data.get_alphafold_structure(uniprot_id)
                        if structure:
                            structures[uniprot_id] = structure
                    except Exception as e:
                        logger.warning(f"Could not get structure for {uniprot_id}: {e}")
                
                if not structures:
                    return render_template('analyze.html', 
                                          error="Could not find structures for any of the proteins")
                
                # Create structure visualization
                structure_html = None
                if len(structures) == 1:
                    # Single structure
                    uniprot_id, structure = next(iter(structures.items()))
                    structure_html = pe.visualization.visualize_structure(structure)
                else:
                    # Compare structures
                    structure_html = pe.visualization.compare_structures(list(structures.values()))
                
                return render_template(
                    'analyze_results.html',
                    analysis_type=analysis_type,
                    proteins=protein_list,
                    structure_html=structure_html,
                    structures_found=list(structures.keys())
                )
            else:
                return render_template('analyze.html', error=f"Unknown analysis type: {analysis_type}")
                
        except Exception as e:
            logger.error(f"Error in analysis: {e}")
            return render_template('analyze.html', error=str(e))
    
    # GET request
    return render_template('analyze.html')

@app.route('/api/phosphosites/<uniprot_id>', methods=['GET'])
def api_phosphosites(uniprot_id):
    """API endpoint for phosphorylation site analysis with supplementary data."""
    try:
        # Get protein data
        protein_data = pe.data.get_protein_by_id(uniprot_id=uniprot_id)
        
        # Get sequence
        sequence = protein_data.get('metadata', {}).get('sequence', {}).get('value')
        if not sequence:
            return jsonify({'error': 'Protein sequence not found'}), 404
            
        # Get structure
        structure = pe.data.get_alphafold_structure(uniprot_id)
        if not structure:
            return jsonify({'error': 'Protein structure not found'}), 404
            
        # Get phosphosites with enhanced data
        try:
            # Use phospho.py for site detection
            phosphosites = pe.analysis.phospho.analyze_phosphosites(sequence, structure, uniprot_id)
            
            # Get site IDs for batch query
            site_ids = [f"{uniprot_id}_{site['resno']}" for site in phosphosites]
            
            # Batch query for supplementary data
            supp_data_dict = get_phosphosites_batch(site_ids)
            
            # Add supplementary data to phosphosites
            for site in phosphosites:
                site_id = f"{uniprot_id}_{site['resno']}"
                supp_data = supp_data_dict.get(site_id)
                
                if supp_data:
                    # Add database fields to site data
                    for key, value in supp_data.items():
                        # Skip None values and fields already present
                        if value is not None and key not in site:
                            # Map database field names to our field names
                            if key == 'SITE_+/-7_AA':
                                site['motif'] = value
                            elif key == 'site_plddt':
                                site['mean_plddt'] = value
                            else:
                                site[key] = value
        except Exception as e:
            # Fall back to basic analysis without supplementary data
            logger.warning(f"Enhanced phosphosite analysis failed, using basic analysis: {e}")
            phosphosites = pe.analysis.phospho.analyze_phosphosites(sequence, structure, uniprot_id)
        
        # Get structural matches for each site
        try:
            # Collect site information for batch query
            sites_info = []
            for site in phosphosites:
                sites_info.append({
                    'site': site['site'],
                    'resno': site['resno'],
                    'siteType': site['site'][0] if 'site' in site else ''
                })
            
            # Batch retrieve structural matches
            site_ids = [f"{uniprot_id}_{site['resno']}" for site in phosphosites]
            print("SITE IDS 1")
            print(site_ids)
            all_matches = find_structural_matches_batch(site_ids, rmsd_threshold=5.0)
            print("ALL MATCHES 1")
            print(all_matches)
            # Add matches to each site
            for site in phosphosites:
                site_id = f"{uniprot_id}_{site['resno']}"
                if site_id in all_matches:
                    # Enhance matches with supplementary data
                    raw_matches = all_matches[site_id]
                    enhanced_matches = enhance_structural_matches(raw_matches, site['site'])
                    site['structural_matches'] = enhanced_matches
        except Exception as e:
            logger.warning(f"Error finding structural matches: {e}")
            # Continue without matches
        
        return jsonify(phosphosites)
    except Exception as e:
        logger.error(f"Error in API endpoint: {e}")
        return jsonify({'error': str(e)}), 400
    
@app.route('/api/site_kinase_analysis/<site_name>')
def site_kinase_analysis(site_name):
    """API endpoint for site-specific kinase analysis with enhanced debugging."""
    try:
        # Get protein ID from session or query parameter
        uniprot_id = request.args.get('uniprot_id') or session.get('current_protein')
        
        if not uniprot_id:
            logger.error("No protein ID specified in request or session")
            return jsonify({'error': 'No protein specified'}), 400
            
        # Generate site ID
        residue_number = site_name.lstrip('STY')
        site_id = f"{uniprot_id}_{residue_number}"
        
        logger.info(f"Processing site_kinase_analysis for {site_id}")
        
        # Get matches from cache instead of session
        cache_data = MATCH_CACHE.get(uniprot_id, {})
        if not cache_data:
            logger.error(f"No cache data found for protein {uniprot_id}")
            return jsonify({'error': f'No data found for protein {uniprot_id}'}), 404
            
        structural_matches = cache_data.get('structural_matches', {}).get(site_name, [])
        sequence_matches = cache_data.get('sequence_matches', {}).get(site_name, [])
        
        logger.info(f"Found {len(structural_matches)} structural matches and {len(sequence_matches)} sequence matches")
        
        # Debug the first few matches to verify structure
        if structural_matches:
            logger.debug(f"First structural match sample: {structural_matches[0]}")
        if sequence_matches:
            logger.debug(f"First sequence match sample: {sequence_matches[0]}")
        
        # Import necessary functions from Cantley_Kinase_Scorer
        from protein_explorer.analysis.Cantley_Kinase_Scorer import (
            process_site_specific_matches,
            generate_site_analysis_visualizations
        )
        
        # Process matches for this site
        logger.info("Processing matches for kinase scoring")
        site_data = process_site_specific_matches(site_name, structural_matches, sequence_matches)
        
        # Log site data structure to verify it's correctly formed
        logger.info(f"Site data keys: {site_data.keys()}")
        logger.info(f"Structural processed: {site_data.get('structural', {}).get('processed', 0)}")
        logger.info(f"Sequence processed: {site_data.get('sequence', {}).get('processed', 0)}")
        
        # Generate visualizations
        logger.info("Generating site analysis visualizations")
        html = generate_site_analysis_visualizations(site_data, site_id, uniprot_id)
        
        # Log HTML length as a sanity check
        logger.info(f"Generated HTML length: {len(html)}")
        
        # Add debugging info to the response
        debug_info = {
            'structural_matches_count': len(structural_matches),
            'sequence_matches_count': len(sequence_matches),
            'structural_processed': site_data.get('structural', {}).get('processed', 0),
            'sequence_processed': site_data.get('sequence', {}).get('processed', 0),
            'html_length': len(html)
        }
        
        # Return data and HTML
        return jsonify({
            'site_name': site_name,
            'site_id': site_id,
            'uniprot_id': uniprot_id,
            'structural_matches': len(structural_matches),
            'sequence_matches': len(sequence_matches),
            'debug': debug_info,
            'html': html
        })
    except Exception as e:
        logger.error(f"Error in site kinase analysis: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return jsonify({
            'error': str(e), 
            'traceback': traceback.format_exc()
        }), 500

# Cache maintenance function to prevent memory issues
def cleanup_match_cache(max_age=3600):  # Default: remove entries older than 1 hour
    """Clean up old entries from the match cache to prevent memory issues."""
    current_time = time.time()
    keys_to_remove = []
    
    for key, data in MATCH_CACHE.items():
        if current_time - data.get('timestamp', 0) > max_age:
            keys_to_remove.append(key)
    
    for key in keys_to_remove:
        del MATCH_CACHE[key]
        
    logger.info(f"Cleaned up {len(keys_to_remove)} old entries from match cache")

@app.route('/api/structure/<uniprot_id>', methods=['GET'])
def api_structure(uniprot_id):
    """API endpoint for protein structure."""
    try:
        structure = pe.data.get_alphafold_structure(uniprot_id)
        if not structure:
            return jsonify({'error': 'Structure not found'}), 404
        
        # Return basic info about the structure
        structure_info = pe.data.parse_pdb_structure(structure)
        return jsonify(structure_info)
    except Exception as e:
        return jsonify({'error': str(e)}), 400

@app.route('/api/network/<uniprot_id>', methods=['GET'])
def api_network(uniprot_id):
    """API endpoint for protein network."""
    try:
        depth = int(request.args.get('depth', 1))
        network = pe.navigation.build_interaction_network([uniprot_id], max_depth=depth)
        
        # Convert network to dict for JSON serialization
        network_data = {
            'nodes': list(network.nodes()),
            'edges': list(network.edges()),
            'node_count': network.number_of_nodes(),
            'edge_count': network.number_of_edges()
        }
        
        return jsonify(network_data)
    except Exception as e:
        return jsonify({'error': str(e)}), 400

    

@app.route('/phosphosite', methods=['GET', 'POST'])
def phosphosite_analysis():
    """Phosphosite structural analysis page leveraging improved database structure."""
    
    # Import enhanced table function
    from protein_explorer.analysis.enhanced_table import enhance_phosphosite_table
    from protein_explorer.analysis.motif_blosum_scorer import find_blosum_matches_batch_optimized
    
    # Import Cantley Kinase Scorer functions
    from protein_explorer.analysis.Cantley_Kinase_Scorer import (
        process_matches_for_kinase_scoring,
        create_protein_kinase_report,
        get_html_protein_kinase_report,
        create_site_selector_ui,
        process_site_specific_matches,
        generate_site_analysis_visualizations
    )
    
    # Initialize variables
    results = None
    error = None
    
    if request.method == 'POST':
        # Get identifier from form
        identifier = request.form.get('identifier', '')
        id_type = request.form.get('id_type', 'uniprot')
        
        if not identifier:
            return render_template('phosphosite.html', error="Please enter an identifier")
        
        try:
            # Get protein info using scaffold
            if id_type.lower() == 'uniprot':
                protein_info = {'uniprot_id': identifier}
                protein_data = pe.data.get_protein_by_id(uniprot_id=identifier)
            else:
                protein_data = pe.data.get_protein_by_id(gene_symbol=identifier)
                
            # Extract protein info
            uniprot_id = protein_data.get('uniprot_id')
            gene_symbol = protein_data.get('gene_symbol', 'Unknown')
            name = protein_data.get('metadata', {}).get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown Protein')
            
            protein_info = {
                'uniprot_id': uniprot_id,
                'gene_symbol': gene_symbol,
                'name': name,
                'full_data': protein_data
            }
            
            # Initialize results structure
            results = {
                'protein_info': protein_info,
                'phosphosites': [],
                'structural_matches': None,
                'sequence_matches': None,
                'error': None
            }
            try:
                # STEP 1: Directly query STY_Structural_Annotations for this protein's phosphosites
                from protein_explorer.db.db import find_structural_annotations_by_uniprot
                
                # Get all phosphosite annotations for this protein
                structural_annotations = find_structural_annotations_by_uniprot(uniprot_id)
                
                if not structural_annotations:
                    logger.warning(f"No structural annotations found for {uniprot_id}")
                    # Fall back to sequence info if available
                    sequence = protein_data.get('metadata', {}).get('sequence', {}).get('value')
                    if sequence:
                        # Get structure for fallback analysis if needed
                        structure = pe.data.get_alphafold_structure(uniprot_id)
                        if structure:
                            # Use traditional analysis as fallback
                            phosphosites = pe.analysis.phospho.analyze_phosphosites(sequence, structure, uniprot_id)
                            results['phosphosites'] = phosphosites
                            results['error'] = "Using basic phosphosite analysis (structural annotations not available)"
                else:
                    logger.info(f"Found {len(structural_annotations)} structural annotations for {uniprot_id}")
                    
                    # STEP 2: Convert structural annotations to phosphosite format
                    phosphosites = []
                    site_ids = []
                    
                    # First, create all phosphosite objects and collect site IDs
                    for annotation in structural_annotations:
                        # Map fields from STY_Structural_Annotations to the expected format
                        phosphosite = {
                            'site': annotation.get('Site', ''),
                            'resno': int(annotation.get('ResidueNumber', 0)),
                            'siteType': annotation.get('ResidueType', ''),
                            'motif': annotation.get('Motif', ''),
                            'mean_plddt': annotation.get('pLDDT', 0),
                            'site_plddt': annotation.get('pLDDT', 0),
                            'meanPLDDT': annotation.get('pLDDT', 0),
                            'nearby_count': annotation.get('NeighborCount', 0),
                            'nearbyCount': annotation.get('NeighborCount', 0),
                            'surface_accessibility': annotation.get('HydroxylExposure', 0) * 100,
                            'surfaceAccessibility': annotation.get('HydroxylExposure', 0) * 100,
                            'acidic_aa_percent': annotation.get('negative_count', 0) / max(1, annotation.get('NeighborCount', 1)) * 100,
                            'acidicPercentage': annotation.get('negative_count', 0) / max(1, annotation.get('NeighborCount', 1)) * 100,
                            'basic_aa_percent': annotation.get('positive_count', 0) / max(1, annotation.get('NeighborCount', 1)) * 100,
                            'basicPercentage': annotation.get('positive_count', 0) / max(1, annotation.get('NeighborCount', 1)) * 100,
                            'hydrophobicityScore': annotation.get('hydrophobic_ratio', 0.5) * 100
                        }
                        
                        # Construct the PhosphositeID properly
                        site_id = f"{uniprot_id}_{annotation.get('ResidueNumber', '')}"
                        phosphosite['site_id'] = site_id  # Store for later reference
                        
                        phosphosites.append(phosphosite)
                        site_ids.append(site_id)
                    
                    # STEP 3: Now do a batch query for all phosphosite supplementary data
                    from protein_explorer.db.db import get_phosphosites_batch
                    supp_data_dict = get_phosphosites_batch(site_ids)
                    logger.info(f"Retrieved supplementary data for {len(supp_data_dict)} out of {len(site_ids)} phosphosites")
                    
                    # Add supplementary data to each phosphosite
                    for phosphosite in phosphosites:
                        site_id = phosphosite['site_id']
                        supp_data = supp_data_dict.get(site_id)
                        
                        if supp_data:
                            phosphosite['is_known'] = supp_data.get('is_known_phosphosite', False)
                            phosphosite['is_known_phosphosite'] = supp_data.get('is_known_phosphosite', 0)
                            
                            # Add kinase information if available
                            known_kinases = []
                            for i in range(1, 6):
                                kinase_field = f"KINASE_{i}"
                                if kinase_field in supp_data and supp_data[kinase_field] and supp_data[kinase_field] != 'unlabeled':
                                    known_kinases.append(supp_data[kinase_field])
                            
                            if known_kinases:
                                phosphosite['known_kinase'] = ', '.join(known_kinases)
                        else:
                            phosphosite['is_known'] = False
                            phosphosite['is_known_phosphosite'] = 0
                    
                    results['phosphosites'] = phosphosites
                
                # STEP 4: Find structural matches for phosphosites using site_ids
                if results['phosphosites']:
                    # Create site_ids for batch query if not already created
                    if not site_ids:
                        site_ids = [f"{uniprot_id}_{site['resno']}" for site in results['phosphosites']]
                    
                    # Use the batch function that leverages the database
                    from protein_explorer.db.db import find_structural_matches_batch
                    all_matches = find_structural_matches_batch(site_ids, rmsd_threshold=10.0)
                    # Create a dictionary of target IDs for batch lookup
                    target_ids = set()
                    for matches in all_matches.values():
                        for match in matches:
                            target_id = match.get('target_id')
                            if target_id:
                                target_ids.add(target_id)
                    
                    # Get structural annotations for all target sites in batch
                    from protein_explorer.db.db import get_structural_annotations_batch
                    target_annotations = {}
                    if target_ids:
                        # First get the raw annotations
                        raw_annotations = get_structural_annotations_batch(list(target_ids))
                        logger.info(f"Retrieved {len(raw_annotations)} structural annotations for target sites")
                        
                        # Create a new dictionary with keys from the Site field
                        transformed_annotations = {}
                        for key, annotation in raw_annotations.items():
                            # Keep the original key to ensure we don't lose data
                            transformed_annotations[key] = annotation
                            # Also store using the Site value as the key if available
                            if 'Site' in annotation:
                                site_value = annotation['Site']
                                transformed_annotations[site_value] = annotation
                        
                        # Replace the original dictionary with the transformed one
                        target_annotations = transformed_annotations
                        
                        logger.info(f"Transformed annotations dictionary with {len(target_annotations)} entries")

                    # Get supplementary data for target sites in batch
                    target_supp_data = get_phosphosites_batch(list(target_ids))
                    logger.info(f"Retrieved supplementary data for {len(target_supp_data)} target sites")
                    
                    # Now enhance the matches with complete structural annotation data
                    structural_matches = {}
                    for site in results['phosphosites']:
                        site_id = f"{uniprot_id}_{site['resno']}"
                        if site_id in all_matches:
                            matches = all_matches[site_id]
                            if matches:
                                enhanced_matches = []
                                for match in matches:
                                    # Create enhanced copy
                                    enhanced_match = match.copy()
                                    
                                    # Get target structural annotation data
                                    target_id = match.get('target_id')
                                    if target_id and target_id in target_annotations:
                                        annotation = target_annotations[target_id]
                                        
                                        # Add ALL fields from structural annotations
                                        for key, value in annotation.items():
                                            # Skip keys that would overwrite existing critical data
                                            if key not in ['target_id', 'target_uniprot', 'target_site']:
                                                enhanced_match[key] = value
                                        
                                        # Ensure motif is uppercase if present
                                        if 'Motif' in enhanced_match and enhanced_match['Motif']:
                                            enhanced_match['motif'] = str(enhanced_match['Motif']).upper()
                                    
                                    # Add ALL fields from supplementary data
                                    if target_id and target_id in target_supp_data:
                                        supp = target_supp_data[target_id]
                                        
                                        for key, value in supp.items():
                                            # Skip keys that would overwrite existing critical data
                                            if key not in ['target_id', 'target_uniprot', 'target_site'] and value is not None:
                                                enhanced_match[key] = value
                                    
                                    enhanced_matches.append(enhanced_match)
                                
                                structural_matches[site['site']] = enhanced_matches

                    results['structural_matches'] = structural_matches

                # Add these two blocks of code in the appropriate places

                # Add after STEP 4 (where structural matches are processed)
                # Create phosphosites_structure object by including target site information
                if results['phosphosites'] and results['structural_matches']:
                    # Create an array to hold all phosphosites (original protein sites + matching sites)
                    phosphosites_structure = []
                    
                    # First, add all original protein phosphosites
                    for site in results['phosphosites']:
                        phosphosites_structure.append(site.copy())
                    
                    # Now add all target sites from structural matches
                    for site_name, matches in results['structural_matches'].items():
                        for match in matches:
                            # Only add if target_id exists
                            if 'target_id' in match:
                                target_id = match['target_id']
                                
                                # Create a new phosphosite object with basic information
                                target_site = {
                                    'site': match.get('target_site', ''),
                                    'site_id': target_id,
                                    'uniprot': match.get('target_uniprot', ''),
                                    'resno': int(target_id.split('_')[-1]) if '_' in target_id else 0,
                                    'is_target': True,  # Flag to identify this as a target site
                                    'source_site': site_name,  # Reference to the source site
                                    'rmsd': match.get('rmsd', 0),
                                    'match_type': 'structural'
                                }
                                
                                # Add ALL fields from the match
                                for key, value in match.items():
                                    # Skip keys we've already processed or that would overwrite critical data
                                    if key not in ['target_id', 'target_uniprot', 'target_site']:
                                        target_site[key] = value
                                
                                # Add to the phosphosites_structure list
                                phosphosites_structure.append(target_site)
                    
                    # Sort by uniprot and then by residue number
                    phosphosites_structure.sort(key=lambda x: (x.get('uniprot', ''), x.get('resno', 0)))
                    
                    # Add to results
                    results['phosphosites_structure'] = phosphosites_structure
                    logger.info(f"Created phosphosites_structure with {len(phosphosites_structure)} entries")
                else:
                    # If no matches or no phosphosites, just use original phosphosites
                    results['phosphosites_structure'] = [site.copy() for site in results.get('phosphosites', [])]
                
                # STEP 5: Find sequence similarity matches for phosphosites
                if results['phosphosites']:
                    # Use site_ids created above for batch queries
                    from protein_explorer.db.db import find_sequence_matches_batch
                    from protein_explorer.analysis.sequence_analyzer_2 import enhance_sequence_matches_batch
                    
                    # Get matches from optimized database query
                    sequence_matches_batch = find_sequence_matches_batch(site_ids, min_similarity=0.3)
                    
                    if sequence_matches_batch:
                        # STEP 1: Get matching phosphosites from find_sequence_matches_batch function
                        logger.info(f"Found sequence matches for {len(sequence_matches_batch)} sites")
                        
                        # Create a set of all target IDs for batch processing
                        target_ids = set()
                        for matches in sequence_matches_batch.values():
                            for match in matches:
                                if 'target_id' in match:
                                    target_ids.add(match['target_id'])
                        
                        # STEP 2: Use get_structural_annotations_batch to enhance the matches
                        from protein_explorer.db.db import get_structural_annotations_batch
                        structural_annotations = {}
                        
                        if target_ids:
                            try:
                                # Get structural annotations for the target IDs
                                raw_annotations = get_structural_annotations_batch(list(target_ids))
                                logger.info(f"Retrieved {len(raw_annotations)} structural annotations")
                                
                                # Create a mapping from target_id to annotation data
                                # Match the 'Site' field with our target_id format
                                for annotation_id, annotation_data in raw_annotations.items():
                                    # Store using original key
                                    structural_annotations[annotation_id] = annotation_data
                                    
                                    # Also store with Site field as key if available
                                    if 'Site' in annotation_data:
                                        site_value = annotation_data['Site']
                                        # Add this entry to our lookup dictionary using the Site value as key
                                        structural_annotations[site_value] = annotation_data
                                
                                logger.info(f"Created mapping for {len(structural_annotations)} structural annotations")
                                
                                # Debug - log the first few keys to verify format
                                sample_keys = list(structural_annotations.keys())[:5]
                                logger.debug(f"Sample structural annotation keys: {sample_keys}")
                                
                            except Exception as e:
                                logger.error(f"Error retrieving structural annotations: {e}")
                                import traceback
                                logger.error(traceback.format_exc())
                        
                        # STEP 3: Add the augmented information to match data
                        enhanced_matches = {}
                        
                        for site_id, matches in sequence_matches_batch.items():
                            enhanced_site_matches = []
                            
                            for match in matches:
                                # Skip matches with missing target_id
                                if 'target_id' not in match:
                                    continue
                                
                                # Skip Tyrosine matches (filter them out here)
                                target_site = match.get('target_site', '')
                                if target_site and target_site[0] == 'Y':
                                    continue
                                
                                target_id = match['target_id']
                                
                                # Create enhanced match (copy original)
                                enhanced_match = match.copy()
                                
                                # Look up the structural annotation directly using target_id
                                if target_id in structural_annotations:
                                    annotation = structural_annotations[target_id]
                                    logger.debug(f"Found annotation for target_id {target_id}")
                                    
                                    # Add structural details
                                    enhanced_match['site_plddt'] = annotation.get('pLDDT', enhanced_match.get('site_plddt', 0))
                                    enhanced_match['mean_plddt'] = annotation.get('MotifPLDDT', enhanced_match.get('mean_plddt', 0))
                                    enhanced_match['nearby_count'] = annotation.get('NeighborCount', enhanced_match.get('nearby_count', 0))
                                    enhanced_match['surface_accessibility'] = annotation.get('HydroxylExposure', 0) * 100  # Convert to percentage
                                    
                                    # Add motif if available
                                    if 'Motif' in annotation and annotation['Motif']:
                                        enhanced_match['motif'] = annotation['Motif'].upper()
                                    
                                    # Add secondary structure and other metrics
                                    enhanced_match['secondary_structure'] = annotation.get('SecondaryStructure', '')
                                    enhanced_match['hydrophobic_ratio'] = annotation.get('hydrophobic_ratio', 0)
                                    enhanced_match['net_charge'] = annotation.get('net_charge', 0)
                                    
                                    # Add all remaining structural annotation data
                                    for key, value in annotation.items():
                                        # Skip keys we've already processed or that would overwrite critical data
                                        if key not in ['pLDDT', 'MotifPLDDT', 'NeighborCount', 'HydroxylExposure', 
                                                    'Motif', 'SecondaryStructure', 'hydrophobic_ratio', 'net_charge',
                                                    'PhosphositeID', 'Site'] and key not in enhanced_match:
                                            enhanced_match[key] = value
                                
                                enhanced_site_matches.append(enhanced_match)
                            
                            # Only add if we have matches after filtering
                            if enhanced_site_matches:
                                # Get the site name to use as the key
                                site_name = None
                                site_resno = None
                                
                                # Extract site name from site_id
                                if '_' in site_id:
                                    uniprot, site_resno = site_id.split('_', 1)
                                    
                                # Find the site in phosphosites data to get site name
                                if site_resno:
                                    for site in results['phosphosites']:
                                        if str(site.get('resno')) == site_resno:
                                            site_name = site.get('site')
                                            break
                                
                                # If site name not found, use residue number as fallback
                                if not site_name and site_resno:
                                    # Try to determine site type
                                    site_type = None
                                    for phosphosite in results['phosphosites']:
                                        if phosphosite.get('resno') == int(site_resno):
                                            site_type = phosphosite.get('siteType') or phosphosite.get('site', [''])[0]
                                            break
                                    
                                    site_name = f"{site_type or 'S'}{site_resno}"
                                
                                # Add to results dictionary
                                if site_name:
                                    enhanced_matches[site_name] = enhanced_site_matches
                        
                        # STEP 4: Get kinase and known phosphorylation information
                        if target_ids:
                            from protein_explorer.db.db import get_phosphosites_batch
                            supp_data_batch = get_phosphosites_batch(list(target_ids))
                            
                            if supp_data_batch:
                                logger.info(f"Retrieved supplementary data for {len(supp_data_batch)} target sites")
                                
                                # Enhance matches with kinase and known site information
                                for site_name, matches in enhanced_matches.items():
                                    for match in matches:
                                        target_id = match.get('target_id')
                                        if target_id and target_id in supp_data_batch:
                                            supp_data = supp_data_batch[target_id]
                                            
                                            # Add is_known status
                                            match['is_known_phosphosite'] = supp_data.get('is_known_phosphosite', 0)
                                            match['is_known'] = bool(match['is_known_phosphosite'])
                                            
                                            # Add kinase information
                                            known_kinases = []
                                            for i in range(1, 6):
                                                kinase_field = f"KINASE_{i}"
                                                if kinase_field in supp_data and supp_data[kinase_field] and supp_data[kinase_field] != 'unlabeled':
                                                    known_kinases.append(supp_data[kinase_field])
                                            
                                            if known_kinases:
                                                match['known_kinase'] = ', '.join(known_kinases)
                                            
                                            # Add all other supplementary data
                                            for key, value in supp_data.items():
                                                # Skip keys we've already processed or that might overwrite critical data
                                                if key not in ['is_known_phosphosite', 'KINASE_1', 'KINASE_2', 'KINASE_3', 'KINASE_4', 'KINASE_5'] and key not in match:
                                                    match[key] = value
                        
                        # STEP 5: Pass this information to network visualization
                        results['sequence_matches'] = enhanced_matches
                        logger.info(f"Prepared sequence network data with {len(enhanced_matches)} sites and all structural annotations")
                    else:
                        logger.info("No sequence matches found for any site")
                        results['sequence_matches'] = {}
                else:
                    results['sequence_matches'] = {}

                if results['phosphosites'] and results['sequence_matches']:
                    # Create an array to hold all phosphosites (original protein sites + matching sites)
                    phosphosites_sequence = []
                    
                    # First, add all original protein phosphosites
                    for site in results['phosphosites']:
                        phosphosites_sequence.append(site.copy())
                    
                    # Now add all target sites from sequence matches
                    for site_name, matches in results['sequence_matches'].items():
                        for match in matches:
                            # Only add if target_id exists
                            if 'target_id' in match:
                                target_id = match['target_id']
                                
                                # Skip if it's a Y site
                                target_site_name = match.get('target_site', '')
                                if target_site_name and target_site_name[0] == 'Y':
                                    continue
                                
                                # Create a new phosphosite object for this target
                                target_site = {
                                    'site': match.get('target_site', ''),
                                    'site_id': target_id,
                                    'uniprot': match.get('target_uniprot', ''),
                                    'resno': int(target_id.split('_')[-1]) if '_' in target_id else 0,
                                    'siteType': match.get('site_type', target_site_name[0] if target_site_name else ''),
                                    'motif': match.get('motif', ''),
                                    'mean_plddt': match.get('mean_plddt', 0),
                                    'site_plddt': match.get('site_plddt', 0),
                                    'meanPLDDT': match.get('mean_plddt', 0),
                                    'nearby_count': match.get('nearby_count', 0),
                                    'nearbyCount': match.get('nearby_count', 0),
                                    'surface_accessibility': match.get('surface_accessibility', 0),
                                    'surfaceAccessibility': match.get('surface_accessibility', 0),
                                    'is_known': match.get('is_known', False),
                                    'is_known_phosphosite': match.get('is_known_phosphosite', 0),
                                    'is_target': True,  # Flag to identify this as a target site
                                    'source_site': site_name,  # Reference to the source site
                                    'similarity': match.get('similarity', 0),
                                    'match_type': 'sequence'
                                }
                                
                                # Add kinase information if available
                                if 'known_kinase' in match:
                                    target_site['known_kinase'] = match['known_kinase']
                                    target_site['known_kinases'] = match['known_kinase']
                                
                                # Add secondary structure and other structural information
                                if 'secondary_structure' in match:
                                    target_site['secondary_structure'] = match['secondary_structure']
                                    
                                if 'hydrophobic_ratio' in match:
                                    target_site['hydrophobic_ratio'] = match['hydrophobic_ratio']
                                    target_site['hydrophobicityScore'] = match['hydrophobic_ratio'] * 100
                                
                                # Add any other fields that are available in the match
                                for key, value in match.items():
                                    if key not in target_site and key not in ['target_id', 'target_uniprot', 'target_site']:
                                        target_site[key] = value
                                
                                # Add to the phosphosites_sequence list
                                phosphosites_sequence.append(target_site)
                    
                    # Sort by uniprot and then by residue number
                    phosphosites_sequence.sort(key=lambda x: (x.get('uniprot', ''), x.get('resno', 0)))
                    
                    # Add to results
                    results['phosphosites_sequence'] = phosphosites_sequence
                    logger.info(f"Created phosphosites_sequence with {len(phosphosites_sequence)} entries")
                else:
                    # If no matches or no phosphosites, just use original phosphosites
                    results['phosphosites_sequence'] = [site.copy() for site in results.get('phosphosites', [])]
                
            except Exception as e:
                logger.error(f"Error analyzing phosphosites: {e}")
                import traceback
                logger.error(traceback.format_exc())
                results['error'] = f"Error analyzing phosphosites: {str(e)}"
            
            # Use enhanced table function for display
            if results['phosphosites']:
                # Create enhanced HTML for the phosphosites
                phosphosites_html = enhance_phosphosite_table(
                    results['phosphosites'], 
                    results['protein_info']['uniprot_id']
                )
                
                # Add the HTML to the results
                results['phosphosites_html'] = phosphosites_html

            fields_to_remove = [
                'DISEASE', 'ALTERATION', 'Disease_DOMAIN', 'DISEASE_PMIDs', 'DISEASE_NOTES', 
                'Regulatory_Domain', 'ON_FUNCTION', 'ON_PROCESS', 'ON_PROT_INTERACT', 
                'ON_OTHER_INTERACT', 'REGULATORY_PMIDs', 'REGULATORY_NOTES', 'SUB_GENE_ID'
            ]

            # Helper function to remove fields from an object
            def remove_unwanted_fields(obj):
                for field in fields_to_remove:
                    if field in obj:
                        del obj[field]

            # Remove unwanted fields from structural_matches dictionary
            if results.get('structural_matches'):
                for site_key, matches in results['structural_matches'].items():
                    for match in matches:
                        remove_unwanted_fields(match)

            # Remove unwanted fields from sequence_matches dictionary
            if results.get('sequence_matches'):
                for site_key, matches in results['sequence_matches'].items():
                    for match in matches:
                        remove_unwanted_fields(match)

            # Remove unwanted fields from phosphosites_structure
            if results.get('phosphosites_structure'):
                for site in results['phosphosites_structure']:
                    remove_unwanted_fields(site)

            # Remove unwanted fields from phosphosites_sequence
            if results.get('phosphosites_sequence'):
                for site in results['phosphosites_sequence']:
                    remove_unwanted_fields(site)

            # Also remove from original phosphosites list just to be thorough
            if results.get('phosphosites'):
                for site in results['phosphosites']:
                    remove_unwanted_fields(site)

            # Log the removal for debugging purposes
            logger.info(f"Removed unwanted fields ({', '.join(fields_to_remove)}) from all result objects")

            # Generate structural network visualization
            network_visualization = create_phosphosite_network_visualization(
                results['protein_info']['uniprot_id'],
                results['phosphosites_structure'],
                results['structural_matches']
            )
            
            # Generate sequence similarity network visualization
            from protein_explorer.visualization.protein_sequence_phosphosite_network import create_sequence_network_visualization
            
            sequence_network_visualization = create_sequence_network_visualization(
                results['protein_info']['uniprot_id'],
                results['phosphosites_sequence'],
                results['sequence_matches']
            )
            
            from itertools import islice
            def pretty_print_first_kvp(dictionary, label, limit=5):
                print(f"\n{label}")
                for k, v in islice(dictionary.items(), limit):
                    print(f"Key: {k}")
                    if isinstance(v, list):
                        print("Value (first item in list):")
                        print(json.dumps(v[0], indent=2))  # only show first item for brevity
                    else:
                        print("Value:")
                        print(json.dumps(v, indent=2))
                    print("-" * 40)

            # Usage
            pretty_print_first_kvp(results['structural_matches'], "STRUCTURAL_MATCH DATA")
            pretty_print_first_kvp(results['sequence_matches'], "SEQUENCE_MATCH DATA")
            print("PHOSPHOSITES DATA")
            print(results['phosphosites'][:5])

            # Generate Cantley Kinase Score visualizations
            if results['phosphosites']:
                # Process structural matches for kinase scoring
                structural_scores = {}
                all_structural_processed = []
                
                if results.get('structural_matches'):
                    for site_name, matches in results['structural_matches'].items():
                        # Process matches for this site
                        scores, processed = process_matches_for_kinase_scoring(
                            matches, match_type='structural'
                        )
                        structural_scores[site_name] = scores
                        all_structural_processed.extend(processed)
                        
                        # Log results
                        logger.info(f"Processed {len(processed)} structural matches for {site_name} with kinase data")
                
                # Process sequence matches for kinase scoring
                sequence_scores = {}
                all_sequence_processed = []
                
                if results.get('sequence_matches'):
                    for site_name, matches in results['sequence_matches'].items():
                        # Process matches for this site
                        scores, processed = process_matches_for_kinase_scoring(
                            matches, match_type='sequence'
                        )
                        sequence_scores[site_name] = scores
                        all_sequence_processed.extend(processed)
                        
                        # Log results
                        logger.info(f"Processed {len(processed)} sequence matches for {site_name} with kinase data")
                
                # Generate kinase report
                kinase_report = create_protein_kinase_report(
                    results['protein_info']['uniprot_id'],
                    results['phosphosites'],
                    results['structural_matches'],
                    results['sequence_matches']
                )
                
                # Generate HTML for kinase report
                results['kinase_html'] = get_html_protein_kinase_report(kinase_report)
                
                # Create site selector UI
                site_selector_html = create_site_selector_ui(
                    results['phosphosites'],
                    results['structural_matches'],
                    results['sequence_matches'],
                    active_site=None  # No site selected initially
                )
                
                # Store in results
                results['site_selector_html'] = site_selector_html
                
                # Store in session for API access
                # With these lines:
                cache_key = results['protein_info']['uniprot_id']
                MATCH_CACHE[cache_key] = {
                    'structural_matches': results['structural_matches'],
                    'sequence_matches': results['sequence_matches'],
                    'timestamp': time.time()  # Add timestamp for potential expiration
                }

                # Store just the protein ID in the session
                session['current_protein'] = results['protein_info']['uniprot_id']
                
                # Initialize empty site analysis container
                results['site_analysis_html'] = """
                <div id="site-analysis">
                    <div class="alert alert-info">
                        Select a phosphosite from the dropdown above to view kinase analysis.
                    </div>
                </div>
                """
            
            return render_template('phosphosite.html', 
                          protein_info=results['protein_info'],
                          phosphosites=results['phosphosites'],
                          phosphosites_html=results.get('phosphosites_html'),
                          structural_matches=results['structural_matches'],
                          sequence_matches=results.get('sequence_matches'),
                          network_visualization=network_visualization,
                          sequence_network_visualization=sequence_network_visualization,
                          kinase_html=results.get('kinase_html'),
                          site_selector_html=results.get('site_selector_html'),
                          site_analysis_html=results.get('site_analysis_html'),
                          error=results.get('error'))
                
        except Exception as e:
            print(f"Error in phosphosite analysis: {e}")
            import traceback
            print(traceback.format_exc())
            error = str(e)
            return render_template('phosphosite.html', error=error)
    
    # GET request - show empty form
    return render_template('phosphosite.html')



# Add this to app.py
@app.route('/api/sequence_matches/<site_id>', methods=['GET'])
def api_sequence_matches(site_id):
    """API endpoint for sequence similarity matches."""
    try:
        # Get query parameters
        top_n = request.args.get('top_n', default=100, type=int)
        min_similarity = request.args.get('min_similarity', default=0.4, type=float)
        
        # Validate site_id format (e.g., UniProtID_ResidueNumber)
        if '_' not in site_id:
            return jsonify({'error': 'Invalid site ID format. Expected: UniProtID_ResidueNumber'}), 400
            
        # Get sequence matches from database
        matches = find_sequence_matches(site_id, min_similarity=min_similarity)
        
        # Sort by similarity (highest first) and take top N
        if matches:
            matches.sort(key=lambda x: x['similarity'], reverse=True)
            if top_n:
                matches = matches[:top_n]
        
        # If no matches, return empty list with a message
        if not matches:
            return jsonify({
                'matches': [],
                'count': 0,
                'message': f'No sequence similarity matches found for {site_id}'
            })
        
        # Create a complete response with metadata
        response = {
            'site_id': site_id,
            'matches': matches,
            'count': len(matches),
            'params': {
                'top_n': top_n,
                'min_similarity': min_similarity
            }
        }
        
        return jsonify(response)
    except Exception as e:
        logger.error(f"Error in sequence matches API: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/sequence_conservation/<site_id>', methods=['GET'])
def api_sequence_conservation(site_id):
    """API endpoint for sequence conservation analysis."""
    try:
        # Get query parameters
        top_n = request.args.get('top_n', default=200, type=int)
        min_similarity = request.args.get('min_similarity', default=0.4, type=float)
        
        # Get the query motif if possible
        query_motif = None
        supp_data = get_phosphosite_data(site_id)
        if supp_data:
            if 'SITE_+/-7_AA' in supp_data and supp_data['SITE_+/-7_AA']:
                query_motif = supp_data['SITE_+/-7_AA']
            elif 'motif' in supp_data and supp_data['motif']:
                query_motif = supp_data['motif']
        
        # Get sequence matches from database
        matches = find_sequence_matches(site_id, min_similarity=min_similarity)
        
        # Sort by similarity (highest first) and take top N
        if matches:
            matches.sort(key=lambda x: x['similarity'], reverse=True)
            if top_n:
                matches = matches[:top_n]
        
        # Analyze conservation
        conservation = analyze_motif_conservation(matches, query_motif=query_motif)
        
        # Add metadata to response
        response = {
            'site_id': site_id,
            'analysis': conservation,
            'match_count': len(matches),
            'query_motif': query_motif
        }
        
        return jsonify(response)
    except Exception as e:
        logger.error(f"Error in sequence conservation API: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/site/<uniprot_id>/<site>')
def site_detail(uniprot_id, site):
    """Display detailed information about a specific phosphorylation site with enhanced supplementary data."""
    try:
        # Get protein data
        protein_data = pe.data.get_protein_by_id(uniprot_id=uniprot_id)
        
        # Get the sequence from metadata
        sequence = protein_data.get('metadata', {}).get('sequence', {}).get('value')
        if not sequence:
            return render_template('error.html', error=f"Protein sequence not found for {uniprot_id}")
            
        # Get structure
        structure = pe.data.get_alphafold_structure(uniprot_id)
        if not structure:
            return render_template('error.html', error=f"Protein structure not found for {uniprot_id}")
        
        # Parse the site string to get type and residue number
        site_match = re.match(r'([A-Z])(\d+)', site)
        if not site_match:
            return render_template('error.html', error=f"Invalid site format: {site}")
            
        site_type = site_match.group(1)
        site_number = int(site_match.group(2))
        
        # Generate site_id in the format UniProtID_ResidueNumber
        site_id = f"{uniprot_id}_{site_number}"
        
        # Get supplementary data for this site from database
        site_data = get_phosphosite_data(site_id)
        
        # If site_data is None, create a minimal site data dictionary
        if not site_data:
            site_data = {
                'site_id': site_id,
                'site': site,
                'resno': site_number,
                'siteType': site_type
            }
        
        # Analyze phosphosites to get site-specific data
        all_phosphosites = analyze_phosphosites(sequence, structure, uniprot_id)
        
        # Find the specific site
        site_data_from_analysis = next((s for s in all_phosphosites if s['site'] == site), None)
        if site_data_from_analysis:
            # Merge with existing site_data
            for key, value in site_data_from_analysis.items():
                if key not in site_data or site_data[key] is None:
                    site_data[key] = value
        
        # Find structural matches
        structural_matches = []
        try:
            # Get raw matches from database
            struct_matches = find_structural_matches(site_id, rmsd_threshold=5.0)
            
            # Filter out self-matches
            struct_matches = [match for match in struct_matches if match.get('rmsd', 0) > 0.01]
            
            # Enhance matches with supplementary data
            structural_matches = enhance_structural_matches(struct_matches, site)
            
            # Sort matches by RMSD
            structural_matches = sorted(structural_matches, key=lambda x: x.get('rmsd', float('inf')))
            
        except Exception as e:
            logger.error(f"Error finding structural matches: {e}")
            structural_matches = []
        
        # Create basic structure visualization
        structure_html = pe.visualization.visualize_structure(
            structure,
            sequence=sequence
        )
        
        # Create specialized visualizations
        try:
            network_html = create_phosphosite_network(site, structural_matches, site_data)
        except Exception as e:
            logger.error(f"Error creating network visualization: {e}")
            network_html = f"<div class='alert alert-warning'>Error creating network visualization: {str(e)}</div>"
            
        # 2. Comparative motif visualization
        try:
            if structural_matches and ('motif' in site_data or 'SITE_+/-7_AA' in site_data):
                # Use SITE_+/-7_AA if motif is not present
                if 'motif' not in site_data and 'SITE_+/-7_AA' in site_data:
                    site_data['motif'] = site_data['SITE_+/-7_AA']
                    
                motif_html = create_comparative_motif_visualization(site_data, structural_matches)
            else:
                motif_html = "<div class='alert alert-info'>No motif data available for comparison.</div>"
        except Exception as e:
            logger.error(f"Error creating motif visualization: {e}")
            motif_html = f"<div class='alert alert-warning'>Error creating motif visualization: {str(e)}</div>"
            
        # 3. Residue distribution analysis
        try:
            if structural_matches:
                distribution_data = analyze_residue_distributions(structural_matches)
            else:
                distribution_data = None
        except Exception as e:
            logger.error(f"Error analyzing residue distributions: {e}")
            distribution_data = None
            
        # 4. Enhanced 3D visualization
        try:
            enhanced_3d_html = enhance_site_visualization(uniprot_id, site, site_data)
        except Exception as e:
            logger.error(f"Error creating enhanced 3D visualization: {e}")
            enhanced_3d_html = f"<div class='alert alert-warning'>Error creating enhanced 3D visualization: {str(e)}</div>"
            
        # 5. Structural context analysis
        try:
            context_data = analyze_phosphosite_context(structure, site_number, site_type)
        except Exception as e:
            logger.error(f"Error analyzing structural context: {e}")
            context_data = None
        
        # Find sequence similarity matches
        sequence_matches = []
        sequence_network_data = None
        sequence_conservation = None
        sequence_motif_html = None
        
        try:
            # Get sequence similarity matches from database
            sequence_matches = find_sequence_matches(site_id, min_similarity=0.4)
            logger.info(f"Found {len(sequence_matches)} sequence matches")
            
            # Sort by similarity (highest first) and take top 200
            if sequence_matches:
                sequence_matches.sort(key=lambda x: x['similarity'], reverse=True)
                sequence_matches = sequence_matches[:200]
                
                # Log a sample of the matches to see what they contain
                logger.info(f"Sample match: {sequence_matches[0]}")
                
                # Create sequence network data
                query_motif = site_data.get('motif') or site_data.get('SITE_+/-7_AA')
                sequence_network_data = create_sequence_network_data(
                    site_id, 
                    sequence_matches,
                    query_motif=query_motif
                )
                logger.info(f"Created network data with {len(sequence_network_data['nodes'])} nodes")
                
                # Create sequence conservation analysis
                sequence_conservation = analyze_motif_conservation(
                    sequence_matches,
                    query_motif=query_motif
                )
                logger.info(f"Created conservation analysis with {sequence_conservation['motif_count']} motifs")
                
                # Create sequence motif visualization
                sequence_motif_html = create_sequence_motif_visualization(
                    site_id,
                    query_motif,
                    sequence_matches
                )
                logger.info(f"Created motif visualization HTML of length {len(sequence_motif_html) if sequence_motif_html else 0}")
            else:
                logger.warning(f"No sequence matches found for {site_id}")
        except Exception as e:
            logger.error(f"Error analyzing sequence similarity: {e}")
            import traceback
            logger.error(traceback.format_exc())
            sequence_matches = []
        
        # Standard kinase prediction data
        structure_kinase_data = {}
        sequence_kinase_data = {}
        
        # Network-based kinase prediction data
        structure_network_kinase_data = {}
        sequence_network_kinase_data = {}
        
        try:
            # Standard structural kinase prediction
            struct_scores = get_site_kinase_scores(site_id, 'structure')
            if struct_scores and 'scores' in struct_scores:
                # Get top kinases
                top_kinases = predict_kinases(site_id, top_n=10, score_type='structure')
                
                # Get kinase families
                kinase_families = categorize_kinases_by_family(top_kinases)
                
                # Get known kinase info
                known_kinase = get_known_kinase_info(site_id, 'structure')
                
                # Prepare data for heatmap visualization
                struct_match_ids = [f"{match['target_uniprot']}_{match['target_site'].replace('S', '').replace('T', '').replace('Y', '')}" 
                                   for match in structural_matches[:20] if match.get('rmsd', 10) < 5.0]
                struct_match_ids.insert(0, site_id)  # Add the query site first
                
                heatmap_data = get_heatmap_data(struct_match_ids, top_n=10, score_type='structure')
                
                # Bundle all data
                structure_kinase_data = {
                    'site_id': site_id,
                    'top_kinases': top_kinases,
                    'kinase_families': kinase_families,
                    'known_kinase': known_kinase,
                    'heatmap': heatmap_data
                }
            
            # Standard sequence kinase prediction
            seq_scores = get_site_kinase_scores(site_id, 'sequence')
            if seq_scores and 'scores' in seq_scores:
                # Get top kinases
                top_kinases = predict_kinases(site_id, top_n=10, score_type='sequence')
                
                # Get kinase families
                kinase_families = categorize_kinases_by_family(top_kinases)
                
                # Get known kinase info
                known_kinase = get_known_kinase_info(site_id, 'sequence')
                
                # Prepare data for heatmap visualization
                seq_match_ids = []
                if sequence_matches:
                    seq_match_ids = [match['target_id'] for match in sequence_matches[:20] 
                                    if match.get('similarity', 0) > 0.6]
                    seq_match_ids.insert(0, site_id)  # Add the query site first
                
                heatmap_data = get_heatmap_data(seq_match_ids, top_n=10, score_type='sequence')
                
                # Bundle all data
                sequence_kinase_data = {
                    'site_id': site_id,
                    'top_kinases': top_kinases,
                    'kinase_families': kinase_families,
                    'known_kinase': known_kinase,
                    'heatmap': heatmap_data
                }
            
            # Network-based kinase predictions
            try:
                # 1. Structure network prediction
                structure_network_kinases = predict_kinases_network(
                    site_id, top_n=10, score_type='structure',
                    similarity_threshold=0.6, rmsd_threshold=3.0
                )
                
                if structure_network_kinases:
                    # Get similar sites
                    struct_similar_sites = get_similar_sites(
                        site_id, uniprot_id, site_data,
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Get heatmap data with enhanced function
                    struct_heatmap_data = get_network_heatmap_data(
                        site_id, top_n=10, score_type='structure',
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Log heatmap data for debugging
                    logger.info(f"Structure heatmap data: {len(struct_heatmap_data['sites'])} sites, " 
                                f"{len(struct_heatmap_data['kinases'])} kinases, "
                                f"{len(struct_heatmap_data['scores'])} scores")
                    
                    # Get family distribution
                    struct_family_dist = get_kinase_family_distribution_network(
                        site_id, score_type='structure',
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Combine into data structure
                    structure_network_kinase_data = {
                        'site_id': site_id,
                        'top_kinases': structure_network_kinases,
                        'heatmap': struct_heatmap_data,
                        'kinase_families': struct_family_dist,
                        'site_count': len(struct_similar_sites),
                        'rmsd_threshold': 3.0,
                        'similarity_threshold': 0.6
                    }
                
                # 2. Sequence network prediction
                sequence_network_kinases = predict_kinases_network(
                    site_id, top_n=10, score_type='sequence',
                    similarity_threshold=0.6, rmsd_threshold=3.0
                )
                
                if sequence_network_kinases:
                    # Get similar sites
                    seq_similar_sites = get_similar_sites(
                        site_id, uniprot_id, site_data,
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Get heatmap data with enhanced function
                    seq_heatmap_data = get_network_heatmap_data(
                        site_id, top_n=10, score_type='sequence',
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Log heatmap data for debugging
                    logger.info(f"Sequence heatmap data: {len(seq_heatmap_data['sites'])} sites, " 
                                f"{len(seq_heatmap_data['kinases'])} kinases, "
                                f"{len(seq_heatmap_data['scores'])} scores")
                    
                    # Get family distribution
                    seq_family_dist = get_kinase_family_distribution_network(
                        site_id, score_type='sequence',
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Combine into data structure
                    sequence_network_kinase_data = {
                        'site_id': site_id,
                        'top_kinases': sequence_network_kinases,
                        'heatmap': seq_heatmap_data,
                        'kinase_families': seq_family_dist,
                        'site_count': len(seq_similar_sites),
                        'rmsd_threshold': 3.0,
                        'similarity_threshold': 0.6
                    }
            except Exception as network_e:
                logger.error(f"Error in network kinase prediction: {network_e}")
                import traceback
                logger.error(traceback.format_exc())
                
                # Provide minimal fallback data to avoid template errors
                structure_network_kinase_data = {
                    'site_id': site_id,
                    'top_kinases': [],
                    'heatmap': {'sites': [], 'kinases': [], 'scores': []},
                    'kinase_families': {},
                    'site_count': 0,
                    'rmsd_threshold': 3.0,
                    'similarity_threshold': 0.6,
                    'error': str(network_e)
                }
                
                sequence_network_kinase_data = {
                    'site_id': site_id,
                    'top_kinases': [],
                    'heatmap': {'sites': [], 'kinases': [], 'scores': []},
                    'kinase_families': {},
                    'site_count': 0,
                    'rmsd_threshold': 3.0,
                    'similarity_threshold': 0.6,
                    'error': str(network_e)
                }

        except Exception as e:
            logger.error(f"Error processing kinase data: {e}")
            import traceback
            logger.error(traceback.format_exc())
        
        # Render the template with all the data
        return render_template(
            'site.html',
            protein=protein_data,
            site=site,
            site_data=site_data,
            site_id=site_id,  # Make sure site_id is passed to the template
            structure_html=structure_html,
            structural_matches=structural_matches,
            supplementary_data=site_data,
            network_html=network_html,
            motif_html=motif_html,
            distribution_data=distribution_data,
            enhanced_3d_html=enhanced_3d_html,
            context_data=context_data,
            # Sequence similarity data:
            sequence_matches=sequence_matches,
            sequence_network_data=sequence_network_data,
            sequence_conservation=sequence_conservation,
            sequence_motif_html=sequence_motif_html,
            # Standard kinase prediction data:
            structure_kinase_data=structure_kinase_data,
            sequence_kinase_data=sequence_kinase_data,
            # Network kinase prediction data:
            structure_network_kinase_data=structure_network_kinase_data,
            sequence_network_kinase_data=sequence_network_kinase_data
        )
    except Exception as e:
        logger.error(f"Error in site detail view: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return render_template('error.html', error=str(e))

@app.route('/site-search', methods=['GET', 'POST'])
def site_search():
    """Search for a specific phosphorylation site."""
    if request.method == 'POST':
        # Get form data
        uniprot_id = request.form.get('uniprot_id', '')
        site = request.form.get('site', '')
        
        if not uniprot_id:
            return render_template('site-search.html', error="Please enter a UniProt ID")
        
        if not site:
            return render_template('site-search.html', error="Please enter a site identifier")
            
        # Validate site format (S, T, or Y followed by a number)
        import re
        if not re.match(r'^[STY]\d+', site):
            return render_template('site-search.html', error="Invalid site format. Expected format: S123, T45, Y678, etc.")
            
        # Redirect to site detail page
        return redirect(url_for('site_detail', uniprot_id=uniprot_id, site=site))
    
    # GET request - show search form
    return render_template('site-search.html')

# Add a new API endpoint for kinase data
@app.route('/api/kinases/<site_id>', methods=['GET'])
def api_kinases(site_id):
    """API endpoint for kinase prediction scores."""
    try:
        score_type = request.args.get('type', 'structure')
        top_n = int(request.args.get('top_n', 10))
        
        # Get kinase scores
        scores = get_site_kinase_scores(site_id, score_type)
        if not scores or 'scores' not in scores:
            return jsonify({'error': f'No {score_type} kinase scores found for {site_id}'}), 404
        
        # Get top predicted kinases
        top_kinases = predict_kinases(site_id, top_n, score_type)
        
        # Get known kinase info
        known_kinase = get_known_kinase_info(site_id, score_type)
        
        # Get comparison with the other score type
        other_type = 'sequence' if score_type == 'structure' else 'structure'
        comparison = get_kinase_comparison_data(site_id, [score_type, other_type], top_n)
        
        # Return all data
        return jsonify({
            'site_id': site_id,
            'score_type': score_type,
            'known_kinase': known_kinase,
            'top_kinases': top_kinases,
            'comparison': comparison
        })
    except Exception as e:
        logger.error(f"Error in kinase API: {e}")
        return jsonify({'error': str(e)}), 500

# Add another API endpoint for comparing kinases between sites
@app.route('/api/kinases/compare', methods=['POST'])
def api_compare_kinases():
    """API endpoint for comparing kinase predictions across multiple sites."""
    try:
        data = request.get_json()
        if not data or 'site_ids' not in data:
            return jsonify({'error': 'Missing site_ids parameter'}), 400
        
        site_ids = data['site_ids']
        score_type = data.get('type', 'structure')
        top_n = int(data.get('top_n', 10))
        
        # Get heatmap data
        heatmap_data = get_heatmap_data(site_ids, top_n, score_type)
        
        # Return the data
        return jsonify({
            'heatmap': heatmap_data
        })
    except Exception as e:
        logger.error(f"Error in kinase comparison API: {e}")
        return jsonify({'error': str(e)}), 500
    
# Add a new API endpoint for network-based kinase predictions
@app.route('/api/network-kinases/<site_id>', methods=['GET'])
def api_network_kinases(site_id):
    """API endpoint for network-based kinase prediction scores."""
    try:
        score_type = request.args.get('type', 'structure')
        top_n = int(request.args.get('top_n', 10))
        similarity_threshold = float(request.args.get('similarity_threshold', 0.6))
        rmsd_threshold = float(request.args.get('rmsd_threshold', 3.0))
        
        # Get network-based kinase predictions
        network_kinases = predict_kinases_network(
            site_id, top_n=top_n, score_type=score_type,
            similarity_threshold=similarity_threshold, 
            rmsd_threshold=rmsd_threshold
        )
        
        if not network_kinases:
            return jsonify({'error': f'No {score_type} network kinase scores found for {site_id}'}), 404
        
        # Get similar sites
        similar_sites = get_similar_sites(
            site_id,
            similarity_threshold=similarity_threshold, 
            rmsd_threshold=rmsd_threshold
        )
        
        # Get heatmap data
        heatmap_data = get_network_heatmap_data(
            site_id, top_n=top_n, score_type=score_type,
            similarity_threshold=similarity_threshold, 
            rmsd_threshold=rmsd_threshold
        )
        
        # Get kinase families
        kinase_families = get_kinase_family_distribution_network(
            site_id, score_type=score_type,
            similarity_threshold=similarity_threshold, 
            rmsd_threshold=rmsd_threshold
        )
        
        # Return all data
        return jsonify({
            'site_id': site_id,
            'score_type': score_type,
            'top_kinases': network_kinases,
            'heatmap': heatmap_data,
            'kinase_families': kinase_families,
            'site_count': len(similar_sites),
            'rmsd_threshold': rmsd_threshold,
            'similarity_threshold': similarity_threshold
        })
    except Exception as e:
        logger.error(f"Error in network kinase API: {e}")
        return jsonify({'error': str(e)}), 500

# Add an endpoint to update network predictions with new thresholds
@app.route('/api/update-network-kinases/<site_id>', methods=['POST'])
def update_network_kinases(site_id):
    """Update network-based kinase predictions with new thresholds with enhanced error handling."""
    try:
        data = request.get_json()
        if not data:
            return jsonify({'error': 'No data provided'}), 400
        
        score_type = data.get('score_type', 'structure')
        top_n = int(data.get('top_n', 10))
        
        # Get appropriate threshold based on score_type
        if score_type == 'structure':
            threshold_param = 'rmsd_threshold'
            threshold_default = 3.0
        else:  # sequence
            threshold_param = 'similarity_threshold'
            threshold_default = 0.6
            
        threshold = float(data.get(threshold_param, threshold_default))
        
        logger.info(f"Updating {score_type} network kinase predictions with {threshold_param}={threshold}")
        
        # Common parameters for all function calls
        kwargs = {
            'site_id': site_id,
            'top_n': top_n,
            'score_type': score_type
        }
        
        # Add the appropriate threshold parameter
        if score_type == 'structure':
            kwargs['rmsd_threshold'] = threshold
            # Use default for similarity
            kwargs['similarity_threshold'] = 0.6
        else:
            kwargs['similarity_threshold'] = threshold
            # Use default for RMSD
            kwargs['rmsd_threshold'] = 3.0
        
        # Get updated predictions
        network_kinases = predict_kinases_network(**kwargs)
        
        # Get updated similar sites
        similar_sites = get_similar_sites(site_id, **{k: v for k, v in kwargs.items() if k != 'top_n'})
        
        # Get updated heatmap data with enhanced function
        try:
            heatmap_data = get_network_heatmap_data(**kwargs)
            
            # Log heatmap data stats to debug any issues
            logger.info(f"Heatmap data: {len(heatmap_data['sites'])} sites, {len(heatmap_data['kinases'])} kinases, {len(heatmap_data['scores'])} scores")
            
            # Add family distribution for radar chart
            family_distribution = get_kinase_family_distribution_network(
                site_id, score_type=score_type,
                similarity_threshold=kwargs.get('similarity_threshold', 0.6),
                rmsd_threshold=kwargs.get('rmsd_threshold', 3.0)
            )
            
            # Return updated data with full response
            result = {
                'site_id': site_id,
                'score_type': score_type,
                'top_kinases': network_kinases,
                'heatmap': heatmap_data,
                'kinase_families': family_distribution,
                'site_count': len(similar_sites),
                'rmsd_threshold': kwargs.get('rmsd_threshold'),
                'similarity_threshold': kwargs.get('similarity_threshold')
            }
            
            return jsonify(result)
            
        except Exception as heatmap_error:
            # Return basic data if heatmap fails
            logger.error(f"Error generating heatmap: {heatmap_error}")
            
            return jsonify({
                'site_id': site_id,
                'score_type': score_type,
                'top_kinases': network_kinases,
                'heatmap': {
                    'sites': [site_id],  # Include at least the query site
                    'kinases': [k['kinase'] for k in network_kinases[:top_n]],
                    'scores': []  # Empty scores as fallback
                },
                'site_count': len(similar_sites),
                'rmsd_threshold': kwargs.get('rmsd_threshold'),
                'similarity_threshold': kwargs.get('similarity_threshold'),
                'error': f"Heatmap generation failed: {str(heatmap_error)}"
            })
            
    except Exception as e:
        logger.error(f"Error updating network kinase predictions: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return jsonify({
            'error': str(e),
            'traceback': traceback.format_exc()
        }), 500

@app.route('/faq')
def faq():
    """Render the FAQ page."""
    return render_template('faq.html')

if __name__ == '__main__':
    app.run(debug=True)

