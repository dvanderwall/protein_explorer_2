"""
Flask web application for the KinoPlex protein explorer application.

This version uses database queries instead of file loading, but maintains
the same functionality and output as the original application.
"""

from flask import Flask, render_template, request, jsonify, redirect, url_for
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

# Import database module
from protein_explorer.db import (
    get_phosphosite_data, get_phosphosites_batch,
    find_structural_matches, find_structural_matches_batch,
    find_sequence_matches, find_sequence_matches_batch,
    get_kinase_scores, get_kinase_scores_batch
)

# Add the parent directory to the path to import protein_explorer
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import protein_explorer as pe

from protein_explorer.data.scaffold import (
    get_uniprot_id_from_gene,
    get_protein_by_id
)

# Import phosphosite analysis functions
from protein_explorer.analysis.phospho import analyze_phosphosites
from protein_explorer.analysis.phospho_analyzer import (
    get_phosphosite_data,
    enhance_phosphosite, 
    enhance_structural_matches,
    analyze_phosphosite_context,
    enhance_site_visualization, 
    create_comparative_motif_visualization,
    analyze_residue_distributions
)

# Import visualization functions
from protein_explorer.visualization.network import create_phosphosite_network

from protein_explorer.analysis.sequence_analyzer import (
    find_sequence_matches,
    analyze_motif_conservation,
    create_sequence_network_data,
    get_motif_enrichment,
    create_sequence_motif_visualization
)

from protein_explorer.analysis.kinase_predictor import (
    load_kinase_scores, 
    get_site_kinase_scores,
    predict_kinases,
    get_heatmap_data,
    get_kinase_comparison_data,
    get_known_kinase_info,
    categorize_kinases_by_family
)

# Import the new network-based functions
from protein_explorer.analysis.network_kinase_predictor import (
    get_similar_sites, compute_aggregated_kinase_scores, predict_kinases_network,
    get_network_heatmap_data, get_network_kinase_comparison, 
    get_kinase_family_distribution_network
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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
            all_matches = find_structural_matches_batch(site_ids, rmsd_threshold=5.0)
            
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
    """Phosphosite structural analysis page with improved error handling and supplementary data."""
    import os
    from flask import request, render_template
    
    # Import enhanced table function
    from protein_explorer.analysis.enhanced_table import enhance_phosphosite_table
    
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
                'error': None
            }
            
            # Try to get phosphosites
            try:
                # Get sequence
                sequence = protein_data.get('metadata', {}).get('sequence', {}).get('value')
                if not sequence:
                    raise ValueError(f"Protein sequence not found for {uniprot_id}")
                    
                # Get structure
                structure = pe.data.get_alphafold_structure(uniprot_id)
                if not structure:
                    raise ValueError(f"Protein structure not found for {uniprot_id}")
                
                # Analyze phosphosites
                phosphosites = pe.analysis.phospho.analyze_phosphosites(sequence, structure, uniprot_id)
                
                # Get site IDs for batch query
                site_ids = [f"{uniprot_id}_{site['resno']}" for site in phosphosites]
                
                # Batch query for supplementary data
                supp_data_dict = get_phosphosites_batch(site_ids)
                
                # Enhance phosphosites with database data
                enhanced_sites = []
                for site in phosphosites:
                    site_id = f"{uniprot_id}_{site['resno']}"
                    supp_data = supp_data_dict.get(site_id)
                    
                    if supp_data:
                        # Create enhanced site
                        enhanced_site = site.copy()
                        
                        # Add supplementary data fields
                        for key, value in supp_data.items():
                            if value is not None:
                                # Map database field names to our field names
                                if key == 'SITE_+/-7_AA' and 'motif' not in enhanced_site:
                                    enhanced_site['motif'] = value
                                elif key == 'site_plddt' and 'mean_plddt' not in enhanced_site:
                                    enhanced_site['mean_plddt'] = value
                                elif key == 'nearby_count' and 'nearby_count' not in enhanced_site:
                                    enhanced_site['nearby_count'] = value
                                elif key not in enhanced_site:
                                    enhanced_site[key] = value
                        
                        enhanced_sites.append(enhanced_site)
                    else:
                        enhanced_sites.append(site)
                
                results['phosphosites'] = enhanced_sites
                
                # Try to find structural matches
                try:
                    # Batch query for structural matches
                    all_matches = find_structural_matches_batch(site_ids, rmsd_threshold=5.0)
                    
                    # Organize by site
                    structural_matches = {}
                    for site in phosphosites:
                        site_id = f"{uniprot_id}_{site['resno']}"
                        if site_id in all_matches:
                            matches = all_matches[site_id]
                            # Skip empty lists
                            if matches:
                                # Enhance with supplementary data
                                enhanced_matches = enhance_structural_matches(matches, site['site'])
                                structural_matches[site['site']] = enhanced_matches
                    
                    results['structural_matches'] = structural_matches
                except Exception as e:
                    logger.error(f"Error analyzing structural matches: {e}")
                    results['error'] = f"Error analyzing structural matches: {str(e)}"
            except Exception as e:
                logger.error(f"Error analyzing phosphosites: {e}")
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
            
            # Return results including enhanced data
            return render_template('phosphosite.html', 
                                  protein_info=results['protein_info'],
                                  phosphosites=results['phosphosites'],
                                  phosphosites_html=results.get('phosphosites_html'),
                                  structural_matches=results['structural_matches'],
                                  error=results.get('error'))
                
        except Exception as e:
            print(f"Error in phosphosite analysis: {e}")
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
            if structural_matches and 'motif' in site_data:
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
            sequence_matches = find_sequence_matches(site_id, top_n=200, min_similarity=0.4)
            logger.info(f"Found {len(sequence_matches)} sequence matches")
            
            # Log a sample of the matches to see what they contain
            if sequence_matches:
                logger.info(f"Sample match: {sequence_matches[0]}")
                
                # Create sequence network data
                sequence_network_data = create_sequence_network_data(
                    site_id, 
                    sequence_matches,
                    query_motif=site_data.get('motif') or site_data.get('SITE_+/-7_AA')
                )
                logger.info(f"Created network data with {len(sequence_network_data['nodes'])} nodes")
                
                # Create sequence conservation analysis
                sequence_conservation = analyze_motif_conservation(