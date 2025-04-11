"""
Cantley_Kinase_Scorer.py

This module processes phosphosite match data to retrieve kinase scores from the Cantley database
and creates visualizations for protein phosphorylation analysis.

Functions:
- process_matches_for_kinase_scoring: Process matches and retrieve kinase scores
- create_heatmap_visualization: Create heatmap of matches and kinases
- create_top_kinases_visualization: Visualize top median kinase scores
- create_kinase_family_visualization: Visualize kinase family distributions
- create_protein_kinase_report: Generate comprehensive kinase report
- get_html_protein_kinase_report: Convert report to HTML for display
- get_proteins_kinase_heatmap: Create multi-protein kinase comparison
"""

import pandas as pd
import numpy as np
import logging
from typing import Dict, List, Tuple, Any, Optional
import json
import re
from collections import Counter

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def process_matches_for_kinase_scoring(matches: List[Dict], match_type: str = 'structural') -> Tuple[Dict, List]:
    """
    Process a list of matches (structural or sequence) for kinase scoring.
    
    Args:
        matches: List of match dictionaries
        match_type: Type of matches ('structural' or 'sequence')
        
    Returns:
        Tuple containing:
        - Dictionary of kinase scores by target_id
        - List of processed matches with kinase score data
    """
    from protein_explorer.db.db import get_cantley_st_kinase_scores_batch, get_cantley_y_kinase_scores_batch
    
    # Skip if no matches
    if not matches:
        logger.warning(f"No {match_type} matches provided")
        return {}, []
    
    # Count residue types to determine majority
    residue_types = []
    for match in matches:
        if 'ResidueType' in match:
            residue_types.append(match['ResidueType'])
        elif 'target_site' in match:
            # Try to extract from target_site (e.g. "S123")
            first_char = match['target_site'][0] if match['target_site'] else None
            if first_char in ['S', 'T', 'Y']:
                residue_types.append(first_char)
    
    # Separate matches by residue type
    st_matches = []
    y_matches = []
    
    # Lists to store target_ids
    st_target_ids = []
    y_target_ids = []
    
    # Process each match
    for match in matches:
        # Skip matches without target_id
        if 'target_id' not in match:
            continue
        
        # Determine residue type
        residue_type = None
        if 'ResidueType' in match:
            residue_type = match['ResidueType']
        elif 'target_site' in match and match['target_site']:
            first_char = match['target_site'][0]
            if first_char in ['S', 'T', 'Y']:
                residue_type = first_char
        
        # If still no residue type, use majority
        if not residue_type:
            counter = Counter(residue_types)
            if counter:
                residue_type = counter.most_common(1)[0][0]
            else:
                # Default to S if we really can't determine
                residue_type = 'S'
        
        # Assign to appropriate list
        if residue_type in ['S', 'T']:
            st_matches.append(match)
            st_target_ids.append(match['target_id'])
        elif residue_type == 'Y':
            y_matches.append(match)
            y_target_ids.append(match['target_id'])
    
    # Log the distribution
    logger.info(f"Processing {len(st_matches)} S/T matches and {len(y_matches)} Y matches")
    
    # Initialize scores dictionary
    all_scores = {}
    
    # Get ST kinase scores if we have ST matches
    if st_target_ids:
        logger.info(f"Getting Cantley ST kinase scores for {len(st_target_ids)} targets")
        st_scores = get_cantley_st_kinase_scores_batch(st_target_ids)
        if st_scores:
            all_scores.update(st_scores)
    
    # Get Y kinase scores if we have Y matches
    if y_target_ids:
        logger.info(f"Getting Cantley Y kinase scores for {len(y_target_ids)} targets")
        y_scores = get_cantley_y_kinase_scores_batch(y_target_ids)
        if y_scores:
            all_scores.update(y_scores)
    
    # Combine matches with scores
    processed_matches = []
    
    # Process ST matches
    for match in st_matches:
        target_id = match['target_id']
        if target_id in all_scores:
            # Create a copy with the scores added
            enhanced_match = match.copy()
            enhanced_match['kinase_scores'] = all_scores[target_id]['scores']
            enhanced_match['score_type'] = 'ST'
            processed_matches.append(enhanced_match)
    
    # Process Y matches
    for match in y_matches:
        target_id = match['target_id']
        if target_id in all_scores:
            # Create a copy with the scores added
            enhanced_match = match.copy()
            enhanced_match['kinase_scores'] = all_scores[target_id]['scores']
            enhanced_match['score_type'] = 'Y'
            processed_matches.append(enhanced_match)
    
    logger.info(f"Processed {len(processed_matches)} matches with kinase scores out of {len(matches)} total matches")
    
    return all_scores, processed_matches

def create_heatmap_visualization(site_id: str, processed_matches: List[Dict], top_n_kinases: int = 20, 
                                max_matches: int = 50) -> str:
    """
    Create a heatmap visualization of kinase scores for a phosphosite and its matches.
    
    Args:
        site_id: The ID of the query phosphosite
        processed_matches: List of matches with kinase scores
        top_n_kinases: Number of top kinases to display
        max_matches: Maximum number of matches to include in the heatmap
        
    Returns:
        HTML string with the heatmap visualization
    """
    from protein_explorer.db.db import get_cantley_st_kinase_scores, get_cantley_y_kinase_scores
    
    # Get residue type from site_id
    residue_type = None
    match = re.search(r'_([A-Z])(\d+)', site_id)
    if match:
        residue_type = match.group(1)
    
    # If no residue type in site_id, try to get from first match
    if not residue_type and processed_matches:
        residue_type = processed_matches[0].get('ResidueType')
    
    # Default to ST if still no residue type
    if not residue_type or residue_type not in ['S', 'T', 'Y']:
        residue_type = 'S'
    
    # Get query site kinase scores
    if residue_type == 'Y':
        query_scores = get_cantley_y_kinase_scores(site_id)
    else:
        query_scores = get_cantley_st_kinase_scores(site_id)
    
    # Prepare data for heatmap
    all_kinase_scores = {}
    site_labels = []
    
    # Add query site scores if available
    if query_scores and 'scores' in query_scores:
        all_kinase_scores[site_id] = query_scores['scores']
        site_labels.append(site_id)
    
    # Add match scores
    for match in processed_matches[:max_matches]:
        if 'target_id' in match and 'kinase_scores' in match:
            target_id = match['target_id']
            all_kinase_scores[target_id] = match['kinase_scores']
            site_labels.append(target_id)
    
    # If no scores, return message
    if not all_kinase_scores:
        return '<div class="alert alert-warning">No kinase score data available for heatmap.</div>'
    
    # Convert to DataFrame for easier processing
    try:
        scores_df = pd.DataFrame(all_kinase_scores).T
        
        # Calculate median scores for each kinase
        median_scores = scores_df.median(axis=0)
        
        # Get top N kinases by median score
        top_kinases = median_scores.nlargest(top_n_kinases).index.tolist()
        
        # Subset the DataFrame to include only top kinases
        top_scores_df = scores_df[top_kinases]
        
        # Prepare data for the heatmap
        heatmap_data = {
            'sites': site_labels,
            'kinases': top_kinases,
            'scores': []
        }
        
        # Populate scores in the format expected by the frontend
        for site in site_labels:
            for kinase in top_kinases:
                score = float(top_scores_df.loc[site, kinase]) if site in top_scores_df.index and kinase in top_scores_df.columns else 0
                heatmap_data['scores'].append({
                    'site': site,
                    'kinase': kinase,
                    'score': score
                })
    except Exception as e:
        logger.error(f"Error creating heatmap data: {e}")
        return f'<div class="alert alert-danger">Error creating heatmap: {str(e)}</div>'
    
    # Create HTML for heatmap visualization
    html = f"""
    <div class="card">
        <div class="card-header">
            <h5 class="mb-0">Kinase Prediction Heatmap</h5>
        </div>
        <div class="card-body">
            <div id="kinase-heatmap" style="height: 600px;"></div>
        </div>
    </div>
    
    <script>
        document.addEventListener('DOMContentLoaded', function() {{
            // Heatmap data
            const heatmapData = {json.dumps(heatmap_data)};
            
            // Create the heatmap
            createKinaseHeatmap('kinase-heatmap', heatmapData);
        }});
        
        function createKinaseHeatmap(elementId, data) {{
            const sites = data.sites;
            const kinases = data.kinases;
            const scores = data.scores;
            
            // Prepare data for Plotly
            const x = [];
            const y = [];
            const z = [];
            const text = [];
            
            // Create a 2D array for the z values
            for (let i = 0; i < sites.length; i++) {{
                const row = [];
                const textRow = [];
                for (let j = 0; j < kinases.length; j++) {{
                    // Find the score for this site and kinase
                    const scoreObj = scores.find(s => s.site === sites[i] && s.kinase === kinases[j]);
                    const score = scoreObj ? scoreObj.score : 0;
                    row.push(score);
                    textRow.push(`Site: ${{sites[i]}}<br>Kinase: ${{kinases[j]}}<br>Score: ${{score.toFixed(3)}}`);
                }}
                z.push(row);
                text.push(textRow);
            }}
            
            // Create the heatmap trace
            const trace = {{
                x: kinases,
                y: sites,
                z: z,
                text: text,
                type: 'heatmap',
                colorscale: 'Viridis',
                hoverinfo: 'text'
            }};
            
            // Layout options
            const layout = {{
                title: 'Kinase Prediction Scores',
                xaxis: {{
                    title: 'Kinase',
                    tickangle: 45
                }},
                yaxis: {{
                    title: 'Phosphosite'
                }},
                margin: {{ l: 150, r: 50, b: 150, t: 100 }}
            }};
            
            // Create the heatmap
            Plotly.newPlot(elementId, [trace], layout);
        }}
    </script>
    """
    
    return html

def create_top_kinases_visualization(site_id: str, processed_matches: List[Dict], top_n: int = 10) -> str:
    """
    Create a visualization of the top kinases for a phosphosite based on median scores from matches.
    
    Args:
        site_id: The ID of the query phosphosite
        processed_matches: List of matches with kinase scores
        top_n: Number of top kinases to display
        
    Returns:
        HTML string with the top kinases visualization
    """
    from protein_explorer.db.db import get_cantley_st_kinase_scores, get_cantley_y_kinase_scores
    
    # Get residue type from site_id
    residue_type = None
    match = re.search(r'_([A-Z])(\d+)', site_id)
    if match:
        residue_type = match.group(1)
    
    # If no residue type in site_id, try to get from first match
    if not residue_type and processed_matches:
        residue_type = processed_matches[0].get('ResidueType')
    
    # Default to S if still no residue type
    if not residue_type or residue_type not in ['S', 'T', 'Y']:
        residue_type = 'S'
    
    # Get query site kinase scores
    if residue_type == 'Y':
        query_scores = get_cantley_y_kinase_scores(site_id)
    else:
        query_scores = get_cantley_st_kinase_scores(site_id)
    
    # Prepare data for visualization
    all_kinase_scores = {}
    
    # Add query site scores if available
    if query_scores and 'scores' in query_scores:
        all_kinase_scores[site_id] = query_scores['scores']
    
    # Add match scores
    for match in processed_matches:
        if 'target_id' in match and 'kinase_scores' in match:
            target_id = match['target_id']
            all_kinase_scores[target_id] = match['kinase_scores']
    
    # If no scores, return message
    if not all_kinase_scores:
        return '<div class="alert alert-warning">No kinase score data available for visualization.</div>'
    
    try:
        # Convert to DataFrame for easier processing
        scores_df = pd.DataFrame(all_kinase_scores).T
        
        # Calculate median scores for each kinase
        median_scores = scores_df.median(axis=0)
        
        # Get top N kinases by median score
        top_kinases_df = median_scores.nlargest(top_n)
        
        # Prepare data for visualization
        kinases = top_kinases_df.index.tolist()
        scores = top_kinases_df.values.tolist()
        
    except Exception as e:
        logger.error(f"Error creating top kinases data: {e}")
        return f'<div class="alert alert-danger">Error creating top kinases visualization: {str(e)}</div>'
    
    # Create HTML for visualization
    html = f"""
    <div class="card">
        <div class="card-header">
            <h5 class="mb-0">Top {top_n} Predicted Kinases</h5>
        </div>
        <div class="card-body">
            <div id="top-kinases-chart" style="height: 400px;"></div>
        </div>
    </div>
    
    <script>
        document.addEventListener('DOMContentLoaded', function() {{
            // Data for chart
            const kinases = {json.dumps(kinases)};
            const scores = {json.dumps(scores)};
            
            // Create the chart
            createTopKinasesChart('top-kinases-chart', kinases, scores);
        }});
        
        function createTopKinasesChart(elementId, kinases, scores) {{
            // Create the bar chart trace
            const trace = {{
                x: kinases,
                y: scores,
                type: 'bar',
                marker: {{
                    color: 'rgba(50, 171, 96, 0.7)',
                    line: {{
                        color: 'rgba(50, 171, 96, 1.0)',
                        width: 1
                    }}
                }}
            }};
            
            // Layout options
            const layout = {{
                title: 'Top Kinases by Median Score',
                xaxis: {{
                    title: 'Kinase',
                    tickangle: 45
                }},
                yaxis: {{
                    title: 'Median Score'
                }},
                margin: {{ l: 50, r: 50, b: 100, t: 100 }}
            }};
            
            // Create the chart
            Plotly.newPlot(elementId, [trace], layout);
        }}
    </script>
    """
    
    return html

def create_protein_kinase_report(uniprot_id: str, phosphosites: List[Dict], 
                                 structural_matches: Dict, sequence_matches: Dict) -> Dict:
    """
    Generate a comprehensive kinase prediction report for a protein.
    
    Args:
        uniprot_id: The UniProt ID of the protein
        phosphosites: List of phosphosite dictionaries
        structural_matches: Dictionary of structural matches (site_name -> list of matches)
        sequence_matches: Dictionary of sequence matches (site_name -> list of matches)
        
    Returns:
        Dictionary containing the report data and visualizations
    """
    # Initialize report data
    report = {
        'uniprot_id': uniprot_id,
        'phosphosite_count': len(phosphosites),
        'has_structural_matches': bool(structural_matches),
        'has_sequence_matches': bool(sequence_matches),
        'site_reports': {},
        'structural_data': None,
        'sequence_data': None
    }
    
    # Process structural matches for all sites
    all_structural_processed = []
    if structural_matches:
        for site_name, matches in structural_matches.items():
            if matches:
                # Generate site ID
                site_id = f"{uniprot_id}_{site_name.lstrip('STY')}"
                
                # Process matches
                _, processed = process_matches_for_kinase_scoring(matches, match_type='structural')
                
                # Store in site_reports
                report['site_reports'][site_name] = {
                    'site_id': site_id,
                    'structural_matches': len(matches),
                    'structural_processed': len(processed),
                    'has_kinase_data': bool(processed)
                }
                
                # Add to all processed matches
                all_structural_processed.extend(processed)
    
    # Process sequence matches for all sites
    all_sequence_processed = []
    if sequence_matches:
        for site_name, matches in sequence_matches.items():
            if matches:
                # Generate site ID
                site_id = f"{uniprot_id}_{site_name.lstrip('STY')}"
                
                # Process matches
                _, processed = process_matches_for_kinase_scoring(matches, match_type='sequence')
                
                # Update or create site report
                if site_name in report['site_reports']:
                    report['site_reports'][site_name].update({
                        'sequence_matches': len(matches),
                        'sequence_processed': len(processed),
                        'has_sequence_kinase_data': bool(processed)
                    })
                else:
                    report['site_reports'][site_name] = {
                        'site_id': site_id,
                        'sequence_matches': len(matches),
                        'sequence_processed': len(processed),
                        'has_sequence_kinase_data': bool(processed)
                    }
                
                # Add to all processed matches
                all_sequence_processed.extend(processed)
    
    # Store all processed matches for visualization
    if all_structural_processed:
        # Choose a representative site for visualization
        rep_site = next(iter(report['site_reports'].keys()))
        rep_site_id = f"{uniprot_id}_{rep_site.lstrip('STY')}"
        
        report['structural_data'] = {
            'site_id': rep_site_id,
            'processed_matches': all_structural_processed
        }
    
    if all_sequence_processed:
        # Choose a representative site for visualization
        rep_site = next(iter(report['site_reports'].keys()))
        rep_site_id = f"{uniprot_id}_{rep_site.lstrip('STY')}"
        
        report['sequence_data'] = {
            'site_id': rep_site_id,
            'processed_matches': all_sequence_processed
        }
    
    return report

def get_html_protein_kinase_report(report: Dict) -> str:
    """
    Generate HTML representation of the protein kinase report.
    
    Args:
        report: Dictionary containing the report data from create_protein_kinase_report
        
    Returns:
        HTML string for rendering the report
    """
    # Create visualizations based on report data
    heatmap_html = ""
    top_kinases_html = ""
    
    if report.get('structural_data'):
        struct_data = report['structural_data']
        
        # Create structural heatmap
        heatmap_html = create_heatmap_visualization(
            struct_data['site_id'], 
            struct_data['processed_matches'], 
            top_n_kinases=20, 
            max_matches=50
        )
        
        # Create top kinases visualization
        top_kinases_html = create_top_kinases_visualization(
            struct_data['site_id'], 
            struct_data['processed_matches'], 
            top_n=10
        )
    
    # Create HTML for site list
    site_list_html = """
    <div class="table-responsive">
        <table class="table table-striped">
            <thead>
                <tr>
                    <th>Site</th>
                    <th>Structural Matches</th>
                    <th>With Kinase Data</th>
                    <th>Sequence Matches</th>
                    <th>With Kinase Data</th>
                    <th>Actions</th>
                </tr>
            </thead>
            <tbody>
    """
    
    # Add site rows
    for site_name, site_data in report['site_reports'].items():
        site_list_html += f"""
            <tr>
                <td>{site_name}</td>
                <td>{site_data.get('structural_matches', 0)}</td>
                <td>{site_data.get('structural_processed', 0)}</td>
                <td>{site_data.get('sequence_matches', 0)}</td>
                <td>{site_data.get('sequence_processed', 0)}</td>
                <td>
                    <a href="/site/{report['uniprot_id']}/{site_name}" class="btn btn-sm btn-primary">View</a>
                </td>
            </tr>
        """
    
    site_list_html += """
            </tbody>
        </table>
    </div>
    """
    
    # Complete HTML
    html = f"""
    <div class="protein-kinase-report">
        <h2>Kinase Prediction Report for {report['uniprot_id']}</h2>
        
        <div class="summary-stats mb-4">
            <div class="row">
                <div class="col-md-4">
                    <div class="card">
                        <div class="card-body">
                            <h5 class="card-title">Phosphosites</h5>
                            <p class="card-text">{report['phosphosite_count']}</p>
                        </div>
                    </div>
                </div>
                <div class="col-md-4">
                    <div class="card">
                        <div class="card-body">
                            <h5 class="card-title">With Structural Matches</h5>
                            <p class="card-text">{sum(1 for site in report['site_reports'].values() if site.get('has_kinase_data', False))}</p>
                        </div>
                    </div>
                </div>
                <div class="col-md-4">
                    <div class="card">
                        <div class="card-body">
                            <h5 class="card-title">With Sequence Matches</h5>
                            <p class="card-text">{sum(1 for site in report['site_reports'].values() if site.get('has_sequence_kinase_data', False))}</p>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        
        <ul class="nav nav-tabs mb-4" id="kinaseTabs" role="tablist">
            <li class="nav-item">
                <a class="nav-link active" id="heatmap-tab" data-toggle="tab" href="#heatmap" role="tab">Kinase Heatmap</a>
            </li>
            <li class="nav-item">
                <a class="nav-link" id="top-kinases-tab" data-toggle="tab" href="#top-kinases" role="tab">Top Kinases</a>
            </li>
            <li class="nav-item">
                <a class="nav-link" id="sites-tab" data-toggle="tab" href="#sites" role="tab">Site Details</a>
            </li>
        </ul>
        
        <div class="tab-content" id="kinaseTabsContent">
            <div class="tab-pane fade show active" id="heatmap" role="tabpanel">
                {heatmap_html if heatmap_html else '<div class="alert alert-info">No kinase data available for heatmap visualization.</div>'}
            </div>
            
            <div class="tab-pane fade" id="top-kinases" role="tabpanel">
                {top_kinases_html if top_kinases_html else '<div class="alert alert-info">No kinase data available for top kinases visualization.</div>'}
            </div>
            
            <div class="tab-pane fade" id="sites" role="tabpanel">
                {site_list_html}
            </div>
        </div>
    </div>
    
    <script>
        // Initialize tabs when document is ready
        $(function() {{
            $('#kinaseTabs a').on('click', function(e) {{
                e.preventDefault();
                $(this).tab('show');
            }});
        }});
    </script>
    """
    
    return html

def get_proteins_kinase_heatmap(reports: Dict[str, Dict]) -> str:
    """
    Generate a multi-protein kinase comparison heatmap.
    
    Args:
        reports: Dictionary mapping UniProt IDs to their kinase reports from create_protein_kinase_report
        
    Returns:
        HTML string with the multi-protein heatmap
    """
    # Implementation details for comparing kinases across proteins...
    # This would create a heatmap showing kinase predictions for multiple proteins
    
    if not reports:
        return '<div class="alert alert-info">No protein reports available for comparison.</div>'
    
    # Collect all kinase scores across proteins
    protein_kinase_scores = {}
    
    for uniprot_id, report in reports.items():
        # If structural data is available, use it
        if report.get('structural_data') and report['structural_data'].get('processed_matches'):
            processed_matches = report['structural_data']['processed_matches']
            
            # Aggregate kinase scores across all sites for this protein
            all_scores = {}
            for match in processed_matches:
                if 'kinase_scores' in match:
                    for kinase, score in match['kinase_scores'].items():
                        if kinase not in all_scores:
                            all_scores[kinase] = []
                        all_scores[kinase].append(score)
            
            # Calculate median score for each kinase
            median_scores = {}
            for kinase, scores in all_scores.items():
                if scores:
                    median_scores[kinase] = np.median(scores)
            
            if median_scores:
                protein_kinase_scores[uniprot_id] = median_scores
    
    # Check if we have data
    if not protein_kinase_scores:
        return '<div class="alert alert-warning">No kinase score data available for proteins.</div>'
    
    # Find common kinases and calculate averages
    all_kinases = set()
    for scores in protein_kinase_scores.values():
        all_kinases.update(scores.keys())
    
    # Get top kinases by average score
    kinase_avg_scores = {}
    for kinase in all_kinases:
        scores = [protein_scores.get(kinase, 0) for protein_scores in protein_kinase_scores.values()]
        kinase_avg_scores[kinase] = sum(scores) / len(scores)
    
    top_kinases = sorted(kinase_avg_scores.keys(), key=lambda k: kinase_avg_scores[k], reverse=True)[:20]
    
    # Create heatmap data
    heatmap_data = {
        'proteins': list(protein_kinase_scores.keys()),
        'kinases': top_kinases,
        'scores': []
    }
    
    # Populate scores
    for protein in heatmap_data['proteins']:
        protein_scores = protein_kinase_scores.get(protein, {})
        for kinase in top_kinases:
            score = protein_scores.get(kinase, 0)
            heatmap_data['scores'].append({
                'protein': protein,
                'kinase': kinase,
                'score': score
            })
    
    # Create HTML for heatmap
    html = f"""
    <div class="card">
        <div class="card-header">
            <h5 class="mb-0">Protein Kinase Comparison</h5>
        </div>
        <div class="card-body">
            <div id="protein-kinase-heatmap" style="height: 600px;"></div>
        </div>
    </div>
    
    <script>
        document.addEventListener('DOMContentLoaded', function() {{
            // Heatmap data
            const heatmapData = {json.dumps(heatmap_data)};
            
            // Create the heatmap
            createProteinKinaseHeatmap('protein-kinase-heatmap', heatmapData);
        }});
        
        function createProteinKinaseHeatmap(elementId, data) {{
            // Implementation similar to the previous heatmap function
            const proteins = data.proteins;
            const kinases = data.kinases;
            const scores = data.scores;
            
            // Prepare data for Plotly
            const z = Array(proteins.length).fill().map(() => Array(kinases.length).fill(0));
            const text = Array(proteins.length).fill().map(() => Array(kinases.length).fill(''));
            
            // Fill in the data
            scores.forEach(item => {{
                const proteinIndex = proteins.indexOf(item.protein);
                const kinaseIndex = kinases.indexOf(item.kinase);
                if (proteinIndex >= 0 && kinaseIndex >= 0) {{
                    z[proteinIndex][kinaseIndex] = item.score;
                    text[proteinIndex][kinaseIndex] = 
                        `Protein: ${{item.protein}}<br>Kinase: ${{item.kinase}}<br>Score: ${{item.score.toFixed(3)}}`;
                }}
            }});
            
            // Create the trace
            const trace = {{
                x: kinases,
                y: proteins,
                z: z,
                text: text,
                type: 'heatmap',
                colorscale: 'Viridis',
                hoverinfo: 'text'
            }};
            
            const layout = {{
                title: 'Protein Kinase Prediction Comparison',
                xaxis: {{ title: 'Kinase', tickangle: 45 }},
                yaxis: {{ title: 'Protein' }},
                margin: {{ l: 150, r: 50, b: 150, t: 100 }}
            }};
            
            Plotly.newPlot(elementId, [trace], layout);
        }}
    </script>
    """
    
    return html


def process_site_specific_matches(site_name: str, 
                                  structural_matches: List[Dict], 
                                  sequence_matches: List[Dict]) -> Dict:
    """
    Process both structural and sequence matches for a specific site.
    
    Args:
        site_name: The name of the phosphosite (e.g., 'S123')
        structural_matches: List of structural matches for this site
        sequence_matches: List of sequence matches for this site
        
    Returns:
        Dictionary with processed matches and kinase data
    """
    # Process structural matches
    structural_scores, processed_structural = process_matches_for_kinase_scoring(
        structural_matches, match_type='structural'
    )
    
    # Process sequence matches
    sequence_scores, processed_sequence = process_matches_for_kinase_scoring(
        sequence_matches, match_type='sequence'
    )
    
    return {
        'site_name': site_name,
        'structural': {
            'matches': len(structural_matches),
            'processed': len(processed_structural),
            'scores': structural_scores,
            'processed_matches': processed_structural
        },
        'sequence': {
            'matches': len(sequence_matches),
            'processed': len(processed_sequence),
            'scores': sequence_scores,
            'processed_matches': processed_sequence
        }
    }

def create_site_selector_ui(phosphosites: List[Dict], 
                           structural_matches: Dict, 
                           sequence_matches: Dict,
                           active_site: str = None) -> str:
    """
    Create a UI for selecting different phosphosites to view their kinase data.
    
    Args:
        phosphosites: List of phosphosite dictionaries
        structural_matches: Dictionary mapping site names to structural matches
        sequence_matches: Dictionary mapping site names to sequence matches
        active_site: Currently selected site (if any)
        
    Returns:
        HTML string with the site selector UI
    """
    # Create options for dropdown
    options_html = ""
    
    for site in phosphosites:
        site_name = site.get('site')
        if not site_name:
            continue
            
        # Check if site has matches
        has_structural = site_name in structural_matches and structural_matches[site_name]
        has_sequence = site_name in sequence_matches and sequence_matches[site_name]
        
        # Skip sites with no matches
        if not has_structural and not has_sequence:
            continue
            
        # Create option with data attributes
        selected = "selected" if site_name == active_site else ""
        options_html += f"""
            <option value="{site_name}" 
                    data-structural="{len(structural_matches.get(site_name, []))}" 
                    data-sequence="{len(sequence_matches.get(site_name, []))}"
                    {selected}>
                {site_name} - Structural: {len(structural_matches.get(site_name, []))}, 
                Sequence: {len(sequence_matches.get(site_name, []))}
            </option>
        """
    
    # Create HTML for selector
    html = f"""
    <div class="card mb-4">
        <div class="card-header">
            <h5 class="mb-0">Select Phosphosite</h5>
        </div>
        <div class="card-body">
            <form id="site-selector-form">
                <div class="form-group">
                    <label for="site-selector">Choose a phosphosite to analyze:</label>
                    <select class="form-control" id="site-selector" name="site">
                        <option value="">Select a site...</option>
                        {options_html}
                    </select>
                </div>
                <button type="submit" class="btn btn-primary">Analyze Site</button>
            </form>
        </div>
    </div>
    
    <script>
        document.addEventListener('DOMContentLoaded', function() {{
            // Handle form submission
            document.getElementById('site-selector-form').addEventListener('submit', function(e) {{
                e.preventDefault();
                const siteSelector = document.getElementById('site-selector');
                const selectedSite = siteSelector.value;
                
                if (selectedSite) {{
                    // Show loading indicator
                    document.getElementById('site-analysis').innerHTML = 
                        '<div class="text-center"><div class="spinner-border"></div><p>Loading analysis...</p></div>';
                    
                    // Load site analysis
                    fetch(`/api/site_kinase_analysis/${{selectedSite}}`)
                        .then(response => response.json())
                        .then(data => {{
                            // Update visualizations
                            updateSiteAnalysis(data);
                        }})
                        .catch(error => {{
                            document.getElementById('site-analysis').innerHTML = 
                                `<div class="alert alert-danger">Error loading analysis: ${{error.message}}</div>`;
                        }});
                }}
            }});
        }});
        
        function updateSiteAnalysis(data) {{
            // Update the site analysis section with new data
            document.getElementById('site-analysis').innerHTML = data.html;
            
            // Initialize any charts or visualizations
            if (typeof initializeCharts === 'function') {{
                initializeCharts();
            }}
        }}
    </script>
    """
    
    return html

def generate_site_analysis_visualizations(site_data: Dict, site_id: str, uniprot_id: str) -> str:
    """
    Generate HTML visualizations for a specific site's kinase analysis.
    
    Args:
        site_data: Processed site data from process_site_specific_matches
        site_id: Full site ID (UniProtID_residue)
        uniprot_id: UniProt ID of the protein
        
    Returns:
        HTML string with all visualizations
    """
    # Create tabs for different visualization types
    tabs_html = """
    <ul class="nav nav-tabs" id="siteAnalysisTabs" role="tablist">
        <li class="nav-item">
            <a class="nav-link active" id="structural-tab" data-bs-toggle="tab" href="#structural" role="tab">
                Structural Similarity
            </a>
        </li>
        <li class="nav-item">
            <a class="nav-link" id="sequence-tab" data-bs-toggle="tab" href="#sequence" role="tab">
                Sequence Similarity
            </a>
        </li>
        <li class="nav-item">
            <a class="nav-link" id="combined-tab" data-bs-toggle="tab" href="#combined" role="tab">
                Combined Analysis
            </a>
        </li>
    </ul>
    """
    
    # Create visualizations for structural matches
    structural_html = ""
    if site_data['structural']['processed'] > 0:
        # Create heatmap with unique div ID
        structural_heatmap_div = f"structural-heatmap-{site_id.replace('_', '-')}"
        structural_top_kinases_div = f"structural-top-kinases-{site_id.replace('_', '-')}"
        
        # Create the data JSON for heatmap
        processed_matches = site_data['structural']['processed_matches']
        
        # Prepare heatmap data
        heatmap_data = prepare_heatmap_data(processed_matches, top_n_kinases=15)
        heatmap_json = json.dumps(heatmap_data)
        
        # Prepare top kinases data
        top_kinases_data = prepare_top_kinases_data(processed_matches, top_n=10)
        top_kinases_json = json.dumps(top_kinases_data)
        
        # Info alert
        structural_html = f"""
        <div class="row">
            <div class="col-12">
                <div class="alert alert-info">
                    Analyzed {site_data['structural']['processed']} structural matches with kinase data 
                    out of {site_data['structural']['matches']} total structural matches.
                </div>
            </div>
            
            <div class="col-12 mb-4">
                <div class="card">
                    <div class="card-header">
                        <h5 class="mb-0">Kinase Prediction Heatmap</h5>
                    </div>
                    <div class="card-body">
                        <div id="{structural_heatmap_div}" style="height: 600px;"></div>
                    </div>
                </div>
                
                <script>
                    (function() {{
                        const heatmapData = {heatmap_json};
                        
                        // Create the heatmap when document is ready
                        function initializeHeatmap() {{
                            const sites = heatmapData.sites;
                            const kinases = heatmapData.kinases;
                            const scores = heatmapData.scores;
                            
                            // Prepare data for Plotly
                            const x = kinases;
                            const y = sites;
                            const z = [];
                            const text = [];
                            
                            // Create z and text arrays for heatmap
                            for (let i = 0; i < sites.length; i++) {{
                                const zRow = [];
                                const textRow = [];
                                for (let j = 0; j < kinases.length; j++) {{
                                    // Find the score for this site and kinase
                                    const scoreObj = scores.find(s => 
                                        s.site === sites[i] && s.kinase === kinases[j]);
                                    const score = scoreObj ? scoreObj.score : 0;
                                    zRow.push(score);
                                    textRow.push(`Site: ${{sites[i]}}<br>Kinase: ${{kinases[j]}}<br>Score: ${{score.toFixed(3)}}`);
                                }}
                                z.push(zRow);
                                text.push(textRow);
                            }}
                            
                            // Create the heatmap trace
                            const trace = {{
                                x: x,
                                y: y,
                                z: z,
                                text: text,
                                type: 'heatmap',
                                colorscale: 'Viridis',
                                hoverinfo: 'text'
                            }};
                            
                            // Layout options
                            const layout = {{
                                title: 'Kinase Prediction Scores',
                                xaxis: {{
                                    title: 'Kinase',
                                    tickangle: 45
                                }},
                                yaxis: {{
                                    title: 'Phosphosite'
                                }},
                                margin: {{ l: 150, r: 50, b: 150, t: 100 }}
                            }};
                            
                            // Create the heatmap
                            Plotly.newPlot('{structural_heatmap_div}', [trace], layout);
                        }}
                        
                        // Wait for the element to exist in the DOM
                        const checkExist = setInterval(function() {{
                            if (document.getElementById('{structural_heatmap_div}')) {{
                                clearInterval(checkExist);
                                initializeHeatmap();
                            }}
                        }}, 100);
                    }})();
                </script>
            </div>
            
            <div class="col-12">
                <div class="card">
                    <div class="card-header">
                        <h5 class="mb-0">Top 10 Predicted Kinases</h5>
                    </div>
                    <div class="card-body">
                        <div id="{structural_top_kinases_div}" style="height: 400px;"></div>
                    </div>
                </div>
                
                <script>
                    (function() {{
                        const topKinasesData = {top_kinases_json};
                        
                        // Create the chart when document is ready
                        function initializeTopKinasesChart() {{
                            const kinases = topKinasesData.kinases;
                            const scores = topKinasesData.scores;
                            
                            // Create the bar chart trace
                            const trace = {{
                                x: kinases,
                                y: scores,
                                type: 'bar',
                                marker: {{
                                    color: 'rgba(50, 171, 96, 0.7)',
                                    line: {{
                                        color: 'rgba(50, 171, 96, 1.0)',
                                        width: 1
                                    }}
                                }}
                            }};
                            
                            // Layout options
                            const layout = {{
                                title: 'Top Kinases by Median Score',
                                xaxis: {{
                                    title: 'Kinase',
                                    tickangle: 45
                                }},
                                yaxis: {{
                                    title: 'Median Score'
                                }},
                                margin: {{ l: 50, r: 50, b: 100, t: 100 }}
                            }};
                            
                            // Create the chart
                            Plotly.newPlot('{structural_top_kinases_div}', [trace], layout);
                        }}
                        
                        // Wait for the element to exist in the DOM
                        const checkExist = setInterval(function() {{
                            if (document.getElementById('{structural_top_kinases_div}')) {{
                                clearInterval(checkExist);
                                initializeTopKinasesChart();
                            }}
                        }}, 100);
                    }})();
                </script>
            </div>
        </div>
        """
    else:
        structural_html = """
        <div class="alert alert-warning">
            No structural matches with kinase data available for this site.
        </div>
        """
    
    # Similar implementation for sequence and combined tabs...
    # Create visualizations for sequence matches
    sequence_html = ""
    if site_data['sequence']['processed'] > 0:
        # Create heatmap with unique div ID
        sequence_heatmap_div = f"sequence-heatmap-{site_id.replace('_', '-')}"
        sequence_top_kinases_div = f"sequence-top-kinases-{site_id.replace('_', '-')}"
        
        # Create the data JSON for heatmap
        processed_matches = site_data['sequence']['processed_matches']
        
        # Prepare heatmap data
        heatmap_data = prepare_heatmap_data(processed_matches, top_n_kinases=15)
        heatmap_json = json.dumps(heatmap_data)
        
        # Prepare top kinases data
        top_kinases_data = prepare_top_kinases_data(processed_matches, top_n=10)
        top_kinases_json = json.dumps(top_kinases_data)
        
        # Info alert
        sequence_html = f"""
        <div class="row">
            <div class="col-12">
                <div class="alert alert-info">
                    Analyzed {site_data['sequence']['processed']} sequence matches with kinase data 
                    out of {site_data['sequence']['matches']} total sequence matches.
                </div>
            </div>
            
            <div class="col-12 mb-4">
                <div class="card">
                    <div class="card-header">
                        <h5 class="mb-0">Kinase Prediction Heatmap</h5>
                    </div>
                    <div class="card-body">
                        <div id="{sequence_heatmap_div}" style="height: 600px;"></div>
                    </div>
                </div>
                
                <script>
                    (function() {{
                        const heatmapData = {heatmap_json};
                        
                        // Create the heatmap when document is ready
                        function initializeHeatmap() {{
                            const sites = heatmapData.sites;
                            const kinases = heatmapData.kinases;
                            const scores = heatmapData.scores;
                            
                            // Prepare data for Plotly
                            const x = kinases;
                            const y = sites;
                            const z = [];
                            const text = [];
                            
                            // Create z and text arrays for heatmap
                            for (let i = 0; i < sites.length; i++) {{
                                const zRow = [];
                                const textRow = [];
                                for (let j = 0; j < kinases.length; j++) {{
                                    // Find the score for this site and kinase
                                    const scoreObj = scores.find(s => 
                                        s.site === sites[i] && s.kinase === kinases[j]);
                                    const score = scoreObj ? scoreObj.score : 0;
                                    zRow.push(score);
                                    textRow.push(`Site: ${{sites[i]}}<br>Kinase: ${{kinases[j]}}<br>Score: ${{score.toFixed(3)}}`);
                                }}
                                z.push(zRow);
                                text.push(textRow);
                            }}
                            
                            // Create the heatmap trace
                            const trace = {{
                                x: x,
                                y: y,
                                z: z,
                                text: text,
                                type: 'heatmap',
                                colorscale: 'Viridis',
                                hoverinfo: 'text'
                            }};
                            
                            // Layout options
                            const layout = {{
                                title: 'Kinase Prediction Scores',
                                xaxis: {{
                                    title: 'Kinase',
                                    tickangle: 45
                                }},
                                yaxis: {{
                                    title: 'Phosphosite'
                                }},
                                margin: {{ l: 150, r: 50, b: 150, t: 100 }}
                            }};
                            
                            // Create the heatmap
                            Plotly.newPlot('{sequence_heatmap_div}', [trace], layout);
                        }}
                        
                        // Wait for the element to exist in the DOM
                        const checkExist = setInterval(function() {{
                            if (document.getElementById('{sequence_heatmap_div}')) {{
                                clearInterval(checkExist);
                                initializeHeatmap();
                            }}
                        }}, 100);
                    }})();
                </script>
            </div>
            
            <div class="col-12">
                <div class="card">
                    <div class="card-header">
                        <h5 class="mb-0">Top 10 Predicted Kinases</h5>
                    </div>
                    <div class="card-body">
                        <div id="{sequence_top_kinases_div}" style="height: 400px;"></div>
                    </div>
                </div>
                
                <script>
                    (function() {{
                        const topKinasesData = {top_kinases_json};
                        
                        // Create the chart when document is ready
                        function initializeTopKinasesChart() {{
                            const kinases = topKinasesData.kinases;
                            const scores = topKinasesData.scores;
                            
                            // Create the bar chart trace
                            const trace = {{
                                x: kinases,
                                y: scores,
                                type: 'bar',
                                marker: {{
                                    color: 'rgba(64, 83, 196, 0.7)',  // Different color from structural
                                    line: {{
                                        color: 'rgba(64, 83, 196, 1.0)',
                                        width: 1
                                    }}
                                }}
                            }};
                            
                            // Layout options
                            const layout = {{
                                title: 'Top Kinases by Median Score',
                                xaxis: {{
                                    title: 'Kinase',
                                    tickangle: 45
                                }},
                                yaxis: {{
                                    title: 'Median Score'
                                }},
                                margin: {{ l: 50, r: 50, b: 100, t: 100 }}
                            }};
                            
                            // Create the chart
                            Plotly.newPlot('{sequence_top_kinases_div}', [trace], layout);
                        }}
                        
                        // Wait for the element to exist in the DOM
                        const checkExist = setInterval(function() {{
                            if (document.getElementById('{sequence_top_kinases_div}')) {{
                                clearInterval(checkExist);
                                initializeTopKinasesChart();
                            }}
                        }}, 100);
                    }})();
                </script>
            </div>
        </div>
        """
    else:
        sequence_html = """
        <div class="alert alert-warning">
            No sequence matches with kinase data available for this site.
        </div>
        """
    combined_html = "..."  # Implement similar to structural_html
    
    # Create tab content
    tab_content_html = f"""
    <div class="tab-content" id="siteAnalysisTabContent">
        <div class="tab-pane fade show active" id="structural" role="tabpanel">
            {structural_html}
        </div>
        <div class="tab-pane fade" id="sequence" role="tabpanel">
            {sequence_html}
        </div>
        <div class="tab-pane fade" id="combined" role="tabpanel">
            {combined_html}
        </div>
    </div>
    """
    
    # Combine all HTML
    final_html = f"""
    <div class="site-analysis" id="site-analysis">
        <h3>Kinase Analysis for Site {site_data['site_name']}</h3>
        
        {tabs_html}
        {tab_content_html}
        
        <div class="mt-4">
            <a href="/site/{uniprot_id}/{site_data['site_name']}" class="btn btn-primary">
                View Detailed Site Analysis
            </a>
        </div>
    </div>
    """
    
    return final_html

def prepare_heatmap_data(processed_matches: List[Dict], top_n_kinases: int = 15) -> Dict:
    """
    Prepare heatmap data from processed matches.
    
    Args:
        processed_matches: List of processed matches with kinase scores
        top_n_kinases: Number of top kinases to include
        
    Returns:
        Dictionary with heatmap data formatted for Plotly
    """
    if not processed_matches:
        return {'sites': [], 'kinases': [], 'scores': []}
    
    # Collect all kinase scores
    all_kinases = set()
    site_kinase_scores = {}
    
    for match in processed_matches:
        if 'kinase_scores' not in match:
            continue
            
        site_id = match.get('target_id', '')
        if not site_id:
            continue
            
        site_kinase_scores[site_id] = {}
        
        for kinase, score in match['kinase_scores'].items():
            all_kinases.add(kinase)
            site_kinase_scores[site_id][kinase] = score
    
    # Calculate median scores for each kinase across all sites
    kinase_median_scores = {}
    for kinase in all_kinases:
        scores = [
            site_scores.get(kinase, 0) 
            for site_scores in site_kinase_scores.values() 
            if kinase in site_scores
        ]
        if scores:
            kinase_median_scores[kinase] = np.median(scores)
    
    # Get top kinases by median score
    top_kinases = sorted(
        kinase_median_scores.keys(),
        key=lambda k: kinase_median_scores[k],
        reverse=True
    )[:top_n_kinases]
    
    # Prepare data for Plotly
    sites = list(site_kinase_scores.keys())
    kinases = top_kinases
    scores = []
    
    # Create scores array
    for site in sites:
        site_scores = site_kinase_scores[site]
        for kinase in kinases:
            scores.append({
                'site': site,
                'kinase': kinase,
                'score': site_scores.get(kinase, 0)
            })
    
    return {
        'sites': sites,
        'kinases': kinases,
        'scores': scores
    }

def prepare_top_kinases_data(processed_matches: List[Dict], top_n: int = 10) -> Dict:
    """
    Prepare top kinases data from processed matches.
    
    Args:
        processed_matches: List of processed matches with kinase scores
        top_n: Number of top kinases to include
        
    Returns:
        Dictionary with top kinases data formatted for Plotly
    """
    if not processed_matches:
        return {'kinases': [], 'scores': []}
    
    # Collect all kinase scores
    all_kinase_scores = {}
    
    for match in processed_matches:
        if 'kinase_scores' not in match:
            continue
            
        for kinase, score in match['kinase_scores'].items():
            if kinase not in all_kinase_scores:
                all_kinase_scores[kinase] = []
            all_kinase_scores[kinase].append(score)
    
    # Calculate median scores for each kinase
    kinase_median_scores = {}
    for kinase, scores in all_kinase_scores.items():
        if scores:
            kinase_median_scores[kinase] = np.median(scores)
    
    # Get top kinases by median score
    top_kinases_with_scores = sorted(
        kinase_median_scores.items(),
        key=lambda x: x[1],
        reverse=True
    )[:top_n]
    
    # Prepare data for Plotly
    kinases = [k for k, s in top_kinases_with_scores]
    scores = [s for k, s in top_kinases_with_scores]
    
    return {
        'kinases': kinases,
        'scores': scores
    }