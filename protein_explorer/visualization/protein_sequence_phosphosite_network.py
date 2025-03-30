from protein_explorer.db.db import get_phosphosites_batch
import json
import logging

logger = logging.getLogger(__name__)

def create_sequence_network_visualization(protein_uniprot_id, phosphosites=None, sequence_matches=None):
    """
    Create a network visualization of phosphosites and their sequence similarity matches.
    With enhanced data handling to ensure proper visualization of known sites.
    
    Args:
        protein_uniprot_id: UniProt ID of the protein
        phosphosites: List of phosphosite dictionaries (optional)
        sequence_matches: Dictionary mapping site names to lists of match dictionaries (optional)
        
    Returns:
        HTML string with network visualization
    """
    # Check if we have any sequence matches to display
    has_matches = False
    if sequence_matches:
        for matches in sequence_matches.values():
            if matches:
                has_matches = True
                break
    
    if not has_matches:
        return """
        <div class="card mb-4">
            <div class="card-header">
                <h5 class="mb-0">Phosphosite Sequence Similarity Network</h5>
            </div>
            <div class="card-body">
                <div class="alert alert-info">
                    No sequence similarity matches found for this protein. Try adjusting similarity thresholds or analyzing a different protein.
                </div>
            </div>
        </div>
        """
    
    # Process the phosphosites and matches to ensure consistent data types
    # This ensures the JS can accurately determine if a site is known
    processed_matches = {}
    
    if sequence_matches:
        # Get all site IDs for the protein nodes
        site_ids = []
        for site in phosphosites:
            if 'resno' in site:
                site_ids.append(f"{protein_uniprot_id}_{site['resno']}")
        
        # Get all target IDs from matches
        target_ids = []
        for site_matches in sequence_matches.values():
            for match in site_matches:
                if 'target_id' in match:
                    target_ids.append(match['target_id'])
        
        # Get phosphosite data for both protein sites and match targets
        all_site_ids = site_ids + target_ids
        logger.info(f"Getting phosphosite data for {len(all_site_ids)} sites")
        supp_data_dict = get_phosphosites_batch(all_site_ids)
        
        # Process the matches with enhanced data
        for site_name, matches in sequence_matches.items():
            processed_site_matches = []
            
            for match in matches:
                # Create an enhanced copy of the match
                processed_match = match.copy()
                
                # Add is_known_phosphosite explicitly to match target
                target_id = match.get('target_id')
                if target_id and target_id in supp_data_dict:
                    target_data = supp_data_dict[target_id]
                    
                    # Add supplementary data
                    if 'is_known_phosphosite' in target_data:
                        # Convert to float for consistent handling in JS
                        processed_match['is_known_phosphosite'] = float(target_data['is_known_phosphosite'])
                    
                    # Add motif if available
                    if 'SITE_+/-7_AA' in target_data and target_data['SITE_+/-7_AA'] and 'motif' not in processed_match:
                        processed_match['motif'] = target_data['SITE_+/-7_AA']
                    
                    # Add other useful fields
                    for field in ['site_plddt', 'mean_plddt', 'nearby_count', 'surface_accessibility']:
                        if field in target_data and target_data[field] is not None:
                            processed_match[field] = target_data[field]
                
                # Add the processed match
                processed_site_matches.append(processed_match)
            
            # Add to the processed matches dictionary
            processed_matches[site_name] = processed_site_matches
    else:
        # If no matches, use empty dictionary
        processed_matches = {}
    
    # Create network visualization HTML
    html = """
    <div class="card mb-4">
        <div class="card-header">
            <h5 class="mb-0">Phosphosite Sequence Similarity Network</h5>
        </div>
        <div class="card-body">
            <div id="similarity-filter-container" class="mb-3">
                <label for="similarity-filter" class="form-label">
                    Similarity Threshold: <span id="similarity-value">0.6</span>
                </label>
                <input 
                    type="range" 
                    class="form-range" 
                    id="similarity-filter" 
                    min="0.4" 
                    max="0.9" 
                    step="0.05" 
                    value="0.6"
                    oninput="document.getElementById('similarity-value').textContent = this.value; updateSequenceNetworkFilter();"
                >
                <div class="d-flex justify-content-between">
                    <small>0.4 (Less similar)</small>
                    <small>0.9 (Very similar)</small>
                </div>
            </div>
            
            <div class="mb-3">
                <div class="d-flex align-items-center flex-wrap">
                    <div class="d-flex align-items-center me-4 mb-2">
                        <div style="width: 16px; height: 16px; background-color: #4CAF50; border-radius: 50%; margin-right: 6px;"></div>
                        <span class="small">Known protein sites</span>
                    </div>
                    <div class="d-flex align-items-center me-4 mb-2">
                        <div style="width: 16px; height: 16px; background-color: #FF9800; border-radius: 50%; margin-right: 6px;"></div>
                        <span class="small">Unknown protein sites</span>
                    </div>
                    <div class="d-flex align-items-center mb-2">
                        <div style="width: 16px; height: 16px; background-color: #9C27B0; border-radius: 50%; margin-right: 6px;"></div>
                        <span class="small">Sequence-similar sites</span>
                    </div>
                </div>
            </div>
            
            <div id="sequence-network-container" style="height: 500px; width: 100%; position: relative; border: 1px solid #ddd; border-radius: 5px;"></div>
            
            <p class="text-muted mt-3 mb-0">
                <small>
                    This network shows the sequence relationships between phosphosites in this protein and
                    similar sites in other proteins. Edges represent sequence similarity above
                    the threshold. Hover over nodes for details and click to view the site page.
                </small>
            </p>
        </div>
    </div>

    <!-- Store match data for the visualization -->
    <div id="sequence-match-data" style="display: none;" data-matches='"""
    
    # Convert processed matches to JSON with proper handling of is_known field
    matches_json = json.dumps(processed_matches)
    html += matches_json + """'></div>

    <!-- Store phosphosites data for better node creation -->
    <div id="phosphosites-data" style="display: none;" data-sites='"""
    
    # Convert phosphosites to JSON with explicit type conversion of is_known fields
    if phosphosites:
        processed_sites = []
        for site in phosphosites:
            site_copy = site.copy()
            
            # Ensure is_known and is_known_phosphosite are properly set
            if 'is_known' in site:
                site_copy['is_known'] = bool(site['is_known'])
            
            if 'is_known_phosphosite' in site:
                # Convert to float for consistent handling
                site_copy['is_known_phosphosite'] = float(site['is_known_phosphosite'])
            elif 'is_known' in site:
                # Derive from is_known if is_known_phosphosite not present
                site_copy['is_known_phosphosite'] = 1.0 if site['is_known'] else 0.0
                
            processed_sites.append(site_copy)
        sites_json = json.dumps(processed_sites)
    else:
        sites_json = '[]'
    
    html += sites_json + """'></div>

    <!-- Inline script to ensure the network visualization works -->
    <script>
    document.addEventListener('DOMContentLoaded', function() {
        console.log('DOM loaded, initializing sequence network for """ + protein_uniprot_id + """');
        // Make sure the visualization function is available
        if (typeof sequenceNetworkVisualization === 'function') {
            // Call the main visualization function
            sequenceNetworkVisualization('""" + protein_uniprot_id + """');
        } else {
            console.error('Sequence network visualization function not found!');
            // Try to load it dynamically
            var script = document.createElement('script');
            script.src = "/static/js/protein-sequence-network-visualization.js";
            script.onload = function() {
                if (typeof sequenceNetworkVisualization === 'function') {
                    sequenceNetworkVisualization('""" + protein_uniprot_id + """');
                } else {
                    console.error('Failed to load sequence network visualization!');
                }
            };
            document.head.appendChild(script);
        }
    });
    </script>
    """
    
    return html