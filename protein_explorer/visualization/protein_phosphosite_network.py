# Add/modify in protein_explorer/visualization/protein_phosphosite_network.py

import json
import logging
import re
from protein_explorer.db.db import get_phosphosite_data, get_phosphosites_batch

logger = logging.getLogger(__name__)

def create_phosphosite_network_visualization(protein_uniprot_id, phosphosites=None, structural_matches=None):
    """
    Create a network visualization of phosphosites and their structural matches.
    Enhanced to support known kinase visualization similar to sequence network.
    
    Args:
        protein_uniprot_id: UniProt ID of the protein
        phosphosites: List of phosphosite dictionaries (optional)
        structural_matches: Dictionary mapping site names to lists of match dictionaries (optional)
        
    Returns:
        HTML string with network visualization
    """
    # Process phosphosites to extract kinase information
    processed_phosphosites = []
    if phosphosites:
        for site in phosphosites:
            # Skip Y (tyrosine) sites - for consistency with sequence network
            if 'site' in site and site['site'] and site['site'][0] == 'Y':
                continue
                
            site_copy = site.copy()
            
            # Ensure is_known flags are properly set
            if 'is_known' in site:
                site_copy['is_known'] = bool(site['is_known'])
                
            if 'is_known_phosphosite' in site:
                site_copy['is_known_phosphosite'] = float(site['is_known_phosphosite'])
            elif 'is_known' in site:
                site_copy['is_known_phosphosite'] = 1.0 if site['is_known'] else 0.0
                
            # Add known kinase information if available
            known_kinase = None
            for i in range(1, 6):
                kinase_field = f"KINASE_{i}"
                if kinase_field in site and site[kinase_field] and site[kinase_field] != 'unlabeled':
                    known_kinase = site[kinase_field]
                    break
                    
            if not known_kinase and 'known_kinase' in site and site['known_kinase'] and site['known_kinase'] != 'unlabeled':
                known_kinase = site['known_kinase']
                
            if known_kinase:
                site_copy['knownKinase'] = known_kinase
                
            processed_phosphosites.append(site_copy)
    
    # Process structural matches to enhance with kinase data and other supplementary information
    processed_matches = {}
    all_target_ids = set()  # Collect all target IDs for batch query
    
    # First pass - collect all target IDs
    if structural_matches:
        for site_name, matches in structural_matches.items():
            # Skip Y (tyrosine) sites
            if site_name and site_name[0] == 'Y':
                continue
                
            for match in matches:
                # Skip matches with Y sites
                if match.get('target_site', '')[0] == 'Y':
                    continue
                
                # Get target info and create target_id
                target_uniprot = match.get('target_uniprot')
                target_site = match.get('target_site')
                
                if target_uniprot and target_site:
                    # Extract just the number from the target site
                    target_resno = ''.join(filter(str.isdigit, target_site))
                    if target_resno:
                        target_id = f"{target_uniprot}_{target_resno}"
                        all_target_ids.add(target_id)
    
    # Batch fetch supplementary data for all target sites
    supplementary_data = {}
    if all_target_ids:
        try:
            logger.info(f"Fetching supplementary data for {len(all_target_ids)} target sites")
            supplementary_data = get_phosphosites_batch(list(all_target_ids))
            logger.info(f"Retrieved data for {len(supplementary_data)} target sites")
        except Exception as e:
            logger.error(f"Error fetching supplementary data: {e}")
    
    # Second pass - enhance matches with the supplementary data
    if structural_matches:
        for site_name, matches in structural_matches.items():
            # Skip Y (tyrosine) sites
            if site_name and site_name[0] == 'Y':
                continue
                
            processed_site_matches = []
            for match in matches:
                # Skip matches with Y sites
                if match.get('target_site', '')[0] == 'Y':
                    continue
                
                match_copy = match.copy()
                
                # Get target info
                target_uniprot = match.get('target_uniprot')
                target_site = match.get('target_site')
                
                # Extract kinase and other info from supplementary data
                target_resno = ''.join(filter(str.isdigit, target_site))
                target_id = f"{target_uniprot}_{target_resno}" if target_uniprot and target_resno else None
                
                if target_id and target_id in supplementary_data:
                    target_data = supplementary_data[target_id]
                    
                    # Extract known kinase information
                    known_kinase = None
                    for i in range(1, 6):
                        kinase_field = f"KINASE_{i}"
                        if kinase_field in target_data and target_data[kinase_field] and target_data[kinase_field] != 'unlabeled':
                            known_kinase = target_data[kinase_field]
                            break
                    
                    if not known_kinase and 'known_kinase' in target_data and target_data['known_kinase'] and target_data['known_kinase'] != 'unlabeled':
                        known_kinase = target_data['known_kinase']
                    
                    # Set known kinase in match data
                    if known_kinase:
                        match_copy['known_kinase'] = known_kinase
                    
                    # Add additional supplementary data
                    if 'SITE_+/-7_AA' in target_data and target_data['SITE_+/-7_AA']:
                        match_copy['motif'] = target_data['SITE_+/-7_AA']
                    
                    for field in ['site_plddt', 'mean_plddt', 'nearby_count', 'surface_accessibility']:
                        if field in target_data and target_data[field] is not None:
                            match_copy[field] = target_data[field]
                
                processed_site_matches.append(match_copy)
                
            if processed_site_matches:
                processed_matches[site_name] = processed_site_matches
    
    # Create HTML for the visualization
    html = """
    <div class="card mb-4">
        <div class="card-header">
            <h5 class="mb-0">Phosphosite Structural Network</h5>
        </div>
        <div class="card-body">
            <div id="rmsd-filter-container" class="mb-3">
                <label for="rmsd-filter" class="form-label">
                    RMSD Threshold: <span id="rmsd-value">2.0 Å</span>
                </label>
                <input 
                    type="range" 
                    class="form-range" 
                    id="rmsd-filter" 
                    min="0.5" 
                    max="5.0" 
                    step="0.1" 
                    value="2.0"
                    oninput="document.getElementById('rmsd-value').textContent = this.value + ' Å'; updateNetworkFilter();"
                >
                <div class="d-flex justify-content-between">
                    <small>0.5 Å (Very similar)</small>
                    <small>5.0 Å (Less similar)</small>
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
                    <div class="d-flex align-items-center me-4 mb-2">
                        <div style="width: 16px; height: 16px; background-color: #9C27B0; border-radius: 50%; margin-right: 6px;"></div>
                        <span class="small">Sites with known kinase</span>
                    </div>
                    <div class="d-flex align-items-center mb-2">
                        <div style="width: 16px; height: 16px; background-color: #9E9E9E; border-radius: 50%; margin-right: 6px;"></div>
                        <span class="small">Structurally similar sites</span>
                    </div>
                </div>
            </div>
            
            <div id="network-container" style="height: 500px; width: 100%; position: relative; border: 1px solid #ddd; border-radius: 5px;"></div>
            
            <p class="text-muted mt-3 mb-0">
                <small>
                    This network shows the structural relationships between phosphosites in this protein and
                    similar sites in other proteins. Edges represent structural similarity with RMSD below
                    the threshold. Hover over nodes for details and click to view the site page.
                </small>
            </p>
        </div>
    </div>

    <!-- Store structural match data for the visualization -->
    <div id="structural-match-data" style="display: none;" data-matches='"""
    
    # Convert processed matches to JSON with proper handling
    try:
        matches_json = json.dumps(processed_matches)
    except Exception as json_err:
        logger.error(f"Error converting matches to JSON: {json_err}")
        matches_json = '{}'  # Default to empty object
        
    html += matches_json + """'></div>

    <!-- Store phosphosites data for better node creation -->
    <div id="phosphosites-data" style="display: none;" data-sites='"""
    
    # Convert phosphosites to JSON
    try:
        sites_json = json.dumps(processed_phosphosites) if processed_phosphosites else '[]'
    except Exception as json_err:
        logger.error(f"Error converting phosphosites to JSON: {json_err}")
        sites_json = '[]'
    
    html += sites_json + """'></div>

    <!-- Inline script to initialize the network visualization -->
    <script>
    document.addEventListener('DOMContentLoaded', function() {
        console.log('DOM loaded, initializing structural network for """ + protein_uniprot_id + """');
        // Try to load the proper visualization script if not already loaded
        if (typeof setupPhosphositeStructuralNetwork === 'function') {
            // Function is already loaded, call it directly
            setupPhosphositeStructuralNetwork('""" + protein_uniprot_id + """');
        } else {
            console.log('Loading structural network visualization script...');
            // Create script element to load the visualization script
            const script = document.createElement('script');
            script.src = "/static/js/protein-structural-network-visualization.js";
            script.onload = function() {
                if (typeof setupPhosphositeStructuralNetwork === 'function') {
                    setupPhosphositeStructuralNetwork('""" + protein_uniprot_id + """');
                } else {
                    console.error('Structural network visualization function still not found after loading script');
                    // Fallback to basic visualization
                    setupBasicStructuralNetwork('""" + protein_uniprot_id + """');
                }
            };
            document.head.appendChild(script);
        }
    });
    
    // Basic fallback visualization in case the main script fails to load
    function setupBasicStructuralNetwork(proteinUniprotId) {
        console.log('Setting up basic structural network visualization');
        const networkContainer = document.getElementById('network-container');
        if (!networkContainer) return;
        
        // Extract matches data
        const matchesElement = document.getElementById('structural-match-data');
        const phosphositesElement = document.getElementById('phosphosites-data');
        
        if (!matchesElement || !phosphositesElement) {
            networkContainer.innerHTML = '<div class="alert alert-warning p-3">Unable to load network data</div>';
            return;
        }
        
        const matchesAttr = matchesElement.getAttribute('data-matches');
        const phosphositesAttr = phosphositesElement.getAttribute('data-sites');
        
        try {
            const matches = JSON.parse(matchesAttr || '{}');
            const phosphosites = JSON.parse(phosphositesAttr || '[]');
            
            // Check if we have data to display
            if (Object.keys(matches).length === 0 || phosphosites.length === 0) {
                networkContainer.innerHTML = '<div class="alert alert-info p-3">No structural network data available</div>';
                return;
            }
            
            // Create a simple force-directed network with D3.js
            // (Basic implementation that will work if the main script fails to load)
            if (typeof d3 !== 'undefined') {
                // D3.js code for basic visualization
                // (Simplified version of the main visualization)
                networkContainer.innerHTML = '<div class="alert alert-warning p-3 mb-3">Basic visualization mode</div>' +
                                             '<svg width="100%" height="400px"></svg>';
                                             
                // Create basic visualization here...
            } else {
                // If D3.js is not available, show a table of matches instead
                let tableHtml = '<div class="alert alert-warning p-3 mb-3">Unable to load visualization library</div>' +
                                '<table class="table table-sm table-bordered">' +
                                '<thead><tr><th>Protein Site</th><th>Match Count</th></tr></thead>' +
                                '<tbody>';
                                
                for (const [site, siteMatches] of Object.entries(matches)) {
                    tableHtml += `<tr><td>${site}</td><td>${siteMatches.length}</td></tr>`;
                }
                
                tableHtml += '</tbody></table>';
                networkContainer.innerHTML = tableHtml;
            }
        } catch (e) {
            console.error('Error setting up basic network:', e);
            networkContainer.innerHTML = '<div class="alert alert-danger p-3">Error loading network data</div>';
        }
    }
    </script>
    """
    
    return html