from protein_explorer.db.db import get_phosphosites_batch
import json
import logging
import re

logger = logging.getLogger(__name__)

def create_sequence_network_visualization(protein_uniprot_id, phosphosites=None, sequence_matches=None):
    """
    Create a network visualization of phosphosites and their sequence similarity matches.
    Filters out Y (tyrosine) sites, ensures motif data is included, and adds known kinase info.
    
    Args:
        protein_uniprot_id: UniProt ID of the protein
        phosphosites: List of phosphosite dictionaries (optional)
        sequence_matches: Dictionary mapping site names to lists of match dictionaries (optional)
        
    Returns:
        HTML string with network visualization
    """
    import logging
    import json
    import re
    
    logger = logging.getLogger(__name__)
    
    logger.info(f"Creating sequence network visualization for {protein_uniprot_id}")
    
    # Check if we have any sequence matches to display
    has_matches = False
    if sequence_matches:
        for matches in sequence_matches.values():
            if matches:
                has_matches = True
                break
    
    if not has_matches:
        logger.info("No sequence matches found, returning empty visualization")
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
    
    # Process the phosphosites and matches to ensure consistent data types and node identification
    # Filter out Y (tyrosine) phosphosites and add known kinase information
    processed_matches = {}
    
    try:
        if sequence_matches:
            # Get all site IDs for known kinase lookup
            site_ids = []
            if phosphosites:
                for site in phosphosites:
                    if 'resno' in site:
                        site_ids.append(f"{protein_uniprot_id}_{site['resno']}")
            
            # Process the matches, filtering out Y sites
            for site_name, matches in sequence_matches.items():
                try:
                    # Skip Y sites
                    if site_name.startswith('Y'):
                        logger.info(f"Skipping tyrosine site: {site_name}")
                        continue
                    
                    processed_site_matches = []
                    
                    # Extract the site_type and site_number for consistent formatting
                    import re
                    site_match = re.match(r'([STY])(\d+)', site_name)
                    site_type = site_match.group(1) if site_match else 'S'  # Default to S if not matched
                    
                    # Skip if this is a Y site
                    if site_type == 'Y':
                        continue
                    
                    site_number = site_match.group(2) if site_match else site_name
                    
                    # Create a standardized site name for this site
                    std_site_name = site_name if site_match else f"{site_type}{site_name}"
                    
                    # Get known kinase information for this site if available
                    site_id = f"{protein_uniprot_id}_{site_number}"
                    known_kinase = None
                    
                    # Extract known kinase from phosphosite data
                    if phosphosites:
                        for site in phosphosites:
                            if 'resno' in site and str(site['resno']) == str(site_number):
                                # Check for KINASE_1 through KINASE_5 fields
                                for i in range(1, 6):
                                    kinase_field = f"KINASE_{i}"
                                    if kinase_field in site and site[kinase_field] and site[kinase_field] != 'unlabeled':
                                        known_kinase = site[kinase_field]
                                        break
                                
                                # Also check for known_kinase field
                                if not known_kinase and 'known_kinase' in site and site['known_kinase']:
                                    known_kinase = site['known_kinase']
                                    
                                break
                    
                    for match in matches:
                        # Skip targets that are Y sites
                        target_site = match.get('target_site', '')
                        if target_site.startswith('Y'):
                            logger.info(f"Skipping tyrosine target: {target_site}")
                            continue
                            
                        # Create a processed copy of the match
                        processed_match = match.copy()
                        
                        # Add source site info for proper edge creation
                        processed_match['source_site'] = std_site_name
                        processed_match['source_uniprot'] = protein_uniprot_id
                        processed_match['source_known_kinase'] = known_kinase
                        
                        # Ensure target_site has consistent format (with residue type letter)
                        if 'target_site' in match:
                            target_site = str(match['target_site'])
                            # Check if target_site doesn't already start with S or T (skip Y)
                            if not re.match(r'^[ST]', target_site):
                                # Extract site type from other fields if available
                                target_site_type = 'S'  # Default to S
                                
                                # Try to get site_type, with safety checks
                                if 'site_type' in match and match['site_type'] is not None:
                                    if isinstance(match['site_type'], str) and match['site_type'] in 'ST':
                                        target_site_type = match['site_type']
                                    
                                # Use the determined type
                                processed_match['target_site'] = f"{target_site_type}{target_site}"
                                # Also store the standardized format for node identification
                                processed_match['std_target_site'] = f"{target_site_type}{target_site}"
                        
                        # Make sure any existing motif is uppercase
                        if 'motif' in processed_match and processed_match['motif'] is not None:
                            if isinstance(processed_match['motif'], str):
                                processed_match['motif'] = processed_match['motif'].upper()
                        
                        # Add kinase information for target site
                        if 'target_id' in processed_match:
                            target_id = processed_match['target_id']
                            
                            # Try to extract kinase info from target data if available
                            try:
                                from protein_explorer.db.db import get_phosphosite_data
                                target_data = get_phosphosite_data(target_id)
                                if target_data:
                                    # Check for KINASE fields
                                    for i in range(1, 6):
                                        kinase_field = f"KINASE_{i}"
                                        if kinase_field in target_data and target_data[kinase_field] and target_data[kinase_field] != 'unlabeled':
                                            processed_match['target_known_kinase'] = target_data[kinase_field]
                                            break
                                    
                                    # Also check for known_kinase field
                                    if 'known_kinase' in target_data and target_data['known_kinase'] and target_data['known_kinase'] != 'unlabeled':
                                        processed_match['target_known_kinase'] = target_data['known_kinase']
                            except Exception as e:
                                logger.warning(f"Error getting kinase info for {target_id}: {e}")
                                
                        processed_site_matches.append(processed_match)
                    
                    # Add to the processed matches dictionary
                    processed_matches[std_site_name] = processed_site_matches
                except Exception as site_err:
                    logger.error(f"Error processing site {site_name}: {site_err}")
        else:
            # If no matches, use empty dictionary
            processed_matches = {}
    except Exception as e:
        logger.error(f"Error processing sequence matches: {e}")
        import traceback
        logger.error(traceback.format_exc())
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
    try:
        matches_json = json.dumps(processed_matches)
    except Exception as json_err:
        logger.error(f"Error converting matches to JSON: {json_err}")
        matches_json = '{}'  # Default to empty object
        
    html += matches_json + """'></div>

    <!-- Store phosphosites data for better node creation -->
    <div id="phosphosites-data" style="display: none;" data-sites='"""
    
    # Convert phosphosites to JSON with explicit type conversion of is_known fields,
    # filtering out Y sites and adding known kinase information
    if phosphosites:
        try:
            processed_sites = []
            for site in phosphosites:
                # Skip Y sites
                if 'site' in site and site['site'] and site['site'][0] == 'Y':
                    logger.info(f"Skipping tyrosine site: {site['site']}")
                    continue
                
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
                
                processed_sites.append(site_copy)
            sites_json = json.dumps(processed_sites)
        except Exception as sites_err:
            logger.error(f"Error converting phosphosites to JSON: {sites_err}")
            sites_json = '[]'  # Default to empty array
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