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
    # Using the optimized batch processing approach
    processed_matches = {}
    
    try:
        if sequence_matches:
            # Build a lookup map for phosphosite data - speeds up repeated lookups
            site_data_map = {}
            if phosphosites:
                for site in phosphosites:
                    if 'resno' in site:
                        # Use both resno and site name as keys for easier lookup
                        site_data_map[str(site['resno'])] = site
                        if 'site' in site:
                            site_data_map[site['site']] = site
            
            # Pre-fetch kinase information for all target sites in one batch
            target_ids = set()
            for matches in sequence_matches.values():
                for match in matches:
                    if 'target_id' in match:
                        target_ids.add(match['target_id'])
            
            # Batch retrieve target site data
            target_data_map = {}
            if target_ids:
                try:
                    target_data_map = get_phosphosites_batch(list(target_ids))
                    logger.info(f"Retrieved data for {len(target_data_map)} target sites")
                except Exception as e:
                    logger.warning(f"Error batch fetching target data: {e}")
            
            # Process each site's matches
            for site_name, matches in sequence_matches.items():
                try:
                    # Skip Y sites
                    if site_name.startswith('Y'):
                        logger.info(f"Skipping tyrosine site: {site_name}")
                        continue
                    
                    # Extract site type and number
                    site_match = re.match(r'([STY])(\d+)', site_name)
                    if not site_match or site_match.group(1) == 'Y':
                        continue
                        
                    site_type = site_match.group(1)
                    site_number = site_match.group(2)
                    std_site_name = site_name
                    site_id = f"{protein_uniprot_id}_{site_number}"
                    
                    # Get site data from our lookup map
                    site_data = site_data_map.get(site_number, site_data_map.get(site_name))
                    
                    # Get kinase information
                    source_known_kinase = extract_kinase_info(site_data)
                    
                    # Process matches for this site
                    processed_site_matches = []
                    for match in matches:
                        # Skip Y targets
                        target_site = match.get('target_site', '')
                        if target_site.startswith('Y'):
                            logger.info(f"Skipping tyrosine target: {target_site}")
                            continue
                        
                        # Create processed match with a copy to avoid modifying original
                        processed_match = match.copy()
                        
                        # Add source site info
                        processed_match['source_site'] = std_site_name
                        processed_match['source_uniprot'] = protein_uniprot_id
                        processed_match['source_known_kinase'] = source_known_kinase
                        
                        # Fix target site format if needed
                        if 'target_site' in processed_match:
                            target_site = str(processed_match['target_site'])
                            if not re.match(r'^[ST]', target_site):
                                target_site_type = (processed_match.get('site_type') 
                                                  if processed_match.get('site_type') in 'ST' 
                                                  else 'S')
                                processed_match['target_site'] = f"{target_site_type}{target_site}"
                                processed_match['std_target_site'] = f"{target_site_type}{target_site}"
                        
                        # Ensure motif is uppercase
                        if 'motif' in processed_match and processed_match['motif']:
                            processed_match['motif'] = str(processed_match['motif']).upper()
                        
                        # Add target kinase info from our pre-fetched data
                        if 'target_id' in processed_match:
                            target_id = processed_match['target_id']
                            target_data = target_data_map.get(target_id)
                            if target_data:
                                target_known_kinase = extract_kinase_info(target_data)
                                if target_known_kinase:
                                    processed_match['target_known_kinase'] = target_known_kinase
                        
                        processed_site_matches.append(processed_match)
                    
                    # Only add sites that have matches after filtering
                    if processed_site_matches:
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
    

    print('JSONIFIED PROCESSED MATCHES')
    print(processed_matches)

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
                    min="0.3" 
                    max="0.99" 
                    step="0.05" 
                    value="0.01"
                    oninput="document.getElementById('similarity-value').textContent = this.value; updateSequenceNetworkFilter();"
                >
                <div class="d-flex justify-content-between">
                    <small>0.3 (Less similar)</small>
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
                    <div class="d-flex align-items-center me-4 mb-2">
                        <div style="width: 16px; height: 16px; background-color: #9C27B0; border-radius: 50%; margin-right: 6px;"></div>
                        <span class="small">Sites with known kinase</span>
                    </div>
                    <div class="d-flex align-items-center mb-2">
                        <div style="width: 16px; height: 16px; background-color: #E91E63; border-radius: 50%; margin-right: 6px;"></div>
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
    
    # Convert processed matches to JSON with proper handling
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
            # Use an optimized approach for processing phosphosites
            processed_sites = []
            for site in phosphosites:
                # Skip Y sites
                if 'site' in site and site['site'] and site['site'][0] == 'Y':
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
                
                # Use extract_kinase_info function for efficient kinase extraction
                kinase_info = extract_kinase_info(site)
                if kinase_info:
                    site_copy['knownKinase'] = kinase_info
                
                processed_sites.append(site_copy)
            print("JSONIFIED PHOSPHOSITES")
            print(processed_sites)
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

def extract_kinase_info(site_data):
    """
    Extract all known kinase information from site data.
    
    Args:
        site_data: Dictionary with site information
        
    Returns:
        Comma-separated string of unique kinases or None
    """
    if not site_data:
        return None
        
    kinases = []
    
    # Check KINASE_1 through KINASE_5 fields
    for i in range(1, 6):
        kinase_field = f"KINASE_{i}"
        if kinase_field in site_data and site_data[kinase_field] and site_data[kinase_field] != 'unlabeled':
            kinases.append(site_data[kinase_field])
    
    # Check known_kinase field
    if 'known_kinase' in site_data and site_data['known_kinase']:
        if isinstance(site_data['known_kinase'], str) and ',' in site_data['known_kinase']:
            kinases.extend([k.strip() for k in site_data['known_kinase'].split(',')])
        else:
            kinases.append(site_data['known_kinase'])
    
    # Remove duplicates while preserving order
    if kinases:
        unique_kinases = []
        seen = set()
        for kinase in kinases:
            if kinase and kinase not in seen:
                seen.add(kinase)
                unique_kinases.append(kinase)
        return ', '.join(unique_kinases)
    
    return None