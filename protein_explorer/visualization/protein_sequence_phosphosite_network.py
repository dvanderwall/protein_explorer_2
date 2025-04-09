from protein_explorer.db.db import get_phosphosites_batch
import json
import logging
import re
import traceback

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
    
    # Function to sanitize string values in dictionaries
    def sanitize_data(data):
        """Recursively clean data structure to ensure serialization safety."""
        if isinstance(data, dict):
            return {k: sanitize_data(v) for k, v in data.items()}
        elif isinstance(data, list):
            return [sanitize_data(item) for item in data]
        elif isinstance(data, str):
            # Just remove the characters; do not replace with others
            cleaned = data.replace("'", "").replace('"', "").replace("\n", "").replace("\r", "")
            return cleaned
        else:
            return data
    
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
    
    # Sanitize sequence matches data before processing
    processed_matches = sanitize_data(sequence_matches) if sequence_matches else {}
    
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
        # Use json.dumps with ensure_ascii=False to preserve special characters
        matches_json = json.dumps(processed_matches, ensure_ascii=False)
        
        # Additional safety check to catch any problematic content
        if '"DISEASE": "B cell lymphoma; non-Hodgkin' in matches_json and not '"DISEASE": "B cell lymphoma; non-Hodgkin"' in matches_json:
            matches_json = matches_json.replace('"DISEASE": "B cell lymphoma; non-Hodgkin', 
                                              '"DISEASE": "B cell lymphoma; non-Hodgkin lymphoma"')
            logger.warning("Fixed unterminated string in DISEASE field")
        
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
                if 'siteType' in site and site['siteType'] and site['siteType'][0] == 'Y':
                    continue
                
                # Create a sanitized copy of the site data
                site_copy = sanitize_data(site.copy())
                
                # Ensure is_known and is_known_phosphosite are properly set
                if 'is_known' in site:
                    site_copy['is_known'] = bool(site['is_known'])
                
                if 'is_known_phosphosite' in site:
                    # Convert to float for consistent handling
                    try:
                        site_copy['is_known_phosphosite'] = float(site['is_known_phosphosite'])
                    except (ValueError, TypeError):
                        site_copy['is_known_phosphosite'] = 0.0
                elif 'is_known' in site:
                    # Derive from is_known if is_known_phosphosite not present
                    site_copy['is_known_phosphosite'] = 1.0 if site['is_known'] else 0.0
                
                # Use extract_kinase_info function for efficient kinase extraction
                kinase_info = extract_kinase_info(site)
                if kinase_info:
                    site_copy['knownKinase'] = kinase_info
                
                processed_sites.append(site_copy)
            
            # Use ensure_ascii=False for better handling of special characters
            sites_json = json.dumps(processed_sites, ensure_ascii=False)
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
    
    # Print first 10 key-value pairs
    #print("PROCESSED SEQUENCE SITES")
    #print(processed_sites[:3] if processed_sites else [])
    #print("PROCESSED SEQUENCE MATCHES")
    #print(list(processed_matches.items())[:2] if processed_matches else [])

    #inspect_data_for_visualization(protein_uniprot_id, processed_sites, processed_matches)

    #print("THE HTML")
    #print(html)
    #with open("test_html_output.html", "w", encoding="utf-8") as f:
    #    f.write(html)
    return html

def extract_kinase_info(site_data):
    """
    Extract all known kinase information from site data.
    Handles multiple data formats and sources with robust error handling.
    
    Args:
        site_data: Dictionary with site information
        
    Returns:
        Comma-separated string of unique kinases or None
    """
    if not site_data or not isinstance(site_data, dict):
        return None
        
    kinases = []
    
    try:
        # Check KINASE_1 through KINASE_5 fields
        for i in range(1, 6):
            kinase_field = f"KINASE_{i}"
            if kinase_field in site_data and site_data[kinase_field] and site_data[kinase_field] != 'unlabeled':
                try:
                    kinases.append(str(site_data[kinase_field]))
                except (TypeError, ValueError):
                    # Skip if can't convert to string
                    pass
        
        # Check known_kinase field with careful string handling
        if 'known_kinase' in site_data and site_data['known_kinase']:
            if isinstance(site_data['known_kinase'], str):
                if ',' in site_data['known_kinase']:
                    # Multiple kinases comma-separated
                    for k in site_data['known_kinase'].split(','):
                        if k.strip() and k.strip() != 'unlabeled':
                            kinases.append(k.strip())
                else:
                    # Single kinase as string
                    if site_data['known_kinase'].strip() != 'unlabeled':
                        kinases.append(site_data['known_kinase'].strip())
            else:
                # Try to convert to string if it's not already
                try:
                    kinase_str = str(site_data['known_kinase'])
                    if kinase_str and kinase_str != 'unlabeled':
                        kinases.append(kinase_str)
                except (TypeError, ValueError):
                    # Skip if can't convert to string
                    pass
        
        # Check target_known_kinase field
        if 'target_known_kinase' in site_data and site_data['target_known_kinase']:
            try:
                if isinstance(site_data['target_known_kinase'], str):
                    if ',' in site_data['target_known_kinase']:
                        # Multiple kinases comma-separated
                        for k in site_data['target_known_kinase'].split(','):
                            if k.strip() and k.strip() != 'unlabeled':
                                kinases.append(k.strip())
                    else:
                        # Single kinase as string
                        if site_data['target_known_kinase'].strip() != 'unlabeled':
                            kinases.append(site_data['target_known_kinase'].strip())
                else:
                    # Try to convert to string if it's not already
                    kinase_str = str(site_data['target_known_kinase'])
                    if kinase_str and kinase_str != 'unlabeled':
                        kinases.append(kinase_str)
            except (TypeError, ValueError):
                # Skip if can't convert to string
                pass
        
        # Check for kinases array
        if 'kinases' in site_data and isinstance(site_data['kinases'], list):
            for kinase in site_data['kinases']:
                if kinase and kinase != 'unlabeled':
                    try:
                        kinases.append(str(kinase))
                    except (TypeError, ValueError):
                        # Skip if can't convert to string
                        pass
        
        # Check top_kinases field for structured kinase data
        if 'top_kinases' in site_data and isinstance(site_data['top_kinases'], list) and site_data['top_kinases']:
            try:
                top_kinase = site_data['top_kinases'][0]
                if isinstance(top_kinase, dict):
                    kinase_name = top_kinase.get('kinase') or top_kinase.get('name')
                    if kinase_name and kinase_name != 'unlabeled':
                        kinases.append(str(kinase_name))
            except (IndexError, TypeError, ValueError):
                # Skip if there's an error accessing or converting
                pass
        
    except Exception as e:
        # Catch-all for any unexpected errors
        logger = logging.getLogger(__name__)
        logger.error(f"Error extracting kinase info: {e}")
        # Try to continue with any kinases we were able to extract
    
    # Remove duplicates while preserving order
    if kinases:
        unique_kinases = []
        seen = set()
        for kinase in kinases:
            if kinase and kinase not in seen and kinase.lower() != 'unlabeled':
                seen.add(kinase)
                unique_kinases.append(kinase)
        return ', '.join(unique_kinases)
    
    return None


def inspect_data_for_visualization(protein_uniprot_id, phosphosites, sequence_matches):
    """
    Utility function to deeply inspect data structure before visualization.
    Logs detailed information about data types, sizes, and sample values.
    
    Args:
        protein_uniprot_id: UniProt ID of the protein
        phosphosites: List of phosphosite dictionaries (from route)
        sequence_matches: Dictionary mapping site names to lists of match dictionaries (from route)
    """
    logger = logging.getLogger(__name__)
    logger.info("=" * 80)
    logger.info(f"DATA INSPECTION FOR {protein_uniprot_id}")
    logger.info("=" * 80)
    
    # Check protein ID
    logger.info(f"Protein UniProt ID: {protein_uniprot_id} (type: {type(protein_uniprot_id).__name__})")
    
    # Inspect phosphosites data
    logger.info("-" * 40)
    logger.info("PHOSPHOSITES INSPECTION")
    logger.info("-" * 40)
    
    if phosphosites is None:
        logger.info("Phosphosites is None")
    else:
        # Check type
        logger.info(f"Phosphosites type: {type(phosphosites).__name__}")
        
        # Check if iterable
        is_iterable = hasattr(phosphosites, '__iter__') and not isinstance(phosphosites, str)
        logger.info(f"Is iterable: {is_iterable}")
        
        if is_iterable:
            try:
                logger.info(f"Length: {len(phosphosites)}")
                
                # Sample values if not empty
                if len(phosphosites) > 0:
                    logger.info("First item type: " + type(phosphosites[0]).__name__)
                    
                    # Check if first item is dict
                    if isinstance(phosphosites[0], dict):
                        logger.info(f"First item keys: {list(phosphosites[0].keys())}")
                        
                        # Check some critical fields
                        for field in ['site', 'is_known', 'is_known_phosphosite', 'motif', 'known_kinase']:
                            if field in phosphosites[0]:
                                value = phosphosites[0][field]
                                logger.info(f"Field '{field}': {value} (type: {type(value).__name__})")
                    
                    # Count site types (S, T, Y)
                    site_types = {}
                    for p in phosphosites:
                        if isinstance(p, dict) and 'site' in p and p['site']:
                            site_type = p['site'][0] if isinstance(p['site'], str) else 'unknown'
                            site_types[site_type] = site_types.get(site_type, 0) + 1
                    
                    logger.info(f"Site type counts: {site_types}")
            except Exception as e:
                logger.error(f"Error examining phosphosites: {e}")
    
    # Inspect sequence matches data
    logger.info("-" * 40)
    logger.info("SEQUENCE MATCHES INSPECTION")
    logger.info("-" * 40)
    
    if sequence_matches is None:
        logger.info("Sequence matches is None")
    else:
        # Check type
        logger.info(f"Sequence matches type: {type(sequence_matches).__name__}")
        
        if isinstance(sequence_matches, dict):
            logger.info(f"Number of keys: {len(sequence_matches)}")
            
            if len(sequence_matches) > 0:
                # Sample first key and its matches
                sample_key = next(iter(sequence_matches.keys()))
                logger.info(f"Sample key: {sample_key} (type: {type(sample_key).__name__})")
                
                matches = sequence_matches[sample_key]
                logger.info(f"Matches for key type: {type(matches).__name__}")
                
                if isinstance(matches, list):
                    logger.info(f"Number of matches for sample key: {len(matches)}")
                    
                    if len(matches) > 0:
                        logger.info(f"First match type: {type(matches[0]).__name__}")
                        
                        if isinstance(matches[0], dict):
                            logger.info(f"First match keys: {list(matches[0].keys())}")
                            
                            # Check critical fields
                            for field in ['target_id', 'similarity', 'is_known_phosphosite', 'motif', 'known_kinase']:
                                if field in matches[0]:
                                    value = matches[0][field]
                                    logger.info(f"Field '{field}': {value} (type: {type(value).__name__})")
                
                # Check for common problems
                problem_keys = []
                for key, value in sequence_matches.items():
                    if not isinstance(value, list):
                        problem_keys.append(f"{key} (type: {type(value).__name__})")
                    elif len(value) > 0 and not isinstance(value[0], dict):
                        problem_keys.append(f"{key} has non-dict items (type: {type(value[0]).__name__})")
                
                if problem_keys:
                    logger.warning(f"Problem keys found: {problem_keys}")
                
                # Count similarity values
                similarity_ranges = {
                    "0.0-0.3": 0,
                    "0.3-0.6": 0,
                    "0.6-0.9": 0,
                    "0.9-1.0": 0,
                    "invalid": 0
                }
                
                for matches_list in sequence_matches.values():
                    if isinstance(matches_list, list):
                        for match in matches_list:
                            if isinstance(match, dict) and 'similarity' in match:
                                try:
                                    sim = float(match['similarity'])
                                    if sim < 0.3:
                                        similarity_ranges["0.0-0.3"] += 1
                                    elif sim < 0.6:
                                        similarity_ranges["0.3-0.6"] += 1
                                    elif sim < 0.9:
                                        similarity_ranges["0.6-0.9"] += 1
                                    else:
                                        similarity_ranges["0.9-1.0"] += 1
                                except (ValueError, TypeError):
                                    similarity_ranges["invalid"] += 1
                
                logger.info(f"Similarity value distribution: {similarity_ranges}")
    
    logger.info("=" * 80)
    logger.info("END DATA INSPECTION")
    logger.info("=" * 80)