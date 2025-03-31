"""
Functions for generating enhanced HTML tables for phosphosite visualization
with improved metric calculations.
"""
from protein_explorer.analysis.phospho_analyzer_2 import get_phosphosite_data
# Modification to enhance_phosphosite_table in protein_explorer/analysis/enhanced_table.py

# Modify enhance_phosphosite_table in protein_explorer/analysis/enhanced_table.py
# to handle NaN values properly

# In protein_explorer/analysis/enhanced_table.py
# Only the relevant changes to add the kinase column

def enhance_phosphosite_table(phosphosites, protein_uniprot_id):
    """
    Add data attributes to the phosphosite table HTML for better visualization.
    
    Args:
        phosphosites: List of phosphosite dictionaries
        protein_uniprot_id: UniProt ID of the protein
        
    Returns:
        HTML string with the enhanced phosphosite table
    """
    import numpy as np
    import math
    
    if not phosphosites:
        return "<div class='alert alert-warning'>No phosphosite data available.</div>"
    
    # Calculate additional metrics for each site if not already present
    for site in phosphosites:
        # Convert motif to uppercase if present
        if 'motif' in site and site['motif']:
            site['motif'] = site['motif'].upper()
            
        if 'motif' in site:
            # Calculate metrics based on motif if available
            if 'acidicPercentage' not in site:
                site['acidicPercentage'] = calculate_acidic_percentage(site['motif'])
            
            if 'basicPercentage' not in site:
                site['basicPercentage'] = calculate_basic_percentage(site['motif'])
            
            if 'aromaticPercentage' not in site:
                site['aromaticPercentage'] = calculate_aromatic_percentage(site['motif'])
            
            if 'hydrophobicityScore' not in site:
                site['hydrophobicityScore'] = calculate_hydrophobicity_score(site['motif'])

        # Extract known kinase information
        known_kinase = None
        for i in range(1, 6):
            kinase_field = f"KINASE_{i}"
            if kinase_field in site and site[kinase_field] and site[kinase_field] != 'unlabeled':
                known_kinase = site[kinase_field]
                break
                
        # Check other possible kinase field names
        if not known_kinase:
            for field in ['known_kinase', 'knownKinase']:
                if field in site and site[field] and site[field] != 'unlabeled':
                    known_kinase = site[field]
                    break
                    
        site['known_kinase'] = known_kinase or ''
    
    # Generate the HTML
    html = """
    <div class="card mt-4">
        <div class="card-header">
            <h5 class="mb-0 d-flex align-items-center">
                Phosphorylation Site Analysis
                <small class="ms-4 text-muted" style="font-size: 0.85rem;">
                    <!-- Green legend box -->
                    <span style="background-color: #c8e6c9; display: inline-block; width: 15px; height: 15px; margin-right: 5px; border: 1px solid #bbb;"></span>
                    Has Structural Similarity Data
                    &nbsp;&nbsp;
                    <!-- Orange legend box -->
                    <span style="background-color: #ffcc80; display: inline-block; width: 15px; height: 15px; margin-right: 5px; border: 1px solid #bbb;"></span>
                    No Structural Similarity Data
                </small>
            </h5>
        </div>
        <div class="card-body p-0">
            <div class="table-responsive">
                <table class="table table-striped table-hover phosphosite-table">
                    <thead class="thead-light">
                        <tr>
                            <th>Site</th>
                            <th>Motif (-7 to +7)</th>
                            <th>Mean pLDDT</th>
                            <th>Site pLDDT</th>
                            <th>Nearby Residues (10Å)</th>
                            <th>Surface Access.</th>
                            <th>Known</th>
                            <th>Kinase</th>
                            <th>Actions</th>
                        </tr>
                    </thead>
                    <tbody id="phosphosite-table">
    """
    
    for site in phosphosites:
        # Extract all the available metrics
        site_type = site.get('siteType', site.get('site', '')[0])
        resno = site.get('resno', 0)
        
        # Handle nearby count - make sure it's a valid number
        nearby_count = site.get('nearbyCount', site.get('nearby_count', 0))
        if isinstance(nearby_count, (str, float)) and (str(nearby_count).lower() == 'nan' or math.isnan(float(nearby_count))):
            nearby_count = 6  # Default value for nearby residues if NaN
        
        is_known = site.get('isKnown', site.get('is_known', False))
        motif = site.get('motif', '')
        
        # Get mean pLDDT - handle string values and NaN properly
        try:
            mean_plddt_value = site.get('meanPLDDT', site.get('mean_plddt', 0))
            
            # Check for NaN values in various formats
            if isinstance(mean_plddt_value, str) and mean_plddt_value.lower() == 'nan':
                mean_plddt = 50.0  # Default value
                mean_plddt_text = f"{mean_plddt:.1f}"
            elif isinstance(mean_plddt_value, float) and math.isnan(mean_plddt_value):
                mean_plddt = 50.0  # Default value
                mean_plddt_text = f"{mean_plddt:.1f}"
            else:
                mean_plddt = float(mean_plddt_value)
                mean_plddt_text = f"{mean_plddt:.1f}"
        except (ValueError, TypeError):
            mean_plddt = 50.0  # Default value
            mean_plddt_text = f"{mean_plddt:.1f}"
        
        # Get metrics with various possible key names - handle string values and NaN properly
        try:
            surface_accessibility_value = site.get('surfaceAccessibility', site.get('surface_accessibility', 0))
            
            # Check for NaN values in various formats
            if isinstance(surface_accessibility_value, str) and surface_accessibility_value.lower() == 'nan':
                surface_accessibility = 50.0  # Default value
            elif isinstance(surface_accessibility_value, float) and math.isnan(surface_accessibility_value):
                surface_accessibility = 50.0  # Default value
            else:
                surface_accessibility = float(surface_accessibility_value)
                
            surface_access_text = f"{surface_accessibility:.1f}%"
        except (ValueError, TypeError):
            # Default value if conversion fails
            surface_accessibility = 50.0
            surface_access_text = f"{surface_accessibility:.1f}%"
        
        # Get site pLDDT - handle NaN properly
        try:
            site_plddt_value = site.get('site_plddt', mean_plddt)
            
            # Check for NaN values in various formats
            if isinstance(site_plddt_value, str) and site_plddt_value.lower() == 'nan':
                site_plddt = mean_plddt  # Use mean pLDDT as fallback
                site_plddt_text = f"{site_plddt:.1f}"
            elif isinstance(site_plddt_value, float) and math.isnan(site_plddt_value):
                site_plddt = mean_plddt  # Use mean pLDDT as fallback
                site_plddt_text = f"{site_plddt:.1f}"
            else:
                site_plddt = float(site_plddt_value)
                site_plddt_text = f"{site_plddt:.1f}"
        except (ValueError, TypeError):
            site_plddt = mean_plddt  # Use mean pLDDT as fallback
            site_plddt_text = f"{site_plddt:.1f}"
        
        # Get the kinase info for display
        known_kinase = site.get('known_kinase', '')
        
        # Determine if the site has complete structural analysis
        site_id = f"{protein_uniprot_id}_{resno}"
        from protein_explorer.analysis.phospho_analyzer_2 import get_phosphosite_data
        supp_data = get_phosphosite_data(site_id)
        if supp_data is not None:
            # Site is analyzed: color green (e.g., light green)
            row_style = 'style="background-color: #c8e6c9;"'
        else:
            # Site not analyzed: color orange
            row_style = 'style="background-color: #ffcc80;"'

        # Add all data attributes including kinase information
        data_attrs = f"""
            data-site="{site.get('site', '')}"
            data-resno="{resno}"
            data-type="{site_type}"
            data-nearby="{nearby_count}"
            data-plddt="{mean_plddt}"
            data-surface="{surface_accessibility}"
            data-acidic="{site.get('acidicPercentage', 0)}"
            data-basic="{site.get('basicPercentage', 0)}"
            data-aromatic="{site.get('aromaticPercentage', 0)}"
            data-bfactor="{site.get('bFactorGradient', 0)}"
            data-hydrophobicity="{site.get('hydrophobicityScore', 0)}"
            data-known="{is_known}"
            data-kinase="{known_kinase}"
            data-site-id="{site_id}"
        """

        # Include the Kinase column in the row rendering
        html += f"""
        <tr {row_style} {data_attrs}>
            <td><a href="/site/{protein_uniprot_id}/{site.get('site', '')}" class="site-link" data-resno="{resno}"><strong id="site-{resno}">{site.get('site', '')}</strong></a></td>
            <td><code class="motif-sequence">{motif}</code></td>
            <td>{mean_plddt_text}</td>
            <td>{site_plddt_text}</td>
            <td>{nearby_count}</td>
            <td>{surface_access_text}</td>
            <td>{"Yes" if is_known else "No"}</td>
            <td class="kinase-col">{known_kinase if known_kinase else "—"}</td>
            <td>
                <a href="/site/{protein_uniprot_id}/{site.get('site', '')}" class="btn btn-sm btn-outline-primary">
                    Details
                </a>
            </td>
        </tr>
        """
    
    html += """
                    </tbody>
                </table>
            </div>
        </div>
    </div>
    """
    
    # Rest of the function remains the same...
    html += """
    <script>
        document.addEventListener('DOMContentLoaded', function() {
            // Add click handlers to site links
            const siteLinks = document.querySelectorAll('.site-link');
            siteLinks.forEach(link => {
                link.addEventListener('click', function(e) {
                    e.preventDefault();
                    const resno = this.getAttribute('data-resno');
                    
                    // Find the span in the sequence viewer
                    const sequenceSpans = document.querySelectorAll('.sequence-viewer span');
                    if (sequenceSpans.length > 0) {
                        // Find and click the span for this residue
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
    
    return html
    
def calculate_acidic_percentage(motif, window=5):
    """Calculate percentage of acidic residues (D, E) in a window around the phosphosite."""
    if not motif or len(motif) < 3:
        return 0
    
    # Find the center position (phosphosite)
    center_pos = len(motif) // 2
    
    # Define the window around the center (-window to +window)
    start_pos = max(0, center_pos - window)
    end_pos = min(len(motif), center_pos + window + 1)
    
    # Extract the window
    window_motif = motif[start_pos:end_pos]
    
    # Count acidic residues
    acidic_residues = 'DE'
    count = 0
    total = 0
    
    for i, aa in enumerate(window_motif):
        # Skip the center residue (phosphosite)
        if i == center_pos - start_pos:
            continue
            
        total += 1
        if aa in acidic_residues:
            count += 1
    
    return (count / total * 100) if total > 0 else 0

def calculate_basic_percentage(motif, window=5):
    """Calculate percentage of basic residues (K, R, H) in a window around the phosphosite."""
    if not motif or len(motif) < 3:
        return 0
    
    # Find the center position (phosphosite)
    center_pos = len(motif) // 2
    
    # Define the window around the center (-window to +window)
    start_pos = max(0, center_pos - window)
    end_pos = min(len(motif), center_pos + window + 1)
    
    # Extract the window
    window_motif = motif[start_pos:end_pos]
    
    # Count basic residues
    basic_residues = 'KRH'
    count = 0
    total = 0
    
    for i, aa in enumerate(window_motif):
        # Skip the center residue (phosphosite)
        if i == center_pos - start_pos:
            continue
            
        total += 1
        if aa in basic_residues:
            count += 1
    
    return (count / total * 100) if total > 0 else 0

def calculate_aromatic_percentage(motif, window=5):
    """Calculate percentage of aromatic residues (F, W, Y) in a window around the phosphosite."""
    if not motif or len(motif) < 3:
        return 0
    
    # Find the center position (phosphosite)
    center_pos = len(motif) // 2
    
    # Define the window around the center (-window to +window)
    start_pos = max(0, center_pos - window)
    end_pos = min(len(motif), center_pos + window + 1)
    
    # Extract the window
    window_motif = motif[start_pos:end_pos]
    
    # Count aromatic residues
    aromatic_residues = 'FWY'
    count = 0
    total = 0
    
    for i, aa in enumerate(window_motif):
        # Skip the center residue (phosphosite)
        if i == center_pos - start_pos:
            continue
            
        total += 1
        if aa in aromatic_residues:
            count += 1
    
    return (count / total * 100) if total > 0 else 0

def calculate_hydrophobicity_score(motif, window=5):
    """Calculate a hydrophobicity score for the region around the phosphosite."""
    if not motif or len(motif) < 3:
        return 50  # Default middle value
    
    # Simple Kyte-Doolittle hydrophobicity scale
    hydrophobicity_scale = {
        'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
        'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
        'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
    }
    
    # Find the center position (phosphosite)
    center_pos = len(motif) // 2
    
    # Define the window around the center (-window to +window)
    start_pos = max(0, center_pos - window)
    end_pos = min(len(motif), center_pos + window + 1)
    
    # Extract the window
    window_motif = motif[start_pos:end_pos]
    
    # Calculate average hydrophobicity
    sum_hydrophobicity = 0
    count = 0
    
    for i, aa in enumerate(window_motif):
        # Skip the center residue (phosphosite)
        if i == center_pos - start_pos:
            continue
            
        if aa in hydrophobicity_scale:
            sum_hydrophobicity += hydrophobicity_scale[aa]
            count += 1
    
    avg_hydrophobicity = sum_hydrophobicity / count if count > 0 else 0
    
    # Normalize to 0-100 scale (from approximately -4.5 to 4.5 range)
    return (avg_hydrophobicity + 4.5) * 100 / 9