"""
Functions for analyzing potential phosphorylation sites in proteins.
"""

import numpy as np
import requests
import json
import os
from typing import Dict, List, Optional, Tuple
import requests

# Additional fix for analyze_phosphosites in protein_explorer/analysis/phospho.py
# to prevent it from returning NaN values

def analyze_phosphosites(sequence: str, 
                        structure_data: str, 
                        uniprot_id: str = None) -> List[Dict]:
    """
    Analyze potential phosphorylation sites (S, T, Y) in a protein sequence.
    With improved surface accessibility calculation and NaN handling.
    
    Args:
        sequence: Protein sequence string
        structure_data: PDB format data as string
        uniprot_id: UniProt ID for checking PhosphositePlus (optional)
        
    Returns:
        List of dictionaries with phosphosite information
    """
    from Bio.PDB import PDBParser, NeighborSearch, Selection
    import io
    import numpy as np
    import math
    import requests
    
    # Initialize known phosphosites dictionary
    known_sites = {}
    
    # Attempt to fetch known phosphosites if UniProt ID is provided
    if uniprot_id:
        try:
            # Use updated UniProt REST API URL
            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"
            
            # Perform the request
            response = requests.get(url)
            
            # Ensure successful response
            if response.status_code == 200:
                # Split the response into lines
                lines = response.text.split('\n')
                
                # Iterate through lines
                for i in range(len(lines)):
                    # Look for MOD_RES lines
                    if lines[i].startswith('FT   MOD_RES'):
                        # Extract site number
                        parts = lines[i].split()
                        
                        if len(parts) >= 3:
                            current_site = int(parts[2])
                            
                            # Check next line for phosphorylation
                            if i + 1 < len(lines):
                                next_line = lines[i + 1]
                                
                                # Check for phosphorylation types
                                if 'Phosphothreonine' in next_line:
                                    known_sites[f'T{current_site}'] = True
                                elif 'Phosphoserine' in next_line:
                                    known_sites[f'S{current_site}'] = True
                                elif 'Phosphotyrosine' in next_line:
                                    known_sites[f'Y{current_site}'] = True
        except Exception as e:
            print(f"Error fetching known phosphosites for {uniprot_id}: {e}")
    
    # Parse PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", io.StringIO(structure_data))
    
    # Extract coordinates and B-factors (PLDDT scores in AlphaFold)
    atoms = list(structure.get_atoms())
    ns = NeighborSearch(atoms)
    
    # Get B-factors (PLDDT) by residue
    plddt_by_residue = {}
    for residue in structure.get_residues():
        b_factors = [atom.get_bfactor() for atom in residue]
        if b_factors:
            plddt_by_residue[residue.get_id()[1]] = np.mean(b_factors)
    
    # Find all S, T, Y residues
    phosphosites = []
    for i, aa in enumerate(sequence):
        if aa in ['S', 'T', 'Y']:
            resno = i + 1  # 1-based residue numbering
            
            # Extract motif (-7 to +7)
            motif_start = max(0, i - 7)
            motif_end = min(len(sequence), i + 8)
            motif = sequence[motif_start:motif_end]
            motif = motif.upper()  # Ensure motif is uppercase
            
            # Calculate mean pLDDT for the motif
            motif_resno_range = range(motif_start + 1, motif_end + 1)
            motif_plddts = [plddt_by_residue.get(j, 0) for j in motif_resno_range]
            
            # Check for NaN or empty values and use defaults if necessary
            valid_plddts = [p for p in motif_plddts if p > 0 and not (isinstance(p, float) and math.isnan(p))]
            if valid_plddts:
                mean_plddt = np.mean(valid_plddts)
            else:
                mean_plddt = 50.0  # Default value if no valid pLDDT values
            
            # Format mean pLDDT to avoid NaN in output
            if math.isnan(mean_plddt):
                mean_plddt = 50.0  # Default value if calculation resulted in NaN
            
            # Count residues within 10Å and calculate surface accessibility
            nearby_count = 0
            surface_accessibility = None
            try:
                # Get center atom for the residue
                target_residue = None
                for res in structure.get_residues():
                    if res.get_id()[1] == resno:
                        target_residue = res
                        break

                if target_residue:
                    center_atom = target_residue['CA'] if 'CA' in target_residue else next(target_residue.get_atoms())
                    nearby_atoms = ns.search(center_atom.get_coord(), 10)  # 10Å radius
                    
                    # Use a set to automatically track unique residues
                    nearby_residues = set()
                    for atom in nearby_atoms:
                        parent = atom.get_parent()
                        full_id = parent.get_id()
                        res_id = full_id[1]
                        if res_id != resno:
                            nearby_residues.add(res_id)
                            
                    nearby_count = len(nearby_residues)
                    
                    # Calculate surface accessibility based on nearby residue count
                    # Fewer nearby residues = more exposed to solvent
                    max_neighbors = 30  # Approximate maximum reasonable number of neighbors in 10Å
                    surface_accessibility = max(0, min(100, (max_neighbors - nearby_count) / max_neighbors * 100))
            except Exception as e:
                print(f"Error calculating neighbors for {aa}{resno}: {str(e)}")
                nearby_count = 6  # Default value if calculation fails
            
            # If we couldn't calculate surface accessibility, estimate it from sequence context
            if surface_accessibility is None or math.isnan(surface_accessibility):
                # Get residues in local environment (-3 to +3)
                local_start = max(0, i - 3)
                local_end = min(len(sequence), i + 4)
                local_sequence = sequence[local_start:local_end]
                
                # Calculate approximate surface accessibility based on amino acid properties
                # Hydrophobic residues tend to be buried, charged/polar tend to be exposed
                hydrophobic = "AVILMFYWC"
                charged = "DEKRH"
                polar = "STNQ"
                
                # Count residue types in local environment
                hydrophobic_count = sum(1 for aa in local_sequence if aa in hydrophobic)
                charged_count = sum(1 for aa in local_sequence if aa in charged)
                polar_count = sum(1 for aa in local_sequence if aa in polar)
                
                # More hydrophobic = less accessible, more charged/polar = more accessible
                local_length = len(local_sequence)
                hydrophobic_fraction = hydrophobic_count / local_length
                charged_polar_fraction = (charged_count + polar_count) / local_length
                
                # Estimate surface accessibility (0-100%)
                surface_accessibility = min(100, max(0, 50 + (charged_polar_fraction - hydrophobic_fraction) * 100))
            
            # Check if known in PhosphositePlus
            site_key = f"{aa}{resno}"
            is_known = site_key in known_sites
            
            # Get site-specific pLDDT if possible
            try:
                site_plddt = plddt_by_residue.get(resno, mean_plddt)
                # Check for NaN
                if math.isnan(site_plddt):
                    site_plddt = mean_plddt  # Fall back to mean pLDDT
            except (TypeError, ValueError):
                site_plddt = mean_plddt  # Fall back to mean pLDDT
            
            # Add to results - ensure no NaN values are stored
            phosphosites.append({
                'site': site_key,
                'resno': resno,
                'motif': motif,
                'mean_plddt': f"{mean_plddt:.1f}",
                'site_plddt': site_plddt,
                'nearby_count': nearby_count,
                'surface_accessibility': surface_accessibility,
                'is_known': is_known
            })
    
    return phosphosites