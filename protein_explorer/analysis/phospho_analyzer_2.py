"""
Phosphosite Analyzer module - Handles phosphosite structural analysis for protein structures.

This module provides functions to analyze phosphorylation sites in proteins
and find structural similarities with other known sites.
"""

import os
import pandas as pd
import requests
from typing import Dict, List, Optional, Union, Tuple
import logging
from protein_explorer.analysis.phospho import analyze_phosphosites
from protein_explorer.data.scaffold import get_protein_by_id, get_alphafold_structure
from protein_explorer.db.db import (
    get_phosphosite_data, get_phosphosites_batch,
    find_structural_matches, find_structural_matches_batch
)
import re

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_phosphosite_supp_data(file_path: str = None) -> pd.DataFrame:
    """
    Load the phosphosite supplementary data.
    
    This function now checks the database instead of file.
    
    Args:
        file_path: Parameter kept for backward compatibility but not used
        
    Returns:
        True if data is available in the database, None otherwise
    """
    # This function is kept for backward compatibility but now uses the database
    logger.info("Using database for phosphosite supplementary data")
    
    # We just return True to indicate we can access the data
    # The actual data access happens via get_phosphosite_data()
    return True

def enhance_phosphosite(phosphosite: Dict, uniprot_id: str) -> Dict:
    """
    Enhance a phosphosite dictionary with supplementary data.
    
    Args:
        phosphosite: Dictionary with phosphosite info
        uniprot_id: UniProt ID
        
    Returns:
        Enhanced phosphosite dictionary
    """
    if 'resno' not in phosphosite:
        return phosphosite
    
    # Get site ID
    site_id = f"{uniprot_id}_{phosphosite['resno']}"
    
    # Get supplementary data
    supp_data = get_phosphosite_data(site_id)
    if not supp_data:
        return phosphosite
    
    # Create a new dictionary to avoid modifying the original
    enhanced_site = phosphosite.copy()
    
    # Enhance with supplementary data
    if 'motif_plddt' in supp_data and supp_data['motif_plddt'] is not None:
        enhanced_site['mean_plddt'] = f"{supp_data['motif_plddt']:.1f}"
    
    if 'nearby_count' in supp_data and supp_data['nearby_count'] is not None:
        enhanced_site['nearby_count'] = supp_data['nearby_count']
    
    if 'SITE_+/-7_AA' in supp_data and supp_data['SITE_+/-7_AA'] is not None:
        enhanced_site['motif'] = supp_data['SITE_+/-7_AA']
    
    # Add any other supplementary fields
    for key in ['site_plddt', 'surface_accessibility', 'polar_aa_percent', 'nonpolar_aa_percent',
                'acidic_aa_percent', 'basic_aa_percent', 'secondary_structure']:
        source_key = key.lower().replace('_', '_')
        if source_key in supp_data and supp_data[source_key] is not None:
            enhanced_site[key] = supp_data[source_key]
    
    return enhanced_site

def enhance_structural_matches(matches, site):
    """
    Enhance structural matches with supplementary data.
    Ensures motif sequences in matches are uppercase.
    
    Args:
        matches: List of structural match dictionaries
        site: Query site string for logging
        
    Returns:
        Enhanced list of matches
    """
    if not matches:
        return []
    
    logger.info(f"Enhancing {len(matches)} structural matches for {site}")
    
    enhanced_matches = []
    for match in matches:
        # Skip self-matches (RMSD ≈ 0)
        if match.get('rmsd', 0) < 0.01:
            continue
        
        # Get target info
        target_uniprot = match.get('target_uniprot')
        target_site = match.get('target_site')
        
        # Create enhanced match (copy the original)
        enhanced_match = match.copy()
        
        # Parse site number from target_site
        site_match = re.match(r'([A-Z]?)(\d+)', target_site)
        if site_match:
            site_type = site_match.group(1) or 'S'  # Default to S if missing
            resno = int(site_match.group(2))
            target_id = f"{target_uniprot}_{resno}"
            
            # Get supplementary data
            target_supp = get_phosphosites_batch(target_id)
            if target_supp:
                # Add supplementary data
                if 'motif_plddt' in target_supp and target_supp['motif_plddt'] is not None:
                    enhanced_match['plddt'] = f"{target_supp['motif_plddt']:.1f}"
                
                if 'nearby_count' in target_supp and target_supp['nearby_count'] is not None:
                    enhanced_match['nearby_count'] = target_supp['nearby_count']
                
                # Try both motif field names
                if 'motif' in target_supp and target_supp['motif'] is not None:
                    enhanced_match['motif'] = target_supp['motif'].upper()
                elif 'SITE_+/-7_AA' in target_supp and target_supp['SITE_+/-7_AA'] is not None:
                    enhanced_match['motif'] = target_supp['SITE_+/-7_AA'].upper()
                
                
                # Add additional fields
                for key in ['site_plddt', 'surface_accessibility', 'secondary_structure']:
                    if key in target_supp and target_supp[key] is not None:
                        enhanced_match[key] = target_supp[key]
        
        # Make sure any existing motif is uppercase
        if 'motif' in enhanced_match and enhanced_match['motif']:
            enhanced_match['motif'] = enhanced_match['motif'].upper()
            
        enhanced_matches.append(enhanced_match)
    
    return enhanced_matches

def enhance_structural_matches_group(matches, site):
    """
    Enhance structural matches with supplementary data using batch processing.
    Ensures motif sequences in matches are uppercase.
    
    Args:
        matches: List of structural match dictionaries
        site: Query site string for logging
        
    Returns:
        Enhanced list of matches
    """
    if not matches:
        return []
    
    logger.info(f"Enhancing {len(matches)} structural matches for {site}")
    
    # Filter out self-matches (RMSD ≈ 0)
    filtered_matches = [match for match in matches if match.get('rmsd', 0) >= 0.01]
    
    # Prepare target_ids for batch processing
    target_ids = []
    match_indices = {}  # To map target_ids back to their matches
    
    for i, match in enumerate(filtered_matches):
        target_uniprot = match.get('target_uniprot')
        target_site = match.get('target_site')
        
        # Parse site number from target_site
        site_match = re.match(r'([A-Z]?)(\d+)', target_site)
        if site_match:
            site_type = site_match.group(1) or 'S'  # Default to S if missing
            resno = int(site_match.group(2))
            target_id = f"{target_uniprot}_{resno}"
            
            target_ids.append(target_id)
            match_indices[target_id] = i
    
    # Get supplementary data for all targets in one batch call
    all_supp_data = get_phosphosites_batch(target_ids)
    
    # Create enhanced matches
    enhanced_matches = []
    
    for i, match in enumerate(filtered_matches):
        # Create enhanced match (copy the original)
        enhanced_match = match.copy()
        
        target_uniprot = match.get('target_uniprot')
        target_site = match.get('target_site')
        
        # Parse site number from target_site
        site_match = re.match(r'([A-Z]?)(\d+)', target_site)
        if site_match:
            site_type = site_match.group(1) or 'S'  # Default to S if missing
            resno = int(site_match.group(2))
            target_id = f"{target_uniprot}_{resno}"
            
            # Get supplementary data from the batch results
            target_supp = all_supp_data.get(target_id)
            if target_supp:
                # Add supplementary data
                if 'motif_plddt' in target_supp and target_supp['motif_plddt'] is not None:
                    enhanced_match['plddt'] = f"{target_supp['motif_plddt']:.1f}"
                
                if 'nearby_count' in target_supp and target_supp['nearby_count'] is not None:
                    enhanced_match['nearby_count'] = target_supp['nearby_count']
                
                # Try both motif field names
                if 'motif' in target_supp and target_supp['motif'] is not None:
                    enhanced_match['motif'] = target_supp['motif'].upper()
                elif 'SITE_+/-7_AA' in target_supp and target_supp['SITE_+/-7_AA'] is not None:
                    enhanced_match['motif'] = target_supp['SITE_+/-7_AA'].upper()
                
                # Add additional fields
                for key in ['site_plddt', 'surface_accessibility', 'secondary_structure']:
                    if key in target_supp and target_supp[key] is not None:
                        enhanced_match[key] = target_supp[key]
        
        # Make sure any existing motif is uppercase
        if 'motif' in enhanced_match and enhanced_match['motif']:
            enhanced_match['motif'] = enhanced_match['motif'].upper()
            
        enhanced_matches.append(enhanced_match)
    
    return enhanced_matches


def enhance_structural_matches_optimized(matches, site):
    """
    Enhance structural matches with supplementary data using the
    new STY_Structural_Annotations table.
    Ensures motif sequences in matches are uppercase.
    
    Args:
        matches: List of structural match dictionaries
        site: Query site string for logging
        
    Returns:
        Enhanced list of matches
    """
    if not matches:
        return []
    
    logger.info(f"Enhancing {len(matches)} structural matches for {site}")
    
    # Get all target IDs for batch query
    target_ids = []
    for match in matches:
        target_uniprot = match.get('target_uniprot')
        target_site = match.get('target_site')
        
        if target_uniprot and target_site:
            # Extract number from the target site
            site_number = ''.join(filter(str.isdigit, target_site))
            if site_number:
                target_id = f"{target_uniprot}_{site_number}"
                target_ids.append(target_id)
    
    # Get comprehensive data for all targets in one batch query
    from protein_explorer.db.db import get_comprehensive_site_data_batch
    all_target_data = get_comprehensive_site_data_batch(target_ids)
    
    # Process each match with the enhanced data
    enhanced_matches = []
    for match in matches:
        # Skip self-matches (RMSD ≈ 0)
        if match.get('rmsd', 0) < 0.01:
            continue
            
        # Create a copy of the match to enhance
        enhanced_match = match.copy()
        
        # Get target info
        target_uniprot = match.get('target_uniprot')
        target_site = match.get('target_site')
        
        # Skip if missing target info
        if not target_uniprot or not target_site:
            continue
            
        # Extract number from the target site
        site_number = ''.join(filter(str.isdigit, target_site))
        if not site_number:
            continue
            
        # Create target ID
        target_id = f"{target_uniprot}_{site_number}"
        
        # Get comprehensive data for this target
        target_data = all_target_data.get(target_id, {})
        
        if target_data:
            # Add data from STY_Structural_Annotations
            if 'pLDDT' in target_data:
                enhanced_match['plddt'] = float(target_data['pLDDT'])
            
            if 'NeighborCount' in target_data:
                enhanced_match['nearby_count'] = int(target_data['NeighborCount'])
            
            if 'HydroxylExposure' in target_data:
                # Convert to percentage for consistency
                enhanced_match['surface_accessibility'] = float(target_data['HydroxylExposure'] * 100)
            
            if 'SecondaryStructure' in target_data:
                enhanced_match['secondary_structure'] = target_data['SecondaryStructure']
            
            # Get motif - try different field names
            if 'Motif' in target_data and target_data['Motif']:
                enhanced_match['motif'] = target_data['Motif'].upper()
            elif 'SITE_+/-7_AA' in target_data and target_data['SITE_+/-7_AA']:
                enhanced_match['motif'] = target_data['SITE_+/-7_AA'].upper()
            
            # Check for kinase information
            known_kinases = []
            for i in range(1, 6):
                kinase_field = f"KINASE_{i}"
                if kinase_field in target_data and target_data[kinase_field] and target_data[kinase_field] != 'unlabeled':
                    known_kinases.append(target_data[kinase_field])
            
            if known_kinases:
                enhanced_match['known_kinase'] = ', '.join(known_kinases)
        
        # Make sure any existing motif is uppercase
        if 'motif' in enhanced_match and enhanced_match['motif']:
            enhanced_match['motif'] = enhanced_match['motif'].upper()
            
        enhanced_matches.append(enhanced_match)
    
    return enhanced_matches

def enhance_structural_matches_group_optimized(matches, site):
    """
    Enhance structural matches with supplementary data using batch processing.
    Optimized version that uses STY_Structural_Annotations and comprehensive site data.
    
    Args:
        matches: List of structural match dictionaries
        site: Query site string for logging
        
    Returns:
        Enhanced list of matches
    """
    if not matches:
        return []
    
    logger.info(f"Enhancing {len(matches)} structural matches for {site}")
    
    # Filter out self-matches (RMSD ≈ 0)
    filtered_matches = [match for match in matches if match.get('rmsd', 0) >= 0.01]
    
    # Prepare target_ids for batch query
    target_ids = []
    match_indices = {}  # To map target_ids back to their matches
    
    for i, match in enumerate(filtered_matches):
        target_uniprot = match.get('target_uniprot')
        target_site = match.get('target_site')
        
        # Parse site number from target_site
        site_match = re.match(r'([A-Z]?)(\d+)', target_site)
        if site_match:
            site_type = site_match.group(1) or 'S'  # Default to S if missing
            resno = int(site_match.group(2))
            target_id = f"{target_uniprot}_{resno}"
            
            target_ids.append(target_id)
            match_indices[target_id] = i
    
    # Get comprehensive data for all targets in one batch query
    from protein_explorer.db.db import get_comprehensive_site_data_batch
    all_target_data = get_comprehensive_site_data_batch(target_ids)
    
    # Create enhanced matches
    enhanced_matches = []
    
    for i, match in enumerate(filtered_matches):
        # Create enhanced match (copy the original)
        enhanced_match = match.copy()
        
        target_uniprot = match.get('target_uniprot')
        target_site = match.get('target_site')
        
        # Parse site number from target_site
        site_match = re.match(r'([A-Z]?)(\d+)', target_site)
        if site_match:
            site_type = site_match.group(1) or 'S'  # Default to S if missing
            resno = int(site_match.group(2))
            target_id = f"{target_uniprot}_{resno}"
            
            # Get comprehensive data from the batch results
            target_data = all_target_data.get(target_id)
            if target_data:
                # Add data from STY_Structural_Annotations
                if 'pLDDT' in target_data:
                    enhanced_match['plddt'] = float(target_data['pLDDT'])
                
                if 'NeighborCount' in target_data:
                    enhanced_match['nearby_count'] = int(target_data['NeighborCount'])
                
                if 'HydroxylExposure' in target_data:
                    # Convert to percentage for consistency
                    enhanced_match['surface_accessibility'] = float(target_data['HydroxylExposure'] * 100)
                
                if 'SecondaryStructure' in target_data:
                    enhanced_match['secondary_structure'] = target_data['SecondaryStructure']
                
                # Get motif - try different field names
                if 'Motif' in target_data and target_data['Motif']:
                    enhanced_match['motif'] = target_data['Motif'].upper()
                elif 'SITE_+/-7_AA' in target_data and target_data['SITE_+/-7_AA']:
                    enhanced_match['motif'] = target_data['SITE_+/-7_AA'].upper()
                
                # Check for kinase information
                known_kinases = []
                for i in range(1, 6):
                    kinase_field = f"KINASE_{i}"
                    if kinase_field in target_data and target_data[kinase_field] and target_data[kinase_field] != 'unlabeled':
                        known_kinases.append(target_data[kinase_field])
                
                if known_kinases:
                    enhanced_match['known_kinase'] = ', '.join(known_kinases)
        
        # Make sure any existing motif is uppercase
        if 'motif' in enhanced_match and enhanced_match['motif']:
            enhanced_match['motif'] = enhanced_match['motif'].upper()
            
        enhanced_matches.append(enhanced_match)
    
    return enhanced_matches

def get_phosphosites(uniprot_id: str) -> List[Dict]:
    """
    Analyze potential phosphorylation sites for a protein, with supplementary data.
    
    Args:
        uniprot_id: UniProt ID of the protein
        
    Returns:
        List of dictionaries with phosphosite information
    """
    logger.info(f"Analyzing phosphosites for {uniprot_id}")
    
    try:
        # Get protein data
        protein_data = get_protein_by_id(uniprot_id=uniprot_id)
        
        # Get sequence
        sequence = protein_data.get('metadata', {}).get('sequence', {}).get('value')
        if not sequence:
            logger.warning(f"Protein sequence not found for {uniprot_id}")
            raise ValueError(f"Protein sequence not found for {uniprot_id}")
            
        # Get structure
        structure = get_alphafold_structure(uniprot_id)
        if not structure:
            logger.warning(f"Protein structure not found for {uniprot_id}. Checking alternative sources...")
            
            # Try a mock structure for testing purposes when no real structure is available
            # This could be replaced with other structure sources (PDB, etc.) in a production environment
            mock_structure = generate_mock_structure(sequence)
            if mock_structure:
                logger.info(f"Using mock structure for {uniprot_id}")
                structure = mock_structure
            else:
                raise ValueError(f"Protein structure not found for {uniprot_id}")
        
        # Analyze phosphosites
        phosphosites = analyze_phosphosites(sequence, structure, uniprot_id)
        
        # Get site IDs for batch lookup
        site_ids = [f"{uniprot_id}_{site['resno']}" for site in phosphosites]
        
        # Batch retrieve supplementary data
        supp_data_dict = get_phosphosites_batch(site_ids)
        
        # Enhance with supplementary data
        enhanced_sites = []
        for site in phosphosites:
            site_id = f"{uniprot_id}_{site['resno']}"
            supp_data = supp_data_dict.get(site_id)
            
            if supp_data:
                # Create enhanced site
                enhanced_site = site.copy()
                
                # Add supplementary data
                if 'motif_plddt' in supp_data and supp_data['motif_plddt'] is not None:
                    enhanced_site['mean_plddt'] = f"{supp_data['motif_plddt']:.1f}"
                
                if 'nearby_count' in supp_data and supp_data['nearby_count'] is not None:
                    enhanced_site['nearby_count'] = supp_data['nearby_count']
                
                if 'SITE_+/-7_AA' in supp_data and supp_data['SITE_+/-7_AA'] is not None:
                    enhanced_site['motif'] = supp_data['SITE_+/-7_AA']
                
                # Add additional fields
                for key in ['site_plddt', 'surface_accessibility', 'polar_aa_percent', 'nonpolar_aa_percent',
                           'acidic_aa_percent', 'basic_aa_percent', 'secondary_structure']:
                    source_key = key.lower().replace('_', '_')
                    if source_key in supp_data and supp_data[source_key] is not None:
                        enhanced_site[key] = supp_data[source_key]
                
                enhanced_sites.append(enhanced_site)
            else:
                enhanced_sites.append(site)
        
        return enhanced_sites
    except Exception as e:
        logger.error(f"Error analyzing phosphosites: {e}")
        raise ValueError(f"Error analyzing phosphosites: {e}")

def generate_mock_structure(sequence: str) -> Optional[str]:
    """
    Generate a mock PDB structure for cases where the AlphaFold structure is not available.
    This is for demonstration purposes only and should be replaced with real structures in production.
    
    Args:
        sequence: Protein sequence
        
    Returns:
        PDB format structure as string, or None if generation fails
    """
    try:
        # Create a very basic linear structure
        # This is extremely simplified and not biologically accurate
        pdb_lines = []
        
        # PDB header
        pdb_lines.append("HEADER    MOCK STRUCTURE")
        pdb_lines.append("TITLE     MOCK STRUCTURE FOR SEQUENCE")
        
        # Add atoms - just alpha carbons in a straight line
        atom_num = 1
        for i, aa in enumerate(sequence):
            x = i * 3.8  # ~3.8Å is typical CA-CA distance
            y = 0
            z = 0
            
            # B-factor (PLDDT) set to 70 (medium confidence) for all residues
            b_factor = 70.0
            
            line = f"ATOM  {atom_num:5d}  CA  {aa}   A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{b_factor:6.2f}           C  "
            pdb_lines.append(line)
            atom_num += 1
        
        # End of file
        pdb_lines.append("END")
        
        return "\n".join(pdb_lines)
    except Exception as e:
        logger.error(f"Error generating mock structure: {e}")
        return None

def find_phosphosite_structural_matches(uniprot_id: str, phosphosites: List[Dict], 
                                      top_n: int = None) -> Dict[str, List[Dict]]:
    """
    Find structural matches for phosphosites in the database.
    
    Args:
        uniprot_id: UniProt ID of the protein
        phosphosites: List of phosphosite dictionaries from analyze_phosphosites
        top_n: Number of top matches to return per site
        
    Returns:
        Dictionary mapping site IDs to lists of match dictionaries
    """
    logger.info(f"Finding structural matches for {uniprot_id}")
    
    try:
        # Create site IDs in the format UniprotID_ResNo
        site_ids = [f"{uniprot_id}_{site['resno']}" for site in phosphosites]
        
        # Get matches for all sites in a batch query
        all_matches = find_structural_matches_batch(site_ids, rmsd_threshold=5.0)
        
        # Organize matches by site
        structural_matches = {}
        
        for site in phosphosites:
            site_id = f"{uniprot_id}_{site['resno']}"
            site_matches = all_matches.get(site_id, [])
            
            # Sort by RMSD and take top N matches if specified
            site_matches.sort(key=lambda x: x.get('rmsd', float('inf')))
            if top_n is not None and len(site_matches) > top_n:
                site_matches = site_matches[:top_n]
            
            if site_matches:
                # Find the corresponding site data
                site_str = site.get('site', f"{site.get('siteType', '')}{site['resno']}")
                structural_matches[site_str] = site_matches
        
        return structural_matches
    except Exception as e:
        logger.error(f"Error finding structural matches: {e}")
        raise ValueError(f"Error finding structural matches: {e}")

def analyze_protein(identifier: str, id_type: str = 'uniprot') -> Dict:
    """
    Complete analysis of phosphosites and structural matches for a protein.
    
    Args:
        identifier: UniProt ID or gene symbol
        id_type: 'uniprot' or 'gene'
        
    Returns:
        Dictionary with protein info, phosphosites, and structural matches
    """
    try:
        # Get protein data
        if id_type.lower() == 'uniprot':
            protein_info = {'uniprot_id': identifier}
            protein_data = get_protein_by_id(uniprot_id=identifier)
        else:
            protein_data = get_protein_by_id(gene_symbol=identifier)
            
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
        
        # Try to get phosphosites
        phosphosites = []
        structural_matches = None
        error_message = None
        
        try:
            # Get phosphosites with supplementary data
            phosphosites = get_phosphosites(uniprot_id)
            
            # Find structural matches
            try:
                structural_matches = find_phosphosite_structural_matches(uniprot_id, phosphosites)
                
                # Enhance with supplementary data
                for site, matches in structural_matches.items():
                    structural_matches[site] = enhance_structural_matches(matches, site)
                    
            except Exception as e:
                error_message = f"Error analyzing structural matches: {str(e)}"
                logger.error(error_message)
        except Exception as e:
            error_message = f"Error analyzing phosphosites: {str(e)}"
            logger.error(f"Error analyzing phosphosites: {e}")
        
        # Compile results
        results = {
            'protein_info': protein_info,
            'phosphosites': phosphosites,
            'structural_matches': structural_matches,
            'error': error_message
        }
        
        return results
    except Exception as e:
        logger.error(f"Error in complete analysis: {e}")
        raise ValueError(f"Error in complete analysis: {e}")

def analyze_phosphosite_context(structure_data, site_number, site_type):
    """
    Analyze structural context around a phosphorylation site.
    
    Args:
        structure_data: PDB format data as string
        site_number: The residue number
        site_type: The residue type (S, T, or Y)
        
    Returns:
        Dictionary with structural context information
    """
    from Bio.PDB import PDBParser, NeighborSearch, Selection, Vector
    import io
    import numpy as np
    print("a")
    # Parse the PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", io.StringIO(structure_data))
    print("b")
    # Get all atoms
    all_atoms = list(structure.get_atoms())
    ns = NeighborSearch(all_atoms)
    # Find the target residue
    target_residue = None
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[1] == site_number:
                    target_residue = residue
                    break
    if not target_residue:
        return {"error": f"Site {site_type}{site_number} not found in structure"}
    # Get center atom for the residue (CA or first atom)
    if 'CA' in target_residue:
        center_atom = target_residue['CA']
    else:
        # Use the first atom if CA not available
        center_atom = next(target_residue.get_atoms())
    # Get the residue coordinates
    center_coords = center_atom.get_coord()
    
    # Find nearby residues (8Å radius)
    nearby_atoms = ns.search(center_coords, 10.0)
    # Group nearby atoms by residue
    nearby_residues = {}
    for atom in nearby_atoms:
        residue = atom.get_parent()
        resno = residue.get_id()[1]
        
        # Skip the target residue itself
        if resno == site_number:
            continue
            
        residue_name = residue.get_resname()
        
        # Add to nearby residues dictionary
        if resno not in nearby_residues:
            nearby_residues[resno] = {
                "resname": residue_name,
                "atoms": [],
                "min_distance": float('inf')
            }
            
        # Store the atom and its distance
        dist = np.linalg.norm(atom.get_coord() - center_coords)
        nearby_residues[resno]["atoms"].append({
            "atom_name": atom.get_name(),
            "distance": dist
        })
        
        # Update minimum distance
        if dist < nearby_residues[resno]["min_distance"]:
            nearby_residues[resno]["min_distance"] = dist
    # Sort nearby residues by distance
    sorted_nearby = sorted(
        nearby_residues.items(), 
        key=lambda x: x[1]["min_distance"]
    )
    # Prepare results
    nearby_info = [
        {
            "resno": resno,
            "resname": data["resname"],
            "min_distance": round(data["min_distance"], 2),
            "atoms": len(data["atoms"])
        }
        for resno, data in sorted_nearby
    ]
    
    # Count amino acid types within contact distance (5Å)
    amino_acid_groups = {
        "polar": ["SER", "THR", "TYR", "CYS", "ASN", "GLN"],
        "nonpolar": ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO", "GLY"],
        "acidic": ["ASP", "GLU"],
        "basic": ["LYS", "ARG", "HIS"]
    }
    
    contact_counts = {group: 0 for group in amino_acid_groups}
    
    for resno, data in sorted_nearby:
        if data["min_distance"] <= 5.0:  # Only count residues within 5Å
            for group, residues in amino_acid_groups.items():
                if data["resname"] in residues:
                    contact_counts[group] += 1
                    break
    # Calculate secondary structure
    #try:
    #    from Bio.PDB.DSSP import DSSP
    #    model = structure[0]  # Use first model
    #    dssp = DSSP(model, io.StringIO(structure_data))
        
        # Get DSSP data for the target residue
    #    chain_id = list(model.child_dict.keys())[0]  # Get first chain ID
    #    dssp_key = (chain_id, site_number)
        
    #    if dssp_key in dssp:
            # DSSP assigns: H (alpha helix), B (beta bridge), E (strand),
            # G (3-10 helix), I (pi helix), T (turn), S (bend), or - (other)
    #        ss_code = dssp[dssp_key][1]
            
            # Simplify to three main categories
    #        if ss_code in ['H', 'G', 'I']:
    #            secondary_structure = 'Helix'
    #        elif ss_code in ['E', 'B']:
    #            secondary_structure = 'Sheet'
    #        else:
    #            secondary_structure = 'Loop'
    #    else:
    #        secondary_structure = 'Unknown'
    #except:
        # If DSSP fails, leave as unknown
    #    secondary_structure = 'Unknown'
    
    # Calculate solvent accessibility
    # We'll use a simple proxy based on neighbor count
    # Fewer neighbors = more exposed
    max_neighbors = 30  # Approximate maximum reasonable number of neighbors in 8Å
    nearby_count = len(nearby_residues)
    solvent_accessibility = max(0, min(100, (max_neighbors - nearby_count) / max_neighbors * 100))
    
    # Extract B-factor (pLDDT in AlphaFold) for the target residue
    b_factors = [atom.get_bfactor() for atom in target_residue]
    mean_plddt = sum(b_factors) / len(b_factors) if b_factors else None
    return {
        "site": f"{site_type}{site_number}",
        "nearby_residues": nearby_info[:10],  # Show top 10
        "nearby_count": len(nearby_info),
        "contact_distribution": contact_counts,
        #"secondary_structure": secondary_structure,
        "solvent_accessibility": round(solvent_accessibility, 1),
        "plddt": round(mean_plddt, 1) if mean_plddt else None
    }

def enhance_site_visualization(uniprot_id, site, supplementary_data=None):
    """
    Create an enhanced visualization of a phosphorylation site.
    Integrates supplementary structural data and highlights structural features.
    
    Args:
        uniprot_id: UniProt ID of the protein
        site: Site identifier (e.g., "S15")
        supplementary_data: Supplementary data for the site
        
    Returns:
        HTML/JavaScript code for the visualization
    """
    from protein_explorer.data.scaffold import get_alphafold_structure
    import re
    import base64
    
    # Parse site to get residue number and type
    site_match = re.match(r'([A-Z])(\d+)', site)
    if not site_match:
        return f"<div class='alert alert-danger'>Invalid site format: {site}</div>"
    
    site_type = site_match.group(1)
    site_number = int(site_match.group(2))
    
    # Get structure
    structure_data = get_alphafold_structure(uniprot_id)
    if not structure_data:
        return f"<div class='alert alert-danger'>Could not retrieve structure for {uniprot_id}</div>"
    
    # Base64 encode the structure
    pdb_base64 = base64.b64encode(structure_data.encode()).decode()
    
    # Analyze structural context if not provided
    if not supplementary_data:
        site_id = f"{uniprot_id}_{site_number}"
        supplementary_data = get_phosphosite_data(site_id)
    
    # Default values if data not available
    site_plddt = supplementary_data.get('site_plddt', 'N/A') if supplementary_data else 'N/A'
    surface_access = supplementary_data.get('surface_accessibility', 'N/A') if supplementary_data else 'N/A'
    nearby_count = supplementary_data.get('nearby_count', 'N/A') if supplementary_data else 'N/A'
    secondary_structure = supplementary_data.get('secondary_structure', 'Unknown') if supplementary_data else 'Unknown'
    
    # Create visualization code
    js_code = f"""
    <div class="card mb-4">
        <div class="card-header">
            <h5 class="mb-0">3D Site Visualization</h5>
        </div>
        <div class="card-body">
            <div class="row">
                <div class="col-md-8">
                    <div id="site-viewer" style="width: 100%; height: 450px; border: 1px solid #ddd; border-radius: 5px;"></div>
                </div>
                <div class="col-md-4">
                    <div class="card">
                        <div class="card-header">
                            <h6 class="mb-0">Site Information</h6>
                        </div>
                        <div class="card-body">
                            <p><strong>Site:</strong> {site}</p>
                            <p><strong>pLDDT:</strong> {site_plddt}</p>
                            <p><strong>Surface Accessibility:</strong> {surface_access}%</p>
                            <p><strong>Nearby Residues:</strong> {nearby_count}</p>
                            <p><strong>Secondary Structure:</strong> {secondary_structure}</p>
                        </div>
                    </div>
                    <div class="mt-3">
                        <button id="reset-view" class="btn btn-sm btn-outline-primary">Reset View</button>
                        <button id="toggle-view" class="btn btn-sm btn-outline-secondary">Full Protein</button>
                        <button id="toggle-color" class="btn btn-sm btn-outline-success">Color by Type</button>
                    </div>
                    <div class="mt-3">
                        <div class="d-flex align-items-center mb-1">
                            <div style="width:12px; height:12px; background-color:#FF4500; border-radius:50%; margin-right:5px;"></div>
                            <small>Target Site</small>
                        </div>
                        <div class="d-flex align-items-center mb-1">
                            <div style="width:12px; height:12px; background-color:#87CEFA; border-radius:50%; margin-right:5px;"></div>
                            <small>Polar Residues</small>
                        </div>
                        <div class="d-flex align-items-center mb-1">
                            <div style="width:12px; height:12px; background-color:#FFD700; border-radius:50%; margin-right:5px;"></div>
                            <small>Non-polar Residues</small>
                        </div>
                        <div class="d-flex align-items-center mb-1">
                            <div style="width:12px; height:12px; background-color:#FF6347; border-radius:50%; margin-right:5px;"></div>
                            <small>Acidic Residues</small>
                        </div>
                        <div class="d-flex align-items-center mb-1">
                            <div style="width:12px; height:12px; background-color:#98FB98; border-radius:50%; margin-right:5px;"></div>
                            <small>Basic Residues</small>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
    
    <script>
    document.addEventListener('DOMContentLoaded', function() {{
        // Initialize NGL viewer
        const viewer = new NGL.Stage('site-viewer', {{backgroundColor: "white"}});
        
        // Handle window resizing
        window.addEventListener('resize', function() {{
            viewer.handleResize();
        }});
        
        // Define residue type groupings
        const aminoAcidGroups = {{
            polar: ["SER", "THR", "TYR", "CYS", "ASN", "GLN"],
            nonpolar: ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO", "GLY"],
            acidic: ["ASP", "GLU"],
            basic: ["LYS", "ARG", "HIS"]
        }};
        
        // Define colors for each group
        const groupColors = {{
            polar: [135/255, 206/255, 250/255],     // Light blue
            nonpolar: [255/255, 215/255, 0/255],    // Gold
            acidic: [255/255, 99/255, 71/255],      // Tomato
            basic: [152/255, 251/255, 152/255]      // Pale green
        }};
        
        // Function to determine AA group
        function getAminoAcidGroup(resname) {{
            for (const [group, residues] of Object.entries(aminoAcidGroups)) {{
                if (residues.includes(resname)) {{
                    return group;
                }}
            }}
            return "other";
        }}
        
        // Color function based on amino acid type
        function colorByType(atom) {{
            // Special color for the target site
            if (atom.resno === {site_number}) {{
                return [1.0, 0.27, 0.0];  // #FF4500 orange-red
            }}
            
            // Color by amino acid type
            const group = getAminoAcidGroup(atom.resname);
            if (group in groupColors) {{
                return groupColors[group];
            }}
            
            // Default grey for others
            return [0.5, 0.5, 0.5];
        }}
        
        // Load structure
        const pdbBlob = new Blob([atob('{pdb_base64}')], {{type: 'text/plain'}});
        
        viewer.loadFile(pdbBlob, {{ext: 'pdb'}}).then(function(component) {{
            // Get target selection
            const siteSelection = "{site_number} and .{site_type}";
            const environmentSelection = siteSelection + " or (" + siteSelection + " around 5)";
            
            // State variables
            let isFullView = false;
            let colorMode = "element";  // "element" or "type"
            
            // Button handlers
            document.getElementById('reset-view').addEventListener('click', function() {{
                updateRepresentations();
            }});
            
            document.getElementById('toggle-view').addEventListener('click', function() {{
                isFullView = !isFullView;
                this.textContent = isFullView ? 'Site Focus' : 'Full Protein';
                updateRepresentations();
            }});
            
            document.getElementById('toggle-color').addEventListener('click', function() {{
                colorMode = colorMode === "element" ? "type" : "element";
                this.textContent = colorMode === "element" ? 'Color by Type' : 'Color by Element';
                updateRepresentations();
            }});
            
            // Update all representations based on current state
            function updateRepresentations() {{
                // Remove all existing representations
                component.removeAllRepresentations();
                
                // Add cartoon representation for entire protein
                component.addRepresentation("cartoon", {{
                    color: colorMode === "type" ? colorByType : "chainid",
                    opacity: 0.7,
                    smoothSheet: true
                }});
                
                // Add ball and stick for target residue
                component.addRepresentation("ball+stick", {{
                    sele: siteSelection,
                    color: colorMode === "type" ? colorByType : "element",
                    aspectRatio: 1.5,
                    scale: 1.2
                }});
                
                // Add licorice for environment (if not full view)
                if (!isFullView) {{
                    component.addRepresentation("licorice", {{
                        sele: environmentSelection + " and not " + siteSelection,
                        color: colorMode === "type" ? colorByType : "element",
                        opacity: 0.8,
                        scale: 0.8
                    }});
                    
                    // Add labels
                    component.addRepresentation("label", {{
                        sele: environmentSelection,
                        color: "#333333",
                        labelType: "format",
                        labelFormat: "{{resname}}{{resno}}",
                        labelGrouping: "residue",
                        attachment: "middle-center",
                        showBackground: true,
                        backgroundColor: "white",
                        backgroundOpacity: 0.5
                    }});
                }}
                
                // Set view
                if (isFullView) {{
                    component.autoView();
                }} else {{
                    component.autoView(environmentSelection, 2000);
                }}
            }}
            
            // Initial setup
            updateRepresentations();
        }}).catch(function(error) {{
            console.error("Error loading structure:", error);
            document.getElementById('site-viewer').innerHTML = 
                '<div class="alert alert-danger mt-3">Error loading structure: ' + error.message + '</div>';
        }});
    }});
    </script>
    """
    
    return js_code

def create_comparative_motif_visualization(primary_site, matches):
    """
    Create a comparative visualization of sequence motifs for the primary site
    and its structural matches, showing -5 to +5 range around the phosphosite.
    
    Args:
        primary_site: Dictionary with primary site information
        matches: List of dictionaries with match information
        
    Returns:
        HTML code for the visualization
    """
    if not primary_site or 'motif' not in primary_site:
        return "<div class='alert alert-warning'>Motif data not available for primary site</div>"

    # Get primary site motif
    primary_motif = primary_site.get('motif', '')
    primary_site_name = primary_site.get('site', 'Unknown')

    # IMPROVED: More aggressive approach to get the UniProt ID
    # First try the direct uniprot_id field
    primary_uniprot = primary_site.get('uniprot_id', '')

    # If that's empty, try other common field names
    if not primary_uniprot:
        for field in ['query_uniprot', 'protein_id', 'uniprotid', 'uniprot']:
            if field in primary_site and primary_site[field]:
                primary_uniprot = primary_site[field]
                break

    # If still empty, try to get from protein dictionary if it exists
    if not primary_uniprot and isinstance(primary_site.get('protein'), dict):
        primary_uniprot = primary_site['protein'].get('uniprot_id', '')

    # If still empty, check if it's in any matches (as query_uniprot)
    if not primary_uniprot and matches:
        for match in matches:
            if 'query_uniprot' in match and match['query_uniprot']:
                primary_uniprot = match['query_uniprot']
                break

    # If still empty, try to infer from site_id if present
    if not primary_uniprot and 'site_id' in primary_site:
        site_id = primary_site['site_id']
        if '_' in site_id:
            primary_uniprot = site_id.split('_')[0]

    # EXPLICIT DEBUG: Print what we found
    print(f"DEBUG - Primary UniProt: {primary_uniprot}")
    if not primary_uniprot:
        print("DEBUG - Failed to find UniProt ID in:", primary_site.keys())


    # Get position information for proper X padding
    primary_site_pos = None
    if primary_site_name:
        import re
        site_match = re.match(r'([STY])(\d+)', primary_site_name)
        if site_match:
            primary_site_pos = int(site_match.group(2))
    
    # Filter matches that have motif data
    valid_matches = [m for m in matches if 'motif' in m and m['motif']]
    
    if not valid_matches:
        return "<div class='alert alert-warning'>No motif data available for matches</div>"
    
    # Sort by RMSD (closest matches first)
    sorted_matches = sorted(valid_matches, key=lambda x: x.get('rmsd', float('inf')))
    
    # Take top N matches
    top_matches = sorted_matches[:10]  # Limit to 10 for visualization
    
    # Create HTML
    html = """
    <style>
        .motif-comparison {
            font-family: monospace;
            margin-bottom: 20px;
        }
        .motif-row {
            display: flex;
            align-items: center;
            margin-bottom: 5px;
        }
        .motif-label {
            width: 130px;
            font-weight: bold;
            text-align: right;
            padding-right: 10px;
            font-size: 0.9rem;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }
        .motif-sequence {
            display: flex;
        }
        .motif-aa {
            width: 24px;
            height: 24px;
            display: flex;
            align-items: center;
            justify-content: center;
            margin: 0 1px;
            border-radius: 3px;
        }
        .motif-aa.highlighted {
            background-color: #ff5722;
            color: white;
            font-weight: bold;
        }
        .motif-aa.sty {
            background-color: #bbdefb;
        }
        .motif-aa.nq {
            background-color: #b39ddb;
        }
        .motif-aa.cys {
            background-color: #ffcc80;
        }
        .motif-aa.proline {
            background-color: #81c784;
        }
        .motif-aa.nonpolar {
            background-color: #ffecb3;
        }
        .motif-aa.acidic {
            background-color: #ffcdd2;
        }
        .motif-aa.basic {
            background-color: #c8e6c9;
        }
        .motif-aa.special {
            background-color: #e1bee7;
        }
        .motif-aa.aa-x {
            background-color: #e0e0e0;
            color: #9e9e9e;
        }
        .match-info {
            margin-left: 10px;
            font-size: 12px;
            color: #333;
        }
        .motif-position {
            display: flex;
            padding-left: 130px;
            margin-bottom: 10px;
        }
        .motif-position span {
            width: 24px;
            text-align: center;
            font-size: 10px;
            color: #666;
        }
    </style>
    
    <div class="motif-comparison">
        <h5 class="mb-3">Motif Comparison</h5>
        
        <!-- Position markers - CHANGED to -5 to +5 -->
        <div class="motif-position">
    """
    
    # CHANGED: Add position markers from -5 to +5 instead of -7 to +7
    for i in range(-5, 6):
        html += f'<span>{i}</span>'
    
    html += """
        </div>
    """
    
    # Helper function to get amino acid class for coloring
    def get_aa_class(aa):
        if aa == 'X':
            return "aa-x"
        elif aa in 'STY':
            return "sty"
        elif aa in 'NQ':
            return "nq"
        elif aa == 'C':
            return "cys"
        elif aa == 'P':
            return "proline"
        elif aa in 'AVILMFWG':
            return "nonpolar"
        elif aa in 'DE':
            return "acidic"
        elif aa in 'KRH':
            return "basic"
        else:
            return "special"
    
    # Full standardization function for -7:+7 motifs
    def standardize_motif(motif, site_position=None):
        """
        Standardize a phosphosite motif to have exactly 7 positions before and after
        the phosphosite, with proper X padding.
        
        Args:
            motif (str): The motif sequence
            site_position (int, optional): The position of the site in the protein sequence
            
        Returns:
            str: The standardized motif
        """
        # Find the center position (phosphosite)
        center_pos = len(motif) // 2
        
        # Get the phosphosite and parts before/after
        site_char = motif[center_pos]
        before_site = motif[:center_pos]
        after_site = motif[center_pos + 1:]
        
        # If we have the absolute site position
        if site_position is not None:
            # Calculate padding needed at beginning based on site position
            aas_before = site_position - 1  # e.g., for S6, this is 5
            padding_needed = max(0, 7 - aas_before)
            
            # Create the before part: add X padding at beginning if needed
            if len(before_site) <= 7:
                # If we have 7 or fewer residues, use all of them with padding
                padded_before = "X" * padding_needed + before_site
            else:
                # If we have more than 7, take the last 7
                padded_before = before_site[-7:]
            
            # Create the after part: take exactly 7 chars, NO padding at end
            padded_after = after_site[:7]
        else:
            # Default behavior when we don't know the site position
            # Ensure we have exactly 7 characters before
            if len(before_site) < 7:
                padded_before = "X" * (7 - len(before_site)) + before_site
            else:
                padded_before = before_site[-7:]
            
            # Ensure we have exactly 7 characters after, no padding unless needed
            padded_after = after_site[:7]
        
        return padded_before + site_char + padded_after
    
    # NEW: Function to trim to just -5:+5 range
    def trim_to_central_range(motif_str):
        """Trim a standardized 15-char motif to just the central 11 positions (-5:+5)"""
        # Assuming the motif is standardized to 15 chars with the phosphosite at position 7 (0-indexed)
        # We want to keep positions 2-12 (0-indexed), which are -5 to +5 around the phosphosite
        return motif_str[2:13]
    
    # Modified helper function to create HTML for a motif
    def create_motif_html(motif, site_pos=None):
        # First standardize motif to full 15 chars (7+1+7)
        std_motif = standardize_motif(motif, site_pos)
        
        # Then trim to just -5:+5 range
        trimmed_motif = trim_to_central_range(std_motif)
        
        # Create HTML for each amino acid in the trimmed motif
        html = '<div class="motif-sequence" style="display: flex; flex-wrap: nowrap;">'
        for i, aa in enumerate(trimmed_motif):
            aa_class = get_aa_class(aa)
            # Center position (phosphosite) is now at position 5 (0-indexed) in the trimmed motif
            highlight_class = "highlighted" if i == 5 else aa_class
            html += f'<div class="motif-aa {highlight_class}" style="width: 24px; height: 24px; display: flex; align-items: center; justify-content: center; margin: 0 1px; border-radius: 3px;">{aa}</div>'
        html += '</div>'
        return html
    
    # Add primary site motif with UniProt ID
    html += f"""
        <div class="motif-row">
            <div class="motif-label">{primary_uniprot}_{primary_site_name}:</div>
            {create_motif_html(primary_motif, primary_site_pos)}
        </div>
    """

    # Add match motifs
    for match in top_matches:
        motif = match.get('motif', '')
        target_site = match.get('target_site', 'Unknown')
        target_uniprot = match.get('target_uniprot', 'Unknown')
        rmsd = match.get('rmsd', 0.0)
        
        # Extract site position from target_site if possible
        target_site_pos = None
        import re
        site_match = re.match(r'([STY])(\d+)', target_site)
        if site_match:
            target_site_pos = int(site_match.group(2))
        else:
            # Try another pattern like just digits
            site_match = re.match(r'(\d+)', target_site)
            if site_match:
                target_site_pos = int(site_match.group(1))
        
        html += f"""
        <div class="motif-row">
            <div class="motif-label">{target_uniprot}_{target_site}:</div>
            {create_motif_html(motif, target_site_pos)}
            <div class="match-info">
                RMSD: {rmsd:.2f}Å | <a href="/site/{target_uniprot}/{target_site}" class="text-decoration-none">View site</a>
            </div>
        </div>
        """
    
    html += """
    </div>
    """
    
    return html


def get_aa_bg_color(aa):
    """Get the background color for an amino acid based on its type."""
    if aa == 'X':
        return "#e0e0e0"  # Light gray for placeholder X
    elif aa in 'STY':
        return "#bbdefb"  # Light blue for STY
    elif aa in 'NQ':
        return "#b39ddb"  # Light purple for NQ
    elif aa == 'C':
        return "#ffcc80"  # Light orange for Cysteine
    elif aa == 'P':
        return "#81c784"  # Light green for Proline
    elif aa in 'AVILMFWG':
        return "#ffecb3"  # Light yellow for other nonpolar
    elif aa in 'DE':
        return "#ffcdd2"  # Light red for acidic
    elif aa in 'KRH':
        return "#c8e6c9"  # Pale green for basic
    else:
        return "#e1bee7"  # Light pink for special cases


def analyze_residue_distributions(structural_matches):
    """
    Analyze the distribution of residues across structural matches to identify
    potential conservation patterns.
    
    Args:
        structural_matches: List of match dictionaries
        
    Returns:
        Dictionary with analysis results
    """
    if not structural_matches:
        return None
        
    # Get motifs from matches
    motifs = []
    for match in structural_matches:
        if 'motif' in match and match['motif']:
            # We assume the phosphosite is in the middle of the motif
            motifs.append(match['motif'])
    
    if not motifs:
        return None
        
    # Determine the motif length (use longest motif)
    motif_length = max(len(m) for m in motifs)
    
    # Calculate the center position (where the phosphosite is)
    center_pos = motif_length // 2
    
    # Count amino acids at each position
    position_counts = []
    for i in range(motif_length):
        counts = {}
        for motif in motifs:
            if i < len(motif):
                aa = motif[i]
                counts[aa] = counts.get(aa, 0) + 1
        position_counts.append(counts)
    
    # Calculate frequencies and identify consensus
    frequencies = []
    consensus = []
    
    for i, counts in enumerate(position_counts):
        total = sum(counts.values())
        freq = {aa: count/total for aa, count in counts.items()}
        frequencies.append(freq)
        
        # Find most common AA
        if counts:
            max_aa = max(counts.items(), key=lambda x: x[1])
            consensus.append(max_aa[0])
        else:
            consensus.append('-')
    
    # Generate relative position labels
    positions = [i - center_pos for i in range(motif_length)]
    
    # Identify conserved positions (>50% same AA)
    conserved = []
    for i, counts in enumerate(position_counts):
        total = sum(counts.values())
        max_count = max(counts.values()) if counts else 0
        
        if max_count / total >= 0.5:
            conserved.append({
                'position': positions[i],
                'amino_acid': consensus[i],
                'frequency': max_count / total * 100
            })
    
    # Group amino acids by type
    aa_groups = {
        'polar': 'STYCNQ',
        'nonpolar': 'AVILMFWPG',
        'acidic': 'DE',
        'basic': 'KRH'
    }
    
    # Count by group at each position
    group_counts = []
    for i in range(motif_length):
        counts = {group: 0 for group in aa_groups}
        for motif in motifs:
            if i < len(motif):
                aa = motif[i]
                for group, aas in aa_groups.items():
                    if aa in aas:
                        counts[group] += 1
                        break
        group_counts.append(counts)
    
    # Generate consensus motif
    consensus_str = ''.join(consensus)
    
    return {
        'motif_count': len(motifs),
        'positions': positions,
        'consensus': consensus_str,
        'frequencies': frequencies,
        'position_counts': position_counts,
        'group_counts': group_counts,
        'conserved': conserved
    }

def find_structural_matches_with_connections(site_ids: List[str], rmsd_threshold: float = 5.0) -> Dict:
    """
    Find structural matches for multiple sites in a batch,
    including connections between the similar sites.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        rmsd_threshold: Maximum RMSD value for matches
        
    Returns:
        Dictionary with nodes, edges and required data for network visualization
    """
    if not site_ids:
        return {'nodes': [], 'links': [], 'site_matches': {}}
    
    # Get batch match results for the sites
    site_matches_dict = find_structural_matches_batch(site_ids, rmsd_threshold)
    
    # Set to track all unique match IDs for the next query
    all_match_ids = set()
    
    # Collect all matches for each site
    for site_id, matches in site_matches_dict.items():
        for match in matches:
            if 'target_id' in match:
                all_match_ids.add(match['target_id'])
    
    # Find interconnections between the matches themselves
    # This requires an additional batch query
    match_interconnections = {}
    if all_match_ids:
        # Convert to list for the query
        match_id_list = list(all_match_ids)
        
        # Find matches among the matches themselves
        match_matches = find_structural_matches_batch(match_id_list, rmsd_threshold)
        
        # Store valid interconnections
        for match_id, connections in match_matches.items():
            valid_connections = []
            for connection in connections:
                # Only include connections to other matches we already know about
                # Filter by our threshold and avoid self-connections
                if ('target_id' in connection and
                    connection['target_id'] in all_match_ids and 
                    connection['target_id'] != match_id and
                    connection['rmsd'] <= rmsd_threshold):
                    valid_connections.append(connection)
            
            if valid_connections:
                match_interconnections[match_id] = valid_connections
    
    # Build network with nodes and links
    nodes = []
    links = []
    node_map = {}  # Track nodes we've added
    link_map = {}  # Track links we've added to avoid duplicates
    
    # Get all supplementary data in one batch query
    all_site_ids = site_ids.copy()
    all_site_ids.extend(list(all_match_ids))
    
    # Get supplementary data for all sites in one batch
    supp_data_batch = get_phosphosites_batch(all_site_ids)
    
    # First add nodes for the query sites
    for site_id in site_ids:
        # Split site_id to get uniprot_id and site number
        parts = site_id.split('_')
        if len(parts) < 2:
            continue
            
        uniprot_id = parts[0]
        site_num = parts[1]
        
        # Extract site type if possible (S, T, Y)
        site_match = re.match(r'([STY])(\d+)', site_num)
        if site_match:
            site_type = site_match.group(1)
            site_num = site_match.group(2)
            site_name = f"{site_type}{site_num}"
        else:
            site_name = site_num
            site_type = None
        
        # Get supplementary data if available
        site_data = supp_data_batch.get(site_id, {})
        
        motif = None
        is_known = False
        known_kinase = None
        
        if site_data:
            # Extract motif
            if 'SITE_+/-7_AA' in site_data and site_data['SITE_+/-7_AA']:
                motif = site_data['SITE_+/-7_AA']
            elif 'motif' in site_data and site_data['motif']:
                motif = site_data['motif']
                
            # Check if known
            if 'is_known_phosphosite' in site_data:
                is_known = site_data['is_known_phosphosite']
                
            # Check for known kinase
            for i in range(1, 6):
                kinase_field = f"KINASE_{i}"
                if kinase_field in site_data and site_data[kinase_field]:
                    known_kinase = site_data[kinase_field]
                    break
        
        # Create node
        node = {
            'id': site_id,
            'name': site_name,
            'uniprot': uniprot_id,
            'type': 'protein',  # Protein site
            'siteType': site_type or 'S',  # Default to S if not specified
            'isKnown': is_known,
            'motif': motif,
            'known_kinase': known_kinase,
            'size': 10  # Slightly larger for protein sites
        }
        
        nodes.append(node)
        node_map[site_id] = node
    
    # Then add nodes and links for the matches
    for site_id, matches in site_matches_dict.items():
        for match in matches:
            target_id = match.get('target_id')
            if not target_id:
                continue
                
            # Skip if already added
            if target_id in node_map:
                continue
                
            # Extract target info
            target_uniprot = match.get('target_uniprot')
            target_site = match.get('target_site')
            rmsd = match.get('rmsd')
            
            if not target_uniprot or not target_site or not rmsd:
                continue
            
            # Get supplementary data for this match
            target_data = supp_data_batch.get(target_id, {})
            
            motif = match.get('motif')
            known_kinase = None
            
            if target_data:
                # Extract motif if not already present
                if not motif:
                    if 'SITE_+/-7_AA' in target_data and target_data['SITE_+/-7_AA']:
                        motif = target_data['SITE_+/-7_AA']
                    elif 'motif' in target_data and target_data['motif']:
                        motif = target_data['motif']
                
                # Get known kinase
                for i in range(1, 6):
                    kinase_field = f"KINASE_{i}"
                    if kinase_field in target_data and target_data[kinase_field]:
                        known_kinase = target_data[kinase_field]
                        break
            
            # Create node for match
            node = {
                'id': target_id,
                'name': target_site,
                'uniprot': target_uniprot,
                'type': 'match',  # Match site
                'siteType': target_site[0] if target_site and target_site[0] in 'STY' else 'S',
                'isKnown': False,  # Assume false for matches
                'rmsd': rmsd,
                'motif': motif,
                'known_kinase': known_kinase,
                'size': 8  # Slightly smaller for match sites
            }
            
            nodes.append(node)
            node_map[target_id] = node
            
            # Create link from site to match
            link_id = f"{site_id}-{target_id}"
            
            # Skip if link already added
            if link_id in link_map:
                continue
                
            links.append({
                'source': site_id,
                'target': target_id,
                'rmsd': rmsd
            })
            
            link_map[link_id] = True
            # Also add reverse direction to avoid duplicate links
            link_map[f"{target_id}-{site_id}"] = True
    
    # Finally add links between matches
    for match_id, connections in match_interconnections.items():
        for connection in connections:
            target_id = connection.get('target_id')
            rmsd = connection.get('rmsd')
            
            if not target_id or not rmsd:
                continue
            
            # Skip if link already added
            link_id = f"{match_id}-{target_id}"
            if link_id in link_map:
                continue
                
            links.append({
                'source': match_id,
                'target': target_id,
                'rmsd': rmsd
            })
            
            link_map[link_id] = True
            # Also add reverse direction to avoid duplicate links
            link_map[f"{target_id}-{match_id}"] = True
    
    return {
        'nodes': nodes,
        'links': links,
        'site_matches': site_matches_dict  # Include original match data for reference
    }