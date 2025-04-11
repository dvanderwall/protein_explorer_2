"""
Cantley Kinase Scorer Module for Protein Explorer

This module analyzes sequence similarity matches to predict kinases for phosphosites
based on the Cantley kinase prediction model scores. It provides functions to:
1. Process sequence and structural similarity matches
2. Score phosphosites based on the kinase scores of their matches 
3. Create visualizations including heatmaps and rankings
4. Generate comprehensive kinase reports for proteins

The module integrates with the Protein Explorer database structure and provides
HTML output for web visualization.
"""

import logging
import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import seaborn as sns
import io
import base64
from typing import Dict, List, Optional, Tuple, Any, Union
import traceback

# Set up logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Import database functions - using relative import for Protein Explorer
try:
    from protein_explorer.db.db import (
        get_cantley_st_kinase_scores, get_cantley_st_kinase_scores_batch,
        get_cantley_y_kinase_scores, get_cantley_y_kinase_scores_batch,
        get_cantley_kinase_names, get_cantley_st_kinase_scores_by_motif
    )
except ImportError:
    logger.warning("Failed to import from protein_explorer.db.db - using direct import")
    try:
        from db.db import (
            get_cantley_st_kinase_scores, get_cantley_st_kinase_scores_batch,
            get_cantley_y_kinase_scores, get_cantley_y_kinase_scores_batch,
            get_cantley_kinase_names, get_cantley_st_kinase_scores_by_motif
        )
    except ImportError:
        logger.error("Failed to import database functions. Functionality will be limited.")

# Define kinase family mapping
# This groups kinases into families for better visualization and analysis
KINASE_FAMILIES = {
    # CMGC family
    "CMGC": [
        "CDK1", "CDK2", "CDK3", "CDK4", "CDK5", "CDK6", "CDK7", "CDK8", "CDK9", "CDK10", 
        "CDK12", "CDK13", "CDK14", "CDK16", "CDK17", "CDK18", "CDK19", "CDKL1", "CDKL5", 
        "GSK3A", "GSK3B", "MAPK1", "MAPK3", "ERK1", "ERK2", "ERK3", "ERK5", "ERK7", "JNK1", 
        "JNK2", "JNK3", "P38A", "P38B", "P38D", "P38G", "DYRK1A", "DYRK1B", "DYRK2", "DYRK3", 
        "DYRK4", "CLK1", "CLK2", "CLK3", "CLK4"
    ],
    
    # AGC family
    "AGC": [
        "AKT1", "AKT2", "AKT3", "PKA", "PKACA", "PKACB", "PKACG", "PKG1", "PKG2", "SGK1", 
        "SGK3", "p70S6K", "P70S6K", "P70S6KB", "p90RSK", "PKC", "PKCA", "PKCB", "PKCD", 
        "PKCE", "PKCG", "PKCH", "PKCI", "PKCT", "PKCZ", "ROCK1", "ROCK2", "PDK1", "RSK2", 
        "RSK3", "RSK4", "GRK1", "GRK2", "GRK3", "GRK4", "GRK5", "GRK6", "GRK7", "PRKD1", 
        "PRKD2", "PRKD3", "PKN1", "PKN2", "PKN3"
    ],
    
    # CAMK family
    "CAMK": [
        "CAMK1A", "CAMK1B", "CAMK1D", "CAMK1G", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G", 
        "CAMK4", "CAMKK1", "CAMKK2", "MARK1", "MARK2", "MARK3", "MARK4", "MELK", "AMPKA1", 
        "AMPKA2", "BRSK1", "BRSK2", "CAMLCK", "DAPK1", "DAPK2", "DAPK3", "SKMLCK", "SMMLCK"
    ],
    
    # CK1 family
    "CK1": [
        "CK1A", "CK1A2", "CK1D", "CK1E", "CK1G1", "CK1G2", "CK1G3"
    ],
    
    # STE family
    "STE": [
        "MEK1", "MEK2", "MEK5", "MKK3", "MKK4", "MKK6", "MKK7", "MAP3K15", "MEKK1", "MEKK2", 
        "MEKK3", "MEKK6", "PAK1", "PAK2", "PAK3", "PAK4", "PAK5", "PAK6", "MST1", "MST2", 
        "MST3", "MST4", "SLK", "LOK", "TAO1", "TAO2", "TAO3", "YSK1", "YSK4"
    ],
    
    # TK family (tyrosine kinases)
    "TK": [
        "ABL", "ARG", "SRC", "FYN", "YES", "LCK", "HCK", "FGR", "BLK", "LYN", "FRK", "CSK", 
        "CTK", "EGFR", "HER2", "HER4", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "PDGFRA", "PDGFRB", 
        "KIT", "CSFR", "FLT3", "VEGFR1", "VEGFR2", "VEGFR3", "TIE2", "MET", "RON", "RET", 
        "ALK", "LTK", "ROS", "TRKA", "TRKB", "TRKC", "MUSK", "AATYK", "DDR1", "DDR2", "FAK", 
        "PYK2", "ACK", "ETK", "TEC", "BTK", "ITK", "TXK", "JAK1", "JAK2", "JAK3", "TYK2", 
        "SYK", "ZAP70"
    ],
    
    # TKL family (tyrosine kinase-like)
    "TKL": [
        "RAF1", "BRAF", "ARAF", "MLK1", "MLK2", "MLK3", "MLK4", "LZK", "DLK", "ZAK", "LIMK1", 
        "LIMK2", "LRRK1", "LRRK2", "IRAK1", "IRAK4", "RIPK1", "RIPK2", "RIPK3", "BMPR1A", 
        "BMPR1B", "ALK2", "ALK4", "TGFBR1", "TGFBR2"
    ],
    
    # Atypical kinases
    "Atypical": [
        "ATM", "ATR", "DNAPK", "SMG1", "mTOR", "PDHK1", "PDHK4", "BCKDK", "GCN2", "HRI", 
        "PERK", "IRE1", "IRE2", "PINK1", "TAK1", "ASK1", "IKKA", "IKKB", "IKKE", "TBK1", 
        "NIK", "COT", "PIM1", "PIM2", "PIM3", "CK2A1", "CK2A2", "HASPIN", "WEE1", "MYT1", 
        "TLK1", "TLK2", "AurA", "AurB", "AurC", "PLK1", "PLK2", "PLK3", "PLK4", "BUB1", 
        "TTK", "CHAK1", "CHAK2", "WNK1", "WNK3", "WNK4", "NEK1", "NEK2", "NEK3", "NEK4", 
        "NEK5", "NEK6", "NEK7", "NEK8", "NEK9", "NEK10", "NEK11"
    ]
}

# Function to get the family for a given kinase
def get_kinase_family(kinase_name: str) -> str:
    """
    Determine the family of a given kinase based on predefined mappings.
    
    Args:
        kinase_name (str): Name of the kinase
        
    Returns:
        str: Family name or "Other" if not found
    """
    # Standardize kinase name by removing common prefixes/suffixes and whitespace
    clean_name = kinase_name.strip().upper()
    
    # Check each family
    for family, kinases in KINASE_FAMILIES.items():
        for k in kinases:
            # Try to account for variations in kinase naming
            if (clean_name == k or 
                clean_name == k.upper() or 
                k.upper() in clean_name or 
                (k.replace('MAPK', 'ERK') in clean_name) or
                (k.replace('ERK', 'MAPK') in clean_name)):
                return family
    
    # Some additional common mappings
    if "CK2" in clean_name:
        return "Atypical"
    if "CDK" in clean_name:
        return "CMGC"
    if "PKC" in clean_name or "PKA" in clean_name or "AKT" in clean_name:
        return "AGC"
    if "CAM" in clean_name:
        return "CAMK"
    
    return "Other"

def process_matches_for_kinase_scoring(sequence_matches: Dict[str, List[Dict]], 
                                        structural_matches: Optional[Dict[str, List[Dict]]] = None,
                                        match_type: str = "sequence") -> Dict[str, Dict]:
    """
    Process sequence or structural matches and prepare them for kinase scoring
    by categorizing matches by residue type.
    
    Args:
        sequence_matches: Dictionary of sequence matches by site ID
        structural_matches: Optional dictionary of structural matches by site ID
        match_type: Type of matches to process - "sequence", "structural", or "combined"
        
    Returns:
        Dictionary with processed match data organized by residue type
    """
    logger.info(f"Processing {match_type} matches for kinase scoring")
    
    # Initialize result dictionary
    result = {
        "proteins": set(),       # Set of unique protein UniProt IDs
        "sites": {},             # All sites from the protein
        "ST_sites": {},          # Sites of type S or T
        "Y_sites": {},           # Sites of type Y
        "ST_matches": {},        # Matches for S/T sites
        "Y_matches": {}          # Matches for Y sites
    }
    
    # Process all the original protein sites
    for site_name, matches in sequence_matches.items():
        # Extract site information
        if not matches:
            continue
            
        # Get the first match to extract query information
        first_match = matches[0] if matches else {}
        
        if not first_match or 'query_id' not in first_match:
            continue
        
        query_id = first_match.get('query_id', '')
        site_parts = query_id.split('_')
        
        if len(site_parts) < 2:
            continue
            
        protein_id = site_parts[0]
        result["proteins"].add(protein_id)
        
        # Determine site type (S, T, or Y)
        site_type = None
        
        # First try direct site_type in the match
        if 'siteType' in first_match:
            site_type = first_match['siteType']
        elif site_name and site_name[0] in ['S', 'T', 'Y']:
            site_type = site_name[0]
            
        # If still no site type, try to infer from site name
        if not site_type and site_name:
            if site_name[0] == 'S':
                site_type = 'S'
            elif site_name[0] == 'T':
                site_type = 'T'
            elif site_name[0] == 'Y':
                site_type = 'Y'
            else:
                # Default to S if unable to determine
                site_type = 'S'
        
        # Store site information
        site_info = {
            'site_name': site_name,
            'site_id': query_id,
            'protein_id': protein_id,
            'site_type': site_type,
            'motif': first_match.get('motif', ''),
            'is_known': first_match.get('is_known', False),
            'matches': []  # Will store match data
        }
        
        # Store site in appropriate category
        result["sites"][site_name] = site_info
        
        if site_type in ['S', 'T']:
            result["ST_sites"][site_name] = site_info
        elif site_type == 'Y':
            result["Y_sites"][site_name] = site_info
        
        # Process all matches for this site based on match_type
        all_matches = []
        
        # Add sequence matches
        if match_type in ["sequence", "combined"]:
            for match in matches:
                # Skip matches without target_id
                if 'target_id' not in match:
                    continue
                
                # Extract residue type - try all possible sources
                residue_type = None
                if 'ResidueType' in match:
                    residue_type = match['ResidueType']
                elif 'site_type' in match:
                    residue_type = match['site_type']
                elif 'target_site' in match and match['target_site'] and match['target_site'][0] in ['S', 'T', 'Y']:
                    residue_type = match['target_site'][0]
                else:
                    # Try to determine from target_id
                    target_parts = match.get('target_id', '').split('_')
                    if len(target_parts) >= 2:
                        target_site = target_parts[1]
                        if target_site and target_site[0] in ['S', 'T', 'Y']:
                            residue_type = target_site[0]
                    
                # If we still don't have a residue type, use the source site's type
                if not residue_type:
                    residue_type = site_type
                
                # Skip self-matches
                if match.get('target_id') == query_id:
                    continue
                    
                # Store structured match information
                match_info = {
                    'target_id': match.get('target_id', ''),
                    'target_uniprot': match.get('target_uniprot', ''),
                    'target_site': match.get('target_site', ''),
                    'residue_type': residue_type,
                    'similarity': match.get('similarity', 0),
                    'rmsd': match.get('rmsd', None),
                    'match_type': 'sequence'
                }
                
                all_matches.append(match_info)
        
        # Add structural matches if provided and requested
        if structural_matches and match_type in ["structural", "combined"] and site_name in structural_matches:
            for match in structural_matches[site_name]:
                # Skip matches without target_id
                if 'target_id' not in match:
                    continue
                
                # Extract residue type - try all possible sources
                residue_type = None
                if 'ResidueType' in match:
                    residue_type = match['ResidueType']
                elif 'site_type' in match:
                    residue_type = match['site_type']
                elif 'target_site' in match and match['target_site'] and match['target_site'][0] in ['S', 'T', 'Y']:
                    residue_type = match['target_site'][0]
                else:
                    # Try to determine from target_id
                    target_parts = match.get('target_id', '').split('_')
                    if len(target_parts) >= 2:
                        target_site = target_parts[1]
                        if target_site and target_site[0] in ['S', 'T', 'Y']:
                            residue_type = target_site[0]
                    
                # If we still don't have a residue type, use the source site's type
                if not residue_type:
                    residue_type = site_type
                    
                # Skip self-matches
                if match.get('target_id') == query_id:
                    continue
                    
                # Store structured match information
                match_info = {
                    'target_id': match.get('target_id', ''),
                    'target_uniprot': match.get('target_uniprot', ''),
                    'target_site': match.get('target_site', ''),
                    'residue_type': residue_type,
                    'similarity': match.get('similarity', 0),
                    'rmsd': match.get('rmsd', 0),
                    'match_type': 'structural'
                }
                
                all_matches.append(match_info)
        
        # Store all matches in site info
        site_info['matches'] = all_matches
        
        # Also organize matches by residue type
        st_matches = []
        y_matches = []
        
        for match in all_matches:
            if match['residue_type'] in ['S', 'T']:
                st_matches.append(match)
            elif match['residue_type'] == 'Y':
                y_matches.append(match)
        
        # Store matches in appropriate category
        if site_type in ['S', 'T']:
            result["ST_matches"][site_name] = st_matches
        elif site_type == 'Y':
            result["Y_matches"][site_name] = y_matches
    
    # Log summary statistics
    logger.info(f"Processed match data: {len(result['proteins'])} proteins, "
                f"{len(result['sites'])} sites, "
                f"{len(result['ST_sites'])} S/T sites, "
                f"{len(result['Y_sites'])} Y sites")
    
    return result

def score_site_kinases(site_name: str, matches: List[Dict], residue_type: str = "ST") -> Dict:
    """
    Score kinases for a specific site based on its matches.
    Updated to work with the Cantley tables where site IDs are in the Motif column.
    
    Args:
        site_name: Name of the site
        matches: List of match dictionaries
        residue_type: Type of residue ("ST" or "Y")
        
    Returns:
        Dictionary with kinase scores and related data
    """
    logger.info(f"Scoring kinases for {site_name} (type: {residue_type})")
    
    # If no matches, return empty results
    if not matches:
        return {
            "site_name": site_name,
            "residue_type": residue_type,
            "match_count": 0,
            "kinase_scores": {},
            "top_kinases": [],
            "kinase_families": {}
        }
    
    # Get list of target IDs
    target_ids = [match['target_id'] for match in matches if 'target_id' in match]
    
    # Skip if no valid target IDs
    if not target_ids:
        return {
            "site_name": site_name,
            "residue_type": residue_type,
            "match_count": 0,
            "kinase_scores": {},
            "top_kinases": [],
            "kinase_families": {}
        }
    
    try:
        # Get kinase scores in batch using the Motif column functions
        if residue_type == "ST":
            scores_batch = get_cantley_st_kinase_scores_batch(target_ids) or {}
        else:  # Y
            scores_batch = get_cantley_y_kinase_scores_batch(target_ids) or {}
        
        # If no scores were found, try with alternative ID formats
        if not scores_batch:
            logger.warning(f"No scores found in Motif column for {site_name}, trying alternative formats")
            
            # Create alternative format IDs to try
            alternative_ids = []
            
            for target_id in target_ids:
                # Original format: "UniProtID_ResidueNumber" (e.g., "P04637_303")
                alternative_ids.append(target_id)
                
                # Split into parts
                parts = target_id.split('_')
                if len(parts) != 2:
                    continue
                    
                uniprot_id = parts[0]
                position = parts[1]
                
                # Generate alternative formats
                # Format: Just the UniProt ID
                alternative_ids.append(uniprot_id)
                
                # Format: Just the position
                alternative_ids.append(position)
                
                # Format: UniProt-Position
                alternative_ids.append(f"{uniprot_id}-{position}")
            
            # Remove duplicates
            alternative_ids = list(set(alternative_ids))
            
            # Try batch query with alternative formats
            if residue_type == "ST":
                scores_batch = get_cantley_st_kinase_scores_batch(alternative_ids) or {}
            else:  # Y
                scores_batch = get_cantley_y_kinase_scores_batch(alternative_ids) or {}
                
            logger.info(f"Found {len(scores_batch)} scores using alternative formats")
        
        # If still no scores, look for partial matches
        if not scores_batch:
            logger.warning(f"No direct matches found for {site_name}, trying partial matches")
            
            # Try a more flexible approach: look for any entries that contain parts of our IDs
            # This would require direct SQL queries with LIKE conditions
            from protein_explorer.db.db import execute_query
            
            partial_matches = {}
            
            for target_id in target_ids[:5]:  # Limit to 5 to avoid excessive queries
                # Split ID
                parts = target_id.split('_')
                if len(parts) != 2:
                    continue
                    
                uniprot_id = parts[0]
                position = parts[1]
                
                # Try to find any entries containing the UniProt ID
                if residue_type == "ST":
                    query = f"SELECT * FROM Cantley_Kinome_Scores_All_STs WHERE Motif LIKE '%{uniprot_id}%' LIMIT 10"
                else:
                    query = f"SELECT * FROM Cantley_Kinome_Scores_All_Ys WHERE Motif LIKE '%{uniprot_id}%' LIMIT 10"
                    
                results = execute_query(query)
                
                if not results.empty:
                    for _, row in results.iterrows():
                        motif = row.get('Motif')
                        if not motif:
                            continue
                            
                        # Record this as a partial match
                        scores = {}
                        for col, val in row.items():
                            if col not in ['SiteID', 'Motif'] and pd.notna(val):
                                scores[col] = float(val)
                                
                        partial_matches[motif] = {
                            'site_id': target_id,  # Original target ID
                            'matched_motif': motif,  # What we found in the database
                            'match_type': 'partial',
                            'scores': scores
                        }
            
            # Use these partial matches if found
            if partial_matches:
                logger.info(f"Found {len(partial_matches)} partial matches")
                # Replace scores_batch with our partial matches
                scores_batch = partial_matches
        
        # If we still have no scores, return empty results
        if not scores_batch:
            logger.warning(f"No scores found for {site_name} after trying all approaches")
            return {
                "site_name": site_name,
                "residue_type": residue_type,
                "match_count": len(target_ids),
                "kinase_scores": {},
                "top_kinases": [],
                "kinase_families": {}
            }
        
        # Combine scores from all matches
        combined_scores = defaultdict(list)
        
        # For each target, extract and store scores
        for target_id in target_ids:
            # Find match info
            match_info = next((m for m in matches if m.get('target_id') == target_id), None)
            if not match_info:
                continue
                
            similarity = match_info.get('similarity', 0)
            rmsd = match_info.get('rmsd', None)
            
            # Check if we have a direct match
            target_data = scores_batch.get(target_id)
            
            # If not, look for any related match
            if not target_data:
                # Check if we have any partial matches related to this target
                for motif, data in scores_batch.items():
                    if data.get('site_id') == target_id or motif in target_id or target_id in motif:
                        target_data = data
                        break
            
            # Skip if still no data
            if not target_data or 'scores' not in target_data:
                continue
            
            # Add each kinase score, weighted by similarity if available
            for kinase, score in target_data['scores'].items():
                # Use similarity for sequence matches, inverse RMSD for structural matches
                weight = 1.0
                if similarity > 0:
                    # Higher similarity = higher weight
                    weight = similarity
                elif rmsd is not None and rmsd > 0:
                    # Lower RMSD = higher weight (max weight 1.0)
                    weight = max(0.1, min(1.0, 1.0 / rmsd))
                
                # Store the raw score, weighted score, and match details
                combined_scores[kinase].append({
                    'target_id': target_id,
                    'raw_score': score,
                    'weighted_score': score * weight,
                    'weight': weight,
                    'similarity': similarity,
                    'rmsd': rmsd
                })
        
        # Calculate aggregate scores for each kinase
        kinase_scores = {}
        
        for kinase, scores in combined_scores.items():
            # Calculate various aggregate metrics
            raw_scores = [s['raw_score'] for s in scores]
            weighted_scores = [s['weighted_score'] for s in scores]
            
            kinase_scores[kinase] = {
                'mean': np.mean(raw_scores),
                'median': np.median(raw_scores),
                'max': np.max(raw_scores),
                'weighted_mean': np.sum(weighted_scores) / np.sum([s['weight'] for s in scores]),
                'count': len(scores),
                'raw_scores': raw_scores,
                'weighted_scores': weighted_scores,
                'matches': scores
            }
        
        # Sort kinases by weighted mean score (descending)
        top_kinases = sorted(
            [(k, v['weighted_mean']) for k, v in kinase_scores.items()],
            key=lambda x: x[1],
            reverse=True
        )
        
        # Categorize kinases by family
        kinase_families = defaultdict(int)
        
        for kinase, score_data in kinase_scores.items():
            family = get_kinase_family(kinase)
            # Use weighted mean as the family score contribution
            kinase_families[family] += score_data['weighted_mean']
        
        # Normalize family scores to percentages
        total_family_score = sum(kinase_families.values())
        if total_family_score > 0:
            for family in kinase_families:
                kinase_families[family] = (kinase_families[family] / total_family_score) * 100
        
        # Return structured results
        return {
            "site_name": site_name,
            "residue_type": residue_type,
            "match_count": len(target_ids),
            "matched_count": len(combined_scores),
            "kinase_scores": kinase_scores,
            "top_kinases": top_kinases[:10],  # Get top 10 kinases
            "kinase_families": dict(kinase_families)
        }
        
    except Exception as e:
        logger.error(f"Error scoring kinases for {site_name}: {e}")
        logger.error(traceback.format_exc())
        
        # Return minimal results on error
        return {
            "site_name": site_name,
            "residue_type": residue_type,
            "match_count": len(target_ids),
            "error": str(e),
            "kinase_scores": {},
            "top_kinases": [],
            "kinase_families": {}
        }

def create_protein_kinase_report(protein_id: str, 
                                phosphosites: List[Dict], 
                                structural_matches: Optional[Dict[str, List[Dict]]] = None,
                                sequence_matches: Optional[Dict[str, List[Dict]]] = None,
                                match_type: str = "sequence",
                                include_heatmap: bool = True,
                                include_family_analysis: bool = True,
                                top_n_kinases: int = 10) -> Dict:
    """
    Create a comprehensive kinase analysis report for a protein based on its phosphosites
    and their matches.
    
    Args:
        protein_id: UniProt ID of the protein
        phosphosites: List of phosphosite dictionaries
        structural_matches: Dictionary of structural matches (optional)
        sequence_matches: Dictionary of sequence matches (optional)
        match_type: Type of matches to analyze - "sequence", "structural", or "combined"
        include_heatmap: Whether to include heatmap data in the report
        include_family_analysis: Whether to include kinase family analysis
        top_n_kinases: Number of top kinases to include in the report
        
    Returns:
        Dictionary containing the kinase report data
    """
    logger.info(f"Creating kinase report for protein {protein_id}")
    
    report = {
        "protein_id": protein_id,
        "phosphosite_count": len(phosphosites),
        "match_type": match_type,
        "site_reports": {},
        "combined_kinases": {},
        "top_kinases": [],
        "kinase_families": {},
        "has_heatmap": False,
        "heatmap_data": {}
    }
    
    # Skip if no phosphosites
    if not phosphosites:
        logger.warning(f"No phosphosites provided for protein {protein_id}")
        return report
        
    # Make sure we have matches to analyze
    if not sequence_matches and not structural_matches:
        logger.warning(f"No matches provided for protein {protein_id}")
        return report
    
    # Process matches for kinase scoring
    if sequence_matches and structural_matches and match_type == "combined":
        # Process both sequence and structural matches
        processed_matches = process_matches_for_kinase_scoring(
            sequence_matches, structural_matches, "combined"
        )
    elif sequence_matches and match_type in ["sequence", "combined"]:
        # Process only sequence matches
        processed_matches = process_matches_for_kinase_scoring(
            sequence_matches, match_type="sequence"
        )
    elif structural_matches and match_type in ["structural", "combined"]:
        # Process only structural matches (convert to expected format)
        processed_matches = process_matches_for_kinase_scoring(
            structural_matches, match_type="structural"
        )
    else:
        logger.warning(f"No matching data available for {match_type} analysis")
        return report
    
    # Store site names for reference
    site_names = list(processed_matches["sites"].keys())
    report["site_names"] = site_names
    
    # Process S/T sites
    st_site_reports = {}
    
    for site_name, site_data in processed_matches["ST_sites"].items():
        # Skip sites without matches
        if site_name not in processed_matches["ST_matches"]:
            continue
            
        matches = processed_matches["ST_matches"][site_name]
        if not matches:
            continue
            
        # Score kinases for this site
        site_report = score_site_kinases(site_name, matches, "ST")
        
        # Store in the report
        if site_report:
            st_site_reports[site_name] = site_report
    
    # Process Y sites
    y_site_reports = {}
    
    for site_name, site_data in processed_matches["Y_sites"].items():
        # Skip sites without matches
        if site_name not in processed_matches["Y_matches"]:
            continue
            
        matches = processed_matches["Y_matches"][site_name]
        if not matches:
            continue
            
        # Score kinases for this site
        site_report = score_site_kinases(site_name, matches, "Y")
        
        # Store in the report
        if site_report:
            y_site_reports[site_name] = site_report
    
    # Combine all site reports
    report["site_reports"] = {**st_site_reports, **y_site_reports}
    
    # Calculate combined kinase scores across all sites
    combined_kinases = defaultdict(list)
    
    for site_name, site_report in report["site_reports"].items():
        for kinase, score in site_report["top_kinases"]:
            combined_kinases[kinase].append(score)
    
    # Calculate average score for each kinase
    kinase_avg_scores = {}
    
    for kinase, scores in combined_kinases.items():
        kinase_avg_scores[kinase] = np.mean(scores)
    
    # Get top kinases overall
    top_kinases = sorted(
        [(k, v) for k, v in kinase_avg_scores.items()],
        key=lambda x: x[1],
        reverse=True
    )
    
    report["combined_kinases"] = kinase_avg_scores
    report["top_kinases"] = top_kinases[:top_n_kinases]
    
    # Analyze kinase families if requested
    if include_family_analysis:
        kinase_families = defaultdict(float)
        
        for kinase, score in kinase_avg_scores.items():
            family = get_kinase_family(kinase)
            kinase_families[family] += score
        
        # Normalize family scores to percentages
        total_family_score = sum(kinase_families.values())
        if total_family_score > 0:
            for family in kinase_families:
                kinase_families[family] = (kinase_families[family] / total_family_score) * 100
        
        report["kinase_families"] = dict(kinase_families)
    
    # Create heatmap data if requested
    if include_heatmap:
        try:
            heatmap_data = create_heatmap_data(report["site_reports"], top_n_kinases)
            report["has_heatmap"] = True
            report["heatmap_data"] = heatmap_data
        except Exception as e:
            logger.error(f"Error creating heatmap data: {e}")
            logger.error(traceback.format_exc())
            report["has_heatmap"] = False
    
    return report

def create_heatmap_data(site_reports: Dict[str, Dict], top_n: int = 10) -> Dict:
    """
    Create data for a heatmap visualization of top kinases across sites.
    
    Args:
        site_reports: Dictionary of site reports from create_protein_kinase_report
        top_n: Number of top kinases to include
        
    Returns:
        Dictionary containing heatmap data
    """
    # Skip if no site reports
    if not site_reports:
        return {
            "sites": [],
            "kinases": [],
            "scores": []
        }
    
    # Get the top kinases across all sites
    all_kinases = []
    for site_name, site_data in site_reports.items():
        top_kinases = site_data.get("top_kinases", [])
        all_kinases.extend([k for k, _ in top_kinases])
    
    # Count frequency of each kinase
    kinase_counts = defaultdict(int)
    for kinase in all_kinases:
        kinase_counts[kinase] += 1
    
    # Get the most frequent kinases
    top_kinases = sorted(
        [(k, v) for k, v in kinase_counts.items()],
        key=lambda x: x[1],
        reverse=True
    )[:top_n]
    
    top_kinase_names = [k for k, _ in top_kinases]
    
    # Create site list
    sites = list(site_reports.keys())
    
    # Create score matrix
    scores = []
    
    for site_name in sites:
        site_data = site_reports[site_name]
        site_scores = []
        
        # Get kinase scores for this site
        for kinase in top_kinase_names:
            # Check if kinase exists in the scores
            kinase_data = site_data.get("kinase_scores", {}).get(kinase, {})
            if kinase_data:
                # Use weighted mean score if available
                score = kinase_data.get("weighted_mean", 0)
            else:
                # Check if kinase is in top kinases
                top_kinases_dict = dict(site_data.get("top_kinases", []))
                score = top_kinases_dict.get(kinase, 0)
            
            site_scores.append(score)
        
        scores.append(site_scores)
    
    return {
        "sites": sites,
        "kinases": top_kinase_names,
        "scores": scores
    }

def create_heatmap_image(heatmap_data: Dict, title: str = "Kinase Prediction Scores") -> str:
    """
    Create a heatmap visualization as a base64-encoded image.
    
    Args:
        heatmap_data: Dictionary containing heatmap data (sites, kinases, scores)
        title: Title for the heatmap
        
    Returns:
        Base64-encoded PNG image as a data URL
    """
    if not heatmap_data or "sites" not in heatmap_data or not heatmap_data["sites"]:
        # Return empty image if no data
        plt.figure(figsize=(6, 3))
        plt.text(0.5, 0.5, "No data available", ha='center', va='center')
        plt.axis('off')
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        plt.close()
        buf.seek(0)
        img_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
        return f"data:image/png;base64,{img_base64}"
    
    # Extract data
    sites = heatmap_data["sites"]
    kinases = heatmap_data["kinases"]
    scores = heatmap_data["scores"]
    
    # Create dataframe from the scores
    df = pd.DataFrame(scores, index=sites, columns=kinases)
    
    # Create figure
    fig_height = max(4, len(sites) * 0.4 + 2)
    fig_width = max(7, len(kinases) * 0.8 + 2)
    plt.figure(figsize=(fig_width, fig_height))
    
    # Plot heatmap with color scaling based on the data range
    cmap = "YlOrRd"  # Yellow-Orange-Red colormap (increases in intensity with value)
    sns.heatmap(df, cmap=cmap, annot=True, fmt=".2f", linewidths=.5, 
                vmin=0, vmax=df.max().max(), square=False)
    
    plt.title(title, pad=20)
    plt.ylabel("Phosphosite")
    plt.xlabel("Kinase")
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    # Adjust layout to make room for labels
    plt.tight_layout()
    
    # Save figure to bytes buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
    plt.close()
    
    # Convert to base64
    buf.seek(0)
    img_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
    
    return f"data:image/png;base64,{img_base64}"

def create_kinase_family_chart(kinase_families: Dict[str, float], 
                              title: str = "Kinase Family Distribution") -> str:
    """
    Create a pie chart of kinase family distribution as a base64-encoded image.
    
    Args:
        kinase_families: Dictionary mapping family names to scores
        title: Title for the chart
        
    Returns:
        Base64-encoded PNG image as a data URL
    """
    if not kinase_families:
        # Return empty image if no data
        plt.figure(figsize=(6, 4))
        plt.text(0.5, 0.5, "No kinase family data available", ha='center', va='center')
        plt.axis('off')
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        plt.close()
        buf.seek(0)
        img_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
        return f"data:image/png;base64,{img_base64}"
    
    # Sort families by score
    sorted_families = sorted(
        [(k, v) for k, v in kinase_families.items() if v > 0],
        key=lambda x: x[1],
        reverse=True
    )
    
    # Create figure
    plt.figure(figsize=(8, 6))
    
    # Extract labels and sizes
    labels = [k for k, v in sorted_families]
    sizes = [v for k, v in sorted_families]
    
    # Create custom colormap for different kinase families
    colors = plt.cm.tab10(np.linspace(0, 1, len(labels)))
    
    # Plot pie chart
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=colors,
            wedgeprops={'linewidth': 1, 'edgecolor': 'white'})
    
    plt.title(title)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
    plt.tight_layout()
    
    # Save figure to bytes buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
    plt.close()
    
    # Convert to base64
    buf.seek(0)
    img_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
    
    return f"data:image/png;base64,{img_base64}"

def create_top_kinases_chart(top_kinases: List[Tuple[str, float]], 
                            title: str = "Top Predicted Kinases") -> str:
    """
    Create a bar chart of top predicted kinases as a base64-encoded image.
    
    Args:
        top_kinases: List of (kinase, score) tuples
        title: Title for the chart
        
    Returns:
        Base64-encoded PNG image as a data URL
    """
    if not top_kinases:
        # Return empty image if no data
        plt.figure(figsize=(6, 4))
        plt.text(0.5, 0.5, "No kinase prediction data available", ha='center', va='center')
        plt.axis('off')
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        plt.close()
        buf.seek(0)
        img_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
        return f"data:image/png;base64,{img_base64}"
    
    # Create figure
    plt.figure(figsize=(10, 6))
    
    # Extract kinases and scores
    kinases = [k for k, _ in top_kinases]
    scores = [s for _, s in top_kinases]
    
    # Determine y-axis positions
    y_pos = np.arange(len(kinases))
    
    # Create horizontal bar chart
    bars = plt.barh(y_pos, scores, align='center', alpha=0.8, 
              color=[plt.cm.tab10(i % 10) for i in range(len(kinases))])
    
    # Add values to the end of each bar
    for bar in bars:
        width = bar.get_width()
        plt.text(width + 0.01, bar.get_y() + bar.get_height()/2, 
                f'{width:.2f}', ha='left', va='center')
    
    # Add labels
    plt.yticks(y_pos, kinases)
    plt.xlabel('Score')
    plt.title(title)
    
    # Remove spines and ticks on the y-axis for cleaner look
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().tick_params(left=False)
    
    plt.tight_layout()
    
    # Save figure to bytes buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
    plt.close()
    
    # Convert to base64
    buf.seek(0)
    img_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
    
    return f"data:image/png;base64,{img_base64}"

def get_html_protein_kinase_report(kinase_report: Dict) -> str:
    """
    Generate HTML representation of a protein kinase report.
    
    Args:
        kinase_report: Kinase report dictionary
        
    Returns:
        HTML string for displaying the report
    """
    protein_id = kinase_report["protein_id"]
    phosphosite_count = kinase_report["phosphosite_count"]
    match_type = kinase_report["match_type"]
    
    # Create top kinases chart
    top_kinases_chart = ""
    if kinase_report.get("top_kinases"):
        top_kinases_img = create_top_kinases_chart(
            kinase_report["top_kinases"], 
            f"Top Predicted Kinases for {protein_id}"
        )
        top_kinases_chart = f'<div class="text-center my-3"><img src="{top_kinases_img}" class="img-fluid" alt="Top Kinases Chart"></div>'
    
    # Create kinase family distribution chart
    family_chart = ""
    if kinase_report.get("kinase_families"):
        family_img = create_kinase_family_chart(
            kinase_report["kinase_families"], 
            f"Kinase Family Distribution for {protein_id}"
        )
        family_chart = f'<div class="text-center my-3"><img src="{family_img}" class="img-fluid" alt="Kinase Family Chart"></div>'
    
    # Create heatmap image if data is available
    heatmap_html = ""
    if kinase_report.get("has_heatmap") and kinase_report.get("heatmap_data"):
        heatmap_img = create_heatmap_image(
            kinase_report["heatmap_data"], 
            f"Kinase Prediction Heatmap for {protein_id}"
        )
        heatmap_html = f'<div class="text-center my-3"><img src="{heatmap_img}" class="img-fluid" alt="Kinase Prediction Heatmap"></div>'
    
    # Create site-specific kinase tables
    site_tables = ""
    for site_name, site_report in kinase_report.get("site_reports", {}).items():
        # Skip sites without top kinases
        if not site_report.get("top_kinases"):
            continue
            
        # Create table rows for top kinases
        kinase_rows = ""
        for i, (kinase, score) in enumerate(site_report["top_kinases"][:10]):
            # Determine family and color
            family = get_kinase_family(kinase)
            # Choose a color based on the family
            color_class = {
                "CMGC": "table-primary",
                "AGC": "table-success",
                "CAMK": "table-info",
                "CK1": "table-warning",
                "STE": "table-danger",
                "TK": "table-secondary",
                "TKL": "table-light",
                "Atypical": "table-dark",
                "Other": ""
            }.get(family, "")
            
            kinase_rows += f"""
            <tr class="{color_class}">
                <td>{i+1}</td>
                <td>{kinase}</td>
                <td>{family}</td>
                <td>{score:.4f}</td>
            </tr>
            """
        
        # Create table for this site
        table = f"""
        <div class="col-md-6 mb-4">
            <div class="card h-100">
                <div class="card-header">
                    <h5 class="mb-0">{site_name} - Top Kinases</h5>
                    <small class="text-muted">
                        Based on {site_report.get("match_count", 0)} {site_report.get("residue_type", "")} phosphosite matches
                    </small>
                </div>
                <div class="card-body p-0">
                    <div class="table-responsive">
                        <table class="table table-sm table-hover table-striped mb-0">
                            <thead>
                                <tr>
                                    <th>Rank</th>
                                    <th>Kinase</th>
                                    <th>Family</th>
                                    <th>Score</th>
                                </tr>
                            </thead>
                            <tbody>
                                {kinase_rows}
                            </tbody>
                        </table>
                    </div>
                </div>
            </div>
        </div>
        """
        
        site_tables += table
    
    # Assemble the full HTML
    html = f"""
    <div class="card mb-4">
        <div class="card-header bg-primary text-white">
            <h5 class="mb-0">Kinase Prediction Analysis for {protein_id}</h5>
        </div>
        <div class="card-body">
            <p>
                This analysis predicts potential kinases for {phosphosite_count} phosphosites in {protein_id}
                based on <strong>{match_type}</strong> similarity to known phosphosites.
            </p>
            
            <div class="row">
                <div class="col-md-6">
                    {top_kinases_chart}
                </div>
                <div class="col-md-6">
                    {family_chart}
                </div>
            </div>
            
            {heatmap_html}
            
            <h4 class="mt-4 mb-3">Site-Specific Kinase Predictions</h4>
            <div class="row">
                {site_tables}
            </div>
            
            <p class="text-muted mt-3">
                <small>
                    Predictions are based on the Cantley kinase specificity model scores of
                    structurally and sequentially similar phosphosites. Higher scores indicate
                    stronger predictions.
                </small>
            </p>
        </div>
    </div>
    """
    
    return html

def get_proteins_kinase_heatmap(protein_id: str, 
                               phosphosites: List[Dict], 
                               structural_matches: Optional[Dict[str, List[Dict]]] = None,
                               sequence_matches: Optional[Dict[str, List[Dict]]] = None,
                               match_type: str = "combined") -> str:
    """
    Generate a standalone heatmap for a protein's phosphosites and their kinase predictions.
    
    Args:
        protein_id: UniProt ID of the protein
        phosphosites: List of phosphosite dictionaries
        structural_matches: Dictionary of structural matches (optional)
        sequence_matches: Dictionary of sequence matches (optional)
        match_type: Type of matches to analyze - "sequence", "structural", or "combined"
        
    Returns:
        HTML string containing the heatmap visualization
    """
    # Create kinase report
    report = create_protein_kinase_report(
        protein_id, 
        phosphosites,
        structural_matches,
        sequence_matches,
        match_type,
        include_heatmap=True
    )
    
    # Generate heatmap HTML
    heatmap_html = ""
    if report.get("has_heatmap") and report.get("heatmap_data"):
        heatmap_img = create_heatmap_image(
            report["heatmap_data"], 
            f"Kinase Prediction Heatmap for {protein_id}"
        )
        
        # Add interactive elements for the heatmap
        heatmap_html = f"""
        <div class="card mb-4">
            <div class="card-header">
                <h5 class="mb-0">Phosphosite Kinase Prediction Heatmap</h5>
            </div>
            <div class="card-body">
                <p>
                    Heatmap showing the top predicted kinases for each phosphosite in {protein_id}.
                    Higher scores (darker colors) indicate stronger predictions.
                </p>
                
                <div class="text-center my-3">
                    <img src="{heatmap_img}" class="img-fluid" alt="Kinase Prediction Heatmap">
                </div>
                
                <p class="text-muted mt-2">
                    <small>
                        Based on {match_type} similarity matches. 
                        Click on a site in the table below to see detailed predictions.
                    </small>
                </p>
            </div>
        </div>
        """
    
    return heatmap_html

def update_site_kinase_predictions(site_name: str, 
                                  sequence_matches: Dict[str, List[Dict]],
                                  structural_matches: Optional[Dict[str, List[Dict]]] = None,
                                  match_type: str = "sequence") -> Dict:
    """
    Update kinase predictions for a specific site based on its matches.
    
    Args:
        site_name: Name of the site
        sequence_matches: Dictionary of sequence matches by site ID
        structural_matches: Dictionary of structural matches by site ID (optional)
        match_type: Type of matches to analyze - "sequence", "structural", or "combined"
        
    Returns:
        Dictionary containing updated kinase predictions
    """
    # Process matches for this specific site
    site_matches = {}
    if site_name in sequence_matches:
        site_matches[site_name] = sequence_matches[site_name]
    
    # Add structural matches if provided and requested
    if structural_matches and match_type in ["structural", "combined"] and site_name in structural_matches:
        if match_type == "structural":
            # Only use structural matches
            site_matches = {site_name: structural_matches[site_name]}
        else:
            # Combined mode - add structural matches to existing sequence matches
            if site_name not in site_matches:
                site_matches[site_name] = []
            site_matches[site_name].extend(structural_matches[site_name])
    
    # Process matches for kinase scoring
    processed_matches = process_matches_for_kinase_scoring(site_matches, match_type=match_type)
    
    # Determine residue type
    
    
    # Score kinases for this site
    st_matches = processed_matches.get("ST_matches", {}).get(site_name, [])
    site_report = score_site_kinases(site_name, st_matches, "ST")
        
    
    return site_report