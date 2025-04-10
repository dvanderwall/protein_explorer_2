"""
Cantley_Kinase_Scorer module for analyzing kinase scores based on phosphosite matches.

This module provides functions to:
1. Query Cantley kinase score tables for phosphosites
2. Aggregate scores across structurally similar sites
3. Generate heatmaps of kinase scores for similar sites
4. Rank kinases for specific phosphosites based on pooled scores
"""

import os
import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Set, Optional, Union, Any
import json
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO
import base64

# Import database functions
try:
    from protein_explorer.db.db import (
        get_cantley_st_kinase_scores, get_cantley_st_kinase_scores_batch,
        get_cantley_y_kinase_scores, get_cantley_y_kinase_scores_batch,
        get_cantley_kinase_names
    )
except ImportError:
    # Fallback import path for different project structures
    from db.db import (
        get_cantley_st_kinase_scores, get_cantley_st_kinase_scores_batch,
        get_cantley_y_kinase_scores, get_cantley_y_kinase_scores_batch,
        get_cantley_kinase_names
    )

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Define kinase family mapping
KINASE_FAMILIES = {
    # CMGC group
    "CDK": ["CDK1", "CDK2", "CDK3", "CDK4", "CDK5", "CDK6", "CDK7", "CDK8", "CDK9", 
            "CDK10", "CDK12", "CDK13", "CDK14", "CDK16", "CDK17", "CDK18", "CDK19"],
    "MAPK": ["ERK1", "ERK2", "ERK3", "ERK5", "ERK7", "JNK1", "JNK2", "JNK3", "P38A", 
             "P38B", "P38D", "P38G"],
    "GSK": ["GSK3A", "GSK3B"],
    "CLK": ["CLK1", "CLK2", "CLK3", "CLK4"],
    "DYRK": ["DYRK1A", "DYRK1B", "DYRK2", "DYRK3", "DYRK4"],
    "HIPK": ["HIPK1", "HIPK2", "HIPK3", "HIPK4"],
    
    # AGC group
    "PKA": ["PKACA", "PKACB", "PKACG"],
    "PKC": ["PKCA", "PKCB", "PKCD", "PKCE", "PKCG", "PKCH", "PKCI", "PKCT", "PKCZ"],
    "PKG": ["PKG1", "PKG2"],
    "PKN": ["PKN1", "PKN2", "PKN3"],
    "RSK": ["p90RSK", "RSK2", "RSK3", "RSK4"],
    "S6K": ["p70S6K", "P70S6KB"],
    "SGK": ["SGK1", "SGK3"],
    "AKT": ["AKT1", "AKT2", "AKT3"],
    "ROCK": ["ROCK1", "ROCK2"],
    "PRKD": ["PRKD1", "PRKD2", "PRKD3"],
    "MAST": ["MASTL"],
    
    # CAMK group
    "CAMK": ["CAMK1A", "CAMK1B", "CAMK1D", "CAMK1G", "CAMK2A", "CAMK2B", "CAMK2D", 
             "CAMK2G", "CAMK4", "CAMKK1", "CAMKK2", "CAMLCK", "smMLCK", "SKMLCK"],
    "MARK": ["MARK1", "MARK2", "MARK3", "MARK4"],
    "DAPK": ["DAPK1", "DAPK2", "DAPK3"],
    "AMPK": ["AMPKA1", "AMPKA2"],
    "BRSK": ["BRSK1", "BRSK2"],
    "NUAK": ["NUAK1", "NUAK2"],
    
    # CK1 group
    "CK1": ["CK1A", "CK1A2", "CK1D", "CK1E", "CK1G1", "CK1G2", "CK1G3"],
    
    # STE group
    "STE7": ["MEK1", "MEK2", "MEK5"],
    "STE11": ["MEKK1", "MEKK2", "MEKK3", "MEKK6"],
    "STE20": ["MST1", "MST2", "MST3", "MST4", "PAK1", "PAK2", "PAK3", "PAK4", "PAK5", "PAK6"],
    
    # Other groups
    "CK2": ["CK2A1", "CK2A2"],
    "PLK": ["PLK1", "PLK2", "PLK3", "PLK4"],
    "ATM": ["ATM", "ATR", "DNAPK"],
    "AURK": ["AURA", "AURB", "AURC"],
    "PIKK": ["MTOR", "SMG1"],
    "NEK": ["NEK1", "NEK2", "NEK3", "NEK4", "NEK5", "NEK6", "NEK7", "NEK8", "NEK9", "NEK10_TYR", "NEK11"],
    "IKK": ["IKKA", "IKKB", "IKKE"],
    "GRK": ["GRK1", "GRK2", "GRK3", "GRK4", "GRK5", "GRK6", "GRK7"],
    "PIM": ["PIM1", "PIM2", "PIM3"],
    "IRAK": ["IRAK1", "IRAK4"],
    
    # Tyrosine kinase families
    "SRC": ["SRC", "LCK", "FYN", "YES", "BLK", "HCK", "FGR", "LYN"],
    "ABL": ["ABL", "ARG"],
    "EGFR": ["EGFR", "HER2", "HER4"],
    "FGFR": ["FGFR1", "FGFR2", "FGFR3", "FGFR4"],
    "JAK": ["JAK1", "JAK2", "JAK3", "TYK2"],
    "EPH": ["EPHA1", "EPHA2", "EPHA3", "EPHA4", "EPHA5", "EPHA6", "EPHA7", "EPHA8", 
            "EPHB1", "EPHB2", "EPHB3", "EPHB4"],
    "PDGFR": ["PDGFRA", "PDGFRB"],
    "VEGFR": ["VEGFR1", "VEGFR2", "VEGFR3"],
    "TRK": ["TRKA", "TRKB", "TRKC"],
    "TEC": ["TEC", "BTK", "ITK", "ETK", "TXK"]
}

# Reverse mapping to find family for any kinase
KINASE_TO_FAMILY = {}
for family, kinase_list in KINASE_FAMILIES.items():
    for kinase in kinase_list:
        KINASE_TO_FAMILY[kinase] = family

def get_residue_type(site_name: str) -> str:
    """
    Determine the residue type (S, T, or Y) from a site name.
    
    Args:
        site_name: Site name like "S123", "T45", or Y27, or "P12345_123"
    
    Returns:
        Single letter residue type code: "S", "T", or "Y"
    """
    logger.info(f"[DEBUG get_residue_type] Determining residue type for site: {site_name}")
    
    # Handle invalid input
    if not site_name or not isinstance(site_name, str):
        logger.error(f"[DEBUG get_residue_type] Invalid site_name: {site_name} (type: {type(site_name)})")
        return "S"  # Default to S for invalid input
    
    # Handle UniProt_ResidueNumber format
    if "_" in site_name:
        # Check if we have a residue type in the second part
        parts = site_name.split("_")
        logger.info(f"[DEBUG get_residue_type] Site has underscore format. Parts: {parts}")
        if len(parts) >= 2:
            if parts[1] and parts[1][0] in "STY":
                logger.info(f"[DEBUG get_residue_type] Found residue type {parts[1][0]} in second part")
                return parts[1][0]
            # If no residue type in second part, default to S (most common)
            logger.info(f"[DEBUG get_residue_type] No residue type in second part, defaulting to S")
            return "S"
    
    # Handle SiteType+Number format (e.g., "S123")
    if site_name and site_name[0] in "STY":
        logger.info(f"[DEBUG get_residue_type] Found residue type {site_name[0]} at start of site name")
        return site_name[0]
    
    # Default to serine if type can't be determined
    logger.warning(f"[DEBUG get_residue_type] Could not determine residue type for {site_name}, defaulting to S")
    return "S"

def get_cantley_kinase_scores(site_id: str) -> Optional[Dict]:
    """
    Get Cantley kinase scores for a specific phosphosite, automatically choosing
    the appropriate table based on residue type.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary with kinase scores or None if not found
    """
    logger.info(f"[DEBUG get_cantley_kinase_scores] Getting kinase scores for site_id: {site_id}")
    
    # Determine residue type
    residue_type = get_residue_type(site_id)
    logger.info(f"[DEBUG get_cantley_kinase_scores] Determined residue type: {residue_type}")
    
    # Get scores from appropriate table
    if residue_type == "Y":
        logger.info(f"[DEBUG get_cantley_kinase_scores] Using Y kinase table for {site_id}")
        result = get_cantley_y_kinase_scores(site_id)
    else:  # S or T
        logger.info(f"[DEBUG get_cantley_kinase_scores] Using S/T kinase table for {site_id}")
        result = get_cantley_st_kinase_scores(site_id)
    
    # Log results summary
    if result is None:
        logger.warning(f"[DEBUG get_cantley_kinase_scores] No scores found for {site_id}")
    else:
        score_count = len(result.get('scores', {}))
        logger.info(f"[DEBUG get_cantley_kinase_scores] Found {score_count} kinase scores for {site_id}")
    
    return result

def get_cantley_kinase_scores_batch(site_ids: List[str]) -> Dict[str, Dict]:
    """
    Get Cantley kinase scores for multiple phosphosites in a batch operation,
    automatically choosing the appropriate table based on residue type.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary mapping site IDs to score dictionaries
    """
    logger.info(f"[DEBUG get_cantley_kinase_scores_batch] Getting kinase scores for {len(site_ids)} sites")
    
    # Check for valid input
    if not site_ids:
        logger.warning("[DEBUG get_cantley_kinase_scores_batch] Empty site_ids list provided")
        return {}
    
    # Log the first few site IDs for debugging
    sample_ids = site_ids[:5]
    logger.info(f"[DEBUG get_cantley_kinase_scores_batch] Sample site IDs: {sample_ids}")
    
    # Split site IDs by residue type
    st_site_ids = []
    y_site_ids = []
    
    for site_id in site_ids:
        try:
            residue_type = get_residue_type(site_id)
            if residue_type == "Y":
                y_site_ids.append(site_id)
            else:  # S or T
                st_site_ids.append(site_id)
        except Exception as e:
            logger.error(f"[DEBUG get_cantley_kinase_scores_batch] Error determining residue type for {site_id}: {e}")
    
    logger.info(f"[DEBUG get_cantley_kinase_scores_batch] Split into {len(st_site_ids)} S/T sites and {len(y_site_ids)} Y sites")
    
    # Get scores from each table
    results = {}
    
    try:
        if st_site_ids:
            logger.info(f"[DEBUG get_cantley_kinase_scores_batch] Querying S/T table for {len(st_site_ids)} sites")
            st_results = get_cantley_st_kinase_scores_batch(st_site_ids)
            logger.info(f"[DEBUG get_cantley_kinase_scores_batch] Got results for {len(st_results)} S/T sites")
            results.update(st_results)
        
        if y_site_ids:
            logger.info(f"[DEBUG get_cantley_kinase_scores_batch] Querying Y table for {len(y_site_ids)} sites")
            y_results = get_cantley_y_kinase_scores_batch(y_site_ids)
            logger.info(f"[DEBUG get_cantley_kinase_scores_batch] Got results for {len(y_results)} Y sites")
            results.update(y_results)
    except Exception as e:
        logger.error(f"[DEBUG get_cantley_kinase_scores_batch] Error getting batch scores: {e}")
        import traceback
        logger.error(traceback.format_exc())
    
    logger.info(f"[DEBUG get_cantley_kinase_scores_batch] Returning results for {len(results)} sites")
    return results

def aggregate_kinase_scores(
    source_site_id: str, 
    match_site_ids: List[str], 
    similarity_weights: Optional[Dict[str, float]] = None,
    min_similarity: float = 0.0
) -> Dict[str, float]:
    """
    Aggregate kinase scores across multiple similar sites, with optional weighting
    based on similarity scores.
    
    Args:
        source_site_id: Source site ID in format 'UniProtID_ResidueNumber'
        match_site_ids: List of target site IDs in format 'UniProtID_ResidueNumber'
        similarity_weights: Optional dict mapping site IDs to similarity weights (0-1)
        min_similarity: Minimum similarity value to include a site (default: 0.0)
        
    Returns:
        Dictionary mapping kinase names to aggregated scores
    """
    logger.info(f"[DEBUG aggregate_kinase_scores] Aggregating scores for source: {source_site_id} and {len(match_site_ids) if match_site_ids else 0} matches")
    
    # Check for valid match_site_ids
    if match_site_ids is None:
        match_site_ids = []
        logger.info(f"[DEBUG aggregate_kinase_scores] match_site_ids is None, using empty list")
    
    # Initialize similarity_weights if not provided
    if similarity_weights is None:
        similarity_weights = {}
        if match_site_ids:
            logger.info(f"[DEBUG aggregate_kinase_scores] No similarity_weights provided, defaulting to 0.5 for all matches")
            similarity_weights = {site_id: 0.5 for site_id in match_site_ids}
    
    # Log min_similarity threshold
    logger.info(f"[DEBUG aggregate_kinase_scores] Using min_similarity threshold: {min_similarity}")
    
    # Get all site IDs including the source
    all_site_ids = [source_site_id] + match_site_ids
    
    # Get kinase scores for all sites
    logger.info(f"[DEBUG aggregate_kinase_scores] Getting kinase scores for {len(all_site_ids)} sites")
    all_scores = get_cantley_kinase_scores_batch(all_site_ids)
    logger.info(f"[DEBUG aggregate_kinase_scores] Got scores for {len(all_scores)} sites")
    
    # Determine whether we're working with S/T or Y sites
    residue_type = get_residue_type(source_site_id)
    logger.info(f"[DEBUG aggregate_kinase_scores] Source site residue type: {residue_type}")
    
    # Get the list of all possible kinases based on residue type
    try:
        kinase_names = get_cantley_kinase_names(residue_type=residue_type)
        logger.info(f"[DEBUG aggregate_kinase_scores] Retrieved kinase names by residue type")
        
        if residue_type == "Y":
            all_kinases = kinase_names["Y_kinases"]
            logger.info(f"[DEBUG aggregate_kinase_scores] Using Y_kinases ({len(all_kinases)} kinases)")
        else:  # S or T
            all_kinases = kinase_names["S/T_kinases"]
            logger.info(f"[DEBUG aggregate_kinase_scores] Using S/T_kinases ({len(all_kinases)} kinases)")
    except Exception as e:
        logger.error(f"[DEBUG aggregate_kinase_scores] Error getting kinase names: {e}")
        logger.info("[DEBUG aggregate_kinase_scores] Falling back to empty list of kinases")
        all_kinases = []
    
    # Initialize result dictionary with zeros for all kinases
    aggregated_scores = {kinase: 0.0 for kinase in all_kinases}
    
    # Track how many sites contributed to each kinase score
    contribution_counts = {kinase: 0 for kinase in all_kinases}
    
    # Process each site's scores
    for site_id in all_site_ids:
        # Get site scores
        site_data = all_scores.get(site_id)
        if not site_data or 'scores' not in site_data:
            logger.info(f"[DEBUG aggregate_kinase_scores] No scores found for {site_id}")
            continue
        
        # Get similarity weight for this site (1.0 for source site)
        weight = 1.0 if site_id == source_site_id else similarity_weights.get(site_id, 0.5)
        logger.info(f"[DEBUG aggregate_kinase_scores] Site {site_id} has weight: {weight}")
        
        # Skip sites below the minimum similarity threshold
        if site_id != source_site_id and weight < min_similarity:
            logger.info(f"[DEBUG aggregate_kinase_scores] Skipping {site_id} - weight {weight} below threshold {min_similarity}")
            continue
        
        # Add weighted scores
        score_count = 0
        for kinase, score in site_data['scores'].items():
            if kinase in aggregated_scores:
                # Only count non-zero scores to avoid diluting the aggregation
                if score > 0:
                    aggregated_scores[kinase] += score * weight
                    contribution_counts[kinase] += 1
                    score_count += 1
        
        logger.info(f"[DEBUG aggregate_kinase_scores] Added {score_count} non-zero scores from {site_id}")
    
    # Normalize scores by the number of contributing sites
    normalized_count = 0
    for kinase in aggregated_scores:
        if contribution_counts[kinase] > 0:
            # Calculate average score, weighted by similarity
            aggregated_scores[kinase] /= contribution_counts[kinase]
            normalized_count += 1
        else:
            # No contributions, set to zero
            aggregated_scores[kinase] = 0.0
    
    logger.info(f"[DEBUG aggregate_kinase_scores] Normalized {normalized_count} kinase scores")
    
    # Log number of non-zero scores in result
    non_zero_count = sum(1 for score in aggregated_scores.values() if score > 0)
    logger.info(f"[DEBUG aggregate_kinase_scores] Final result has {non_zero_count} non-zero scores out of {len(aggregated_scores)} kinases")
    
    return aggregated_scores

def create_kinase_heatmap(
    site_ids: List[str],
    top_n_kinases: int = 20,
    score_type: str = "auto",
    similarity_weights: Optional[Dict[str, float]] = None,
    color_map: str = "viridis"
) -> Dict:
    """
    Create a heatmap visualization of kinase scores for multiple sites.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        top_n_kinases: Number of top kinases to include in heatmap (default: 20)
        score_type: "auto", "ST", or "Y" - determines which score table to use
        similarity_weights: Optional dict mapping site IDs to similarity weights
        color_map: Matplotlib colormap name to use
        
    Returns:
        Dictionary with heatmap data and top kinases
    """
    logger.info(f"[DEBUG create_kinase_heatmap] Creating heatmap for {len(site_ids) if site_ids else 0} sites")
    logger.info(f"[DEBUG create_kinase_heatmap] Parameters: top_n_kinases={top_n_kinases}, score_type={score_type}")
    
    if not site_ids:
        logger.warning("[DEBUG create_kinase_heatmap] No site IDs provided")
        return {"error": "No site IDs provided"}
    
    # Log a sample of the site IDs
    if site_ids:
        sample = site_ids[:5]
        logger.info(f"[DEBUG create_kinase_heatmap] Sample site IDs: {sample}")
    
    # Initialize similarity_weights if not provided
    if similarity_weights is None:
        logger.info("[DEBUG create_kinase_heatmap] No similarity weights provided, using default weight of 1.0")
        similarity_weights = {}
    
    # Determine residue type if auto mode
    try:
        if score_type == "auto":
            # Use type of first site
            residue_type = get_residue_type(site_ids[0])
            if residue_type == "Y":
                score_type = "Y"
                logger.info("[DEBUG create_kinase_heatmap] Auto mode determined Y residue type")
            else:
                score_type = "ST"
                logger.info("[DEBUG create_kinase_heatmap] Auto mode determined S/T residue type")
    except Exception as e:
        logger.error(f"[DEBUG create_kinase_heatmap] Error determining residue type: {e}")
        score_type = "ST"  # Default to ST if there's an error
        logger.info("[DEBUG create_kinase_heatmap] Defaulting to S/T score type after error")
    
    # Get scores for all sites
    logger.info(f"[DEBUG create_kinase_heatmap] Getting scores for {len(site_ids)} sites with score_type: {score_type}")
    all_scores = {}
    try:
        if score_type == "Y":
            all_scores = get_cantley_y_kinase_scores_batch(site_ids)
        else:  # ST
            all_scores = get_cantley_st_kinase_scores_batch(site_ids)
        
        logger.info(f"[DEBUG create_kinase_heatmap] Retrieved scores for {len(all_scores)} sites")
    except Exception as e:
        logger.error(f"[DEBUG create_kinase_heatmap] Error getting scores: {e}")
        return {"error": f"Error retrieving kinase scores: {str(e)}"}
    
    # Extract scores for each site
    site_scores = {}
    all_kinases = set()
    
    for site_id, data in all_scores.items():
        if data and 'scores' in data:
            logger.info(f"[DEBUG create_kinase_heatmap] Processing scores for site: {site_id}")
            # Get similarity weight or default to 1.0
            weight = 1.0
            if similarity_weights and site_id in similarity_weights:
                weight = similarity_weights[site_id]
                logger.info(f"[DEBUG create_kinase_heatmap] Using similarity weight: {weight} for site: {site_id}")
            
            # Extract scores and add to site_scores
            site_kinase_count = len(data['scores'])
            logger.info(f"[DEBUG create_kinase_heatmap] Site {site_id} has {site_kinase_count} kinase scores")
            
            site_scores[site_id] = {k: v * weight for k, v in data['scores'].items()}
            # Track all kinases
            all_kinases.update(data['scores'].keys())
    
    # Log the number of unique kinases found
    logger.info(f"[DEBUG create_kinase_heatmap] Total unique kinases across all sites: {len(all_kinases)}")
    
    # If no scores found, return error
    if not site_scores:
        logger.warning("[DEBUG create_kinase_heatmap] No kinase scores found for any sites")
        return {"error": "No kinase scores found for the provided sites"}
    
    # Create a matrix of all scores
    logger.info("[DEBUG create_kinase_heatmap] Creating score matrix")
    score_matrix = []
    ordered_sites = []
    ordered_site_labels = []
    
    for site_id, scores in site_scores.items():
        row = []
        for kinase in all_kinases:
            row.append(scores.get(kinase, 0.0))
        
        score_matrix.append(row)
        ordered_sites.append(site_id)
        
        # Create shorter label with UniProt ID and site
        if "_" in site_id:
            parts = site_id.split("_")
            if len(parts) >= 2:
                uniprot = parts[0]
                residue = parts[1]
                # Extract site type if not present in residue
                if not residue[0].isalpha():
                    residue_type = "S" if score_type == "ST" else "Y"
                    residue = f"{residue_type}{residue}"
                ordered_site_labels.append(f"{uniprot}_{residue}")
                logger.info(f"[DEBUG create_kinase_heatmap] Created site label: {uniprot}_{residue} for {site_id}")
            else:
                ordered_site_labels.append(site_id)
        else:
            ordered_site_labels.append(site_id)
    
    # Convert to numpy array
    try:
        logger.info(f"[DEBUG create_kinase_heatmap] Converting matrix of size {len(score_matrix)}x{len(all_kinases)} to numpy array")
        score_array = np.array(score_matrix)
        
        # Calculate sum of scores for each kinase
        kinase_totals = score_array.sum(axis=0)
        
        # Convert kinases to list and sort by total score
        kinase_list = list(all_kinases)
        kinase_scores = [(kinase, kinase_totals[i]) for i, kinase in enumerate(kinase_list)]
        kinase_scores.sort(key=lambda x: x[1], reverse=True)
        
        # Get top N kinases
        top_kinases = [k for k, s in kinase_scores[:top_n_kinases]]
        top_kinase_scores = [s for k, s in kinase_scores[:top_n_kinases]]
        
        logger.info(f"[DEBUG create_kinase_heatmap] Top kinase: {top_kinases[0] if top_kinases else 'None'}")
        
        # Extract just the columns for top kinases
        top_indices = [kinase_list.index(k) for k in top_kinases]
        top_score_array = score_array[:, top_indices]
        
        # Log the shape of the final score array
        logger.info(f"[DEBUG create_kinase_heatmap] Final score array shape: {top_score_array.shape}")
        
        # Create heatmap data structure
        heatmap_data = {
            "sites": ordered_site_labels,
            "site_ids": ordered_sites,
            "kinases": top_kinases, 
            "scores": top_score_array.tolist(),
            "top_kinase_scores": top_kinase_scores,
            "all_kinase_scores": kinase_scores[:min(100, len(kinase_scores))],  # Limit to top 100 to avoid huge payloads
            "score_type": score_type
        }
        
        logger.info(f"[DEBUG create_kinase_heatmap] Successfully created heatmap data with {len(ordered_sites)} sites and {len(top_kinases)} kinases")
        return heatmap_data
        
    except Exception as e:
        logger.error(f"[DEBUG create_kinase_heatmap] Error creating heatmap data: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return {"error": f"Error creating heatmap data: {str(e)}"}


def render_kinase_heatmap(
    heatmap_data: Dict, 
    title: str = "Kinase Score Heatmap",
    color_map: str = "viridis",
    fig_width: int = 10,
    fig_height: int = 8,
    as_base64: bool = True
) -> Union[plt.Figure, str]:
    """
    Render a kinase score heatmap as a figure or base64 encoded image.
    
    Args:
        heatmap_data: Heatmap data dict from create_kinase_heatmap
        title: Title for the heatmap
        color_map: Matplotlib colormap name
        fig_width: Figure width in inches
        fig_height: Figure height in inches
        as_base64: If True, return base64 encoded PNG, otherwise return Figure
        
    Returns:
        Matplotlib Figure object or base64 encoded PNG string
    """
    if "error" in heatmap_data:
        logger.error(f"Cannot render heatmap: {heatmap_data['error']}")
        return "" if as_base64 else None
    
    # Extract data
    sites = heatmap_data["sites"]
    kinases = heatmap_data["kinases"]
    scores = heatmap_data["scores"]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        scores, 
        annot=False, 
        cmap=color_map,
        xticklabels=kinases,
        yticklabels=sites,
        ax=ax
    )
    
    # Add title
    ax.set_title(title)
    
    # Set labels
    ax.set_xlabel("Kinases")
    ax.set_ylabel("Phosphosites")
    
    # Rotate x labels for better readability
    plt.xticks(rotation=90)
    
    # Adjust layout
    plt.tight_layout()
    
    if as_base64:
        # Convert to base64 encoded PNG
        buf = BytesIO()
        plt.savefig(buf, format='png', dpi=100)
        buf.seek(0)
        img_str = base64.b64encode(buf.read()).decode('utf-8')
        plt.close(fig)
        return img_str
    else:
        return fig

def predict_top_kinases(
    source_site_id: str,
    match_site_ids: List[str] = None,
    similarity_weights: Dict[str, float] = None,
    min_similarity: float = 0.0,
    top_n: int = 10
) -> List[Dict[str, Any]]:
    """
    Predict the top kinases for a phosphosite based on its own scores
    and those of structurally or sequence similar sites.
    
    Args:
        source_site_id: Source site ID in format 'UniProtID_ResidueNumber'
        match_site_ids: List of target site IDs to consider (if None, use only source)
        similarity_weights: Dict mapping site IDs to similarity weights (0-1)
        min_similarity: Minimum similarity value to include a site
        top_n: Number of top kinases to return
        
    Returns:
        List of dict with kinase info: {'kinase': name, 'score': score, 'family': family}
    """
    logger.info(f"[DEBUG predict_top_kinases] Predicting top kinases for site: {source_site_id}")
    logger.info(f"[DEBUG predict_top_kinases] Match count: {len(match_site_ids) if match_site_ids else 0}")
    logger.info(f"[DEBUG predict_top_kinases] Parameters: min_similarity={min_similarity}, top_n={top_n}")
    
    # If no match_site_ids provided, use only the source site
    if match_site_ids is None:
        match_site_ids = []
        logger.info("[DEBUG predict_top_kinases] No match_site_ids provided, using only source site")
    
    # If no similarity_weights provided, use default of 0.5 for all matches
    if similarity_weights is None and match_site_ids:
        logger.info("[DEBUG predict_top_kinases] No similarity_weights provided, defaulting to 0.5 for all matches")
        similarity_weights = {site_id: 0.5 for site_id in match_site_ids}
    
    # Log first few matches and weights for debugging
    if match_site_ids and similarity_weights:
        sample_matches = match_site_ids[:3]
        sample_weights = {site_id: similarity_weights.get(site_id, 0.5) for site_id in sample_matches}
        logger.info(f"[DEBUG predict_top_kinases] Sample matches and weights: {sample_weights}")
    
    # Get aggregated scores
    try:
        logger.info(f"[DEBUG predict_top_kinases] Aggregating scores for {source_site_id} with {len(match_site_ids)} matches")
        aggregated_scores = aggregate_kinase_scores(
            source_site_id, 
            match_site_ids, 
            similarity_weights,
            min_similarity
        )
        
        # Log number of scores
        score_count = sum(1 for score in aggregated_scores.values() if score > 0)
        logger.info(f"[DEBUG predict_top_kinases] Aggregated {score_count} non-zero scores")
        
        # Sort kinases by score
        sorted_kinases = sorted(
            [(k, s) for k, s in aggregated_scores.items() if s > 0], 
            key=lambda x: x[1], 
            reverse=True
        )
        
        logger.info(f"[DEBUG predict_top_kinases] Sorted {len(sorted_kinases)} kinases by score")
        
        # Create result list with family information
        result = []
        for kinase, score in sorted_kinases[:top_n]:
            family = KINASE_TO_FAMILY.get(kinase, "Other")
            result.append({
                "kinase": kinase,
                "score": score,
                "family": family
            })
        
        logger.info(f"[DEBUG predict_top_kinases] Returning {len(result)} top kinases")
        if result:
            top_kinase = result[0]['kinase']
            top_score = result[0]['score']
            logger.info(f"[DEBUG predict_top_kinases] Top kinase: {top_kinase} with score: {top_score}")
        
        return result
    
    except Exception as e:
        logger.error(f"[DEBUG predict_top_kinases] Error predicting top kinases: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return []  # Return empty list on error


def get_kinase_family_distribution(
    source_site_id: str,
    match_site_ids: List[str] = None,
    similarity_weights: Dict[str, float] = None,
    min_similarity: float = 0.0,
    top_n_families: int = 10
) -> Dict[str, float]:
    """
    Calculate kinase family distribution for a phosphosite, based on aggregated scores.
    
    Args:
        source_site_id: Source site ID in format 'UniProtID_ResidueNumber'
        match_site_ids: List of target site IDs to consider (if None, use only source)
        similarity_weights: Dict mapping site IDs to similarity weights (0-1)
        min_similarity: Minimum similarity value to include a site
        top_n_families: Number of top kinase families to return
        
    Returns:
        Dict mapping family names to scores
    """
    # Get all kinase scores
    aggregated_scores = aggregate_kinase_scores(
        source_site_id, 
        match_site_ids or [], 
        similarity_weights,
        min_similarity
    )
    
    # Create family-to-score mapping
    family_scores = {}
    
    for kinase, score in aggregated_scores.items():
        if score <= 0:
            continue
            
        family = KINASE_TO_FAMILY.get(kinase, "Other")
        if family in family_scores:
            family_scores[family] += score
        else:
            family_scores[family] = score
    
    # Sort families by score
    sorted_families = sorted(
        [(fam, score) for fam, score in family_scores.items()],
        key=lambda x: x[1],
        reverse=True
    )
    
    # Take top N families
    top_families = sorted_families[:top_n_families]
    
    # Return as dict
    return {fam: score for fam, score in top_families}

def create_family_pie_chart(
    family_scores: Dict[str, float],
    title: str = "Kinase Family Distribution",
    as_base64: bool = True
) -> Union[plt.Figure, str]:
    """
    Create a pie chart visualization of kinase family distribution.
    
    Args:
        family_scores: Dict mapping family names to scores
        title: Title for the chart
        as_base64: If True, return base64 encoded PNG, otherwise return Figure
        
    Returns:
        Matplotlib Figure object or base64 encoded PNG string
    """
    if not family_scores:
        logger.error("Cannot create pie chart: No family scores provided")
        return "" if as_base64 else None
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Create pie chart
    wedges, texts, autotexts = ax.pie(
        family_scores.values(), 
        labels=family_scores.keys(),
        autopct='%1.1f%%',
        startangle=90, 
        shadow=False
    )
    
    # Enhance text visibility
    for text in texts:
        text.set_fontsize(10)
    for autotext in autotexts:
        autotext.set_fontsize(9)
        autotext.set_color('white')
    
    # Add title
    ax.set_title(title)
    
    # Make circular
    ax.axis('equal')
    
    # Add legend outside of pie
    ax.legend(
        wedges, family_scores.keys(),
        title="Kinase Families",
        loc="center left",
        bbox_to_anchor=(1, 0, 0.5, 1)
    )
    
    # Adjust layout
    plt.tight_layout()
    
    if as_base64:
        # Convert to base64 encoded PNG
        buf = BytesIO()
        plt.savefig(buf, format='png', dpi=100)
        buf.seek(0)
        img_str = base64.b64encode(buf.read()).decode('utf-8')
        plt.close(fig)
        return img_str
    else:
        return fig

def convert_to_html_heatmap(heatmap_data: Dict) -> str:
    """
    Generate an HTML visualization of a kinase score heatmap without using matplotlib.
    This is useful for web-based rendering where direct D3.js or HTML/CSS visualization
    is preferred over static images.
    
    Args:
        heatmap_data: Heatmap data dict from create_kinase_heatmap
        
    Returns:
        HTML string with interactive heatmap visualization
    """
    if "error" in heatmap_data:
        return f'<div class="alert alert-warning">Error: {heatmap_data["error"]}</div>'
    
    # Extract data
    sites = heatmap_data["sites"]
    kinases = heatmap_data["kinases"]
    scores = heatmap_data["scores"]
    
    # Serialize data for JavaScript
    scores_json = json.dumps(scores)
    sites_json = json.dumps(sites)
    kinases_json = json.dumps(kinases)
    
    html = f"""
    <div class="card mb-4">
        <div class="card-header">
            <h5 class="mb-0">Kinase Score Heatmap ({heatmap_data["score_type"]} kinases)</h5>
        </div>
        <div class="card-body">
            <div id="kinase-heatmap-container" style="height: 500px; width: 100%; position: relative;"></div>
            
            <p class="text-muted mt-3 mb-0">
                <small>
                    This heatmap shows kinase prediction scores for each phosphosite (rows) and kinase (columns).
                    Darker colors indicate higher prediction scores.
                </small>
            </p>
        </div>
    </div>
    
    <script>
    document.addEventListener('DOMContentLoaded', function() {{
        // Get data from backend
        const sites = {sites_json};
        const kinases = {kinases_json};
        const scores = {scores_json};
        
        // Set up dimensions
        const margin = {{top: 50, right: 70, bottom: 150, left: 120}};
        const width = document.getElementById('kinase-heatmap-container').offsetWidth - margin.left - margin.right;
        const height = 500 - margin.top - margin.bottom;
        
        // Create SVG
        const svg = d3.select('#kinase-heatmap-container')
            .append('svg')
            .attr('width', width + margin.left + margin.right)
            .attr('height', height + margin.top + margin.bottom)
            .append('g')
            .attr('transform', `translate(${{margin.left}},${{margin.top}})`);
        
        // Create scales
        const xScale = d3.scaleBand()
            .domain(kinases)
            .range([0, width])
            .padding(0.05);
        
        const yScale = d3.scaleBand()
            .domain(sites)
            .range([0, height])
            .padding(0.05);
        
        // Find min and max values in the data
        let minValue = Infinity;
        let maxValue = -Infinity;
        for (let i = 0; i < scores.length; i++) {{
            for (let j = 0; j < scores[i].length; j++) {{
                if (scores[i][j] < minValue) minValue = scores[i][j];
                if (scores[i][j] > maxValue) maxValue = scores[i][j];
            }}
        }}
        
        // Create color scale
        const colorScale = d3.scaleSequential()
            .interpolator(d3.interpolateViridis)
            .domain([minValue, maxValue]);
        
        // Create tooltip
        const tooltip = d3.select("body")
            .append("div")
            .attr("class", "tooltip")
            .style("position", "absolute")
            .style("background-color", "white")
            .style("border", "1px solid #ddd")
            .style("border-radius", "4px")
            .style("padding", "10px")
            .style("z-index", "10")
            .style("visibility", "hidden");
        
        // Create cells
        svg.selectAll()
            .data(scores.flatMap((row, i) => row.map((value, j) => ({{
                site: sites[i],
                kinase: kinases[j],
                value: value,
                i: i,
                j: j
            }}))))
            .enter()
            .append("rect")
            .attr("x", d => xScale(d.kinase))
            .attr("y", d => yScale(d.site))
            .attr("width", xScale.bandwidth())
            .attr("height", yScale.bandwidth())
            .style("fill", d => colorScale(d.value))
            .on("mouseover", function(event, d) {{
                d3.select(this)
                    .style("stroke", "black")
                    .style("stroke-width", 2);
                
                tooltip
                    .style("visibility", "visible")
                    .html(`<strong>Site:</strong> ${{d.site}}<br><strong>Kinase:</strong> ${{d.kinase}}<br><strong>Score:</strong> ${{d.value.toFixed(3)}}`);
            }})
            .on("mousemove", function(event) {{
                tooltip
                    .style("top", (event.pageY - 10) + "px")
                    .style("left", (event.pageX + 10) + "px");
            }})
            .on("mouseout", function() {{
                d3.select(this)
                    .style("stroke", "none");
                
                tooltip
                    .style("visibility", "hidden");
            }});
        
        // Add x-axis
        svg.append("g")
            .style("font-size", "10px")
            .attr("transform", `translate(0,${{height}})`)
            .call(d3.axisBottom(xScale))
            .selectAll("text")
            .attr("transform", "rotate(-65)")
            .style("text-anchor", "end")
            .attr("dx", "-.8em")
            .attr("dy", ".15em");
        
        // Add y-axis
        svg.append("g")
            .style("font-size", "10px")
            .call(d3.axisLeft(yScale));
        
        // Add title
        svg.append("text")
            .attr("x", width / 2)
            .attr("y", -margin.top / 2)
            .attr("text-anchor", "middle")
            .style("font-size", "16px")
            .text("Kinase Prediction Scores");
        
        // Add legend
        const legendWidth = 200;
        const legendHeight = 20;
        
        const legend = svg.append("g")
            .attr("transform", `translate(${{width - legendWidth - 20}},-30)`);
        
        const defs = svg.append("defs");
        
        const linearGradient = defs.append("linearGradient")
            .attr("id", "linear-gradient");
        
        linearGradient.selectAll("stop")
            .data([
                {{offset: "0%", color: colorScale(minValue)}},
                {{offset: "100%", color: colorScale(maxValue)}}
            ])
            .enter().append("stop")
            .attr("offset", d => d.offset)
            .attr("stop-color", d => d.color);
        
        legend.append("rect")
            .attr("width", legendWidth)
            .attr("height", legendHeight)
            .style("fill", "url(#linear-gradient)");
        
        legend.append("text")
            .attr("x", 0)
            .attr("y", legendHeight + 15)
            .text(minValue.toFixed(2))
            .style("font-size", "10px");
        
        legend.append("text")
            .attr("x", legendWidth)
            .attr("y", legendHeight + 15)
            .attr("text-anchor", "end")
            .text(maxValue.toFixed(2))
            .style("font-size", "10px");
        
        legend.append("text")
            .attr("x", legendWidth / 2)
            .attr("y", legendHeight + 30)
            .attr("text-anchor", "middle")
            .text("Prediction Score")
            .style("font-size", "10px");
    }});
    </script>
    """
    
    return html

def process_matches_for_kinase_scoring(
    protein_uniprot_id: str,
    phosphosites: List[Dict],
    structural_matches: Dict[str, List[Dict]],
    sequence_matches: Dict[str, List[Dict]],
    match_type: str = "structural",  # "structural", "sequence", or "combined"
    max_matches_per_site: int = 20,  # Limit number of matches per site for efficiency
    min_structural_similarity: float = 0.2,  # min threshold for structural similarity (inverse of RMSD)
    min_sequence_similarity: float = 0.6,  # min threshold for sequence similarity
    rm_threshold: float = 10.0  # Maximum RMSD value to include
) -> Dict:
    """
    Process structural and/or sequence matches to prepare for kinase scoring analysis.
    This function maps match data to weight values and site IDs for kinase score aggregation.
    
    Args:
        protein_uniprot_id: UniProt ID of the source protein
        phosphosites: List of phosphosite dicts from the source protein
        structural_matches: Dict mapping site names to lists of structural match dicts
        sequence_matches: Dict mapping site names to lists of sequence match dicts
        match_type: Which type of matches to process - "structural", "sequence", or "combined"
        max_matches_per_site: Maximum number of matches to include per site
        min_structural_similarity: Minimum structural similarity threshold (inverse of RMSD)
        min_sequence_similarity: Minimum sequence similarity threshold
        rm_threshold: Maximum RMSD value to include
        
    Returns:
        Dict with processed match data ready for kinase analysis
    """
    logger.info(f"[DEBUG process_matches] Starting match processing for protein: {protein_uniprot_id}")
    logger.info(f"[DEBUG process_matches] Phosphosites: {len(phosphosites) if phosphosites else 0}, "
                f"Match type: {match_type}")
    logger.info(f"[DEBUG process_matches] Structural matches: {len(structural_matches) if structural_matches else 0}, "
                f"Sequence matches: {len(sequence_matches) if sequence_matches else 0}")
    
    # Debug info for structural_matches
    if structural_matches:
        logger.info(f"[DEBUG process_matches] Structural match keys: {list(structural_matches.keys())[:5]}")
        if structural_matches and list(structural_matches.keys()):
            first_key = list(structural_matches.keys())[0]
            first_matches = structural_matches[first_key]
            logger.info(f"[DEBUG process_matches] Sample structural match for {first_key}: {first_matches[0] if first_matches else 'No matches'}")
    
    # Debug info for sequence_matches
    if sequence_matches:
        logger.info(f"[DEBUG process_matches] Sequence match keys: {list(sequence_matches.keys())[:5]}")
        if sequence_matches and list(sequence_matches.keys()):
            first_key = list(sequence_matches.keys())[0]
            first_matches = sequence_matches[first_key]
            logger.info(f"[DEBUG process_matches] Sample sequence match for {first_key}: {first_matches[0] if first_matches else 'No matches'}")
        
    # Debug info for phosphosites
    if phosphosites:
        logger.info(f"[DEBUG process_matches] Sample phosphosite: {phosphosites[0]}")
        keys_in_phosphosites = list(phosphosites[0].keys())
        logger.info(f"[DEBUG process_matches] Keys in phosphosite: {keys_in_phosphosites}")

    result = {
        "protein_id": protein_uniprot_id,
        "sites": [],  # List of phosphosites with matches
        "site_matches": {},  # Dict mapping site IDs to lists of match site IDs
        "site_weights": {},  # Dict mapping site IDs to dicts of match weights
        "match_counts": {}  # Dict with match count statistics
    }
    
    # Count matches
    total_structural_matches = 0
    total_sequence_matches = 0
    
    # Process each phosphosite
    for site in phosphosites:
        # Extract site info - try multiple field names
        site_name = site.get('site') or site.get('name') or site.get('site_name')
        if not site_name:
            logger.warning(f"[DEBUG process_matches] Skipping site without name: {site}")
            continue
            
        # Skip tyrosine sites - they use different kinases
        if site_name and site_name[0] == 'Y':
            logger.info(f"[DEBUG process_matches] Skipping tyrosine site: {site_name}")
            continue
        
        # Get residue number - try multiple field names
        resno = site.get('resno') or site.get('residue_number')
        if resno is None:
            # Try to extract from site name
            if site_name and site_name[0] in 'STY' and site_name[1:].isdigit():
                resno = int(site_name[1:])
            else:
                logger.warning(f"[DEBUG process_matches] Skipping site without residue number: {site_name}")
                continue
            
        # Create canonical site ID
        site_id = f"{protein_uniprot_id}_{resno}"
        logger.info(f"[DEBUG process_matches] Processing site: {site_name}, ID: {site_id}")
        
        # Initialize match collections for this site
        site_match_ids = []
        site_match_weights = {}
        used_match_ids = set()  # Track matches we've already processed
        
        # Process structural matches if requested
        if match_type in ["structural", "combined"] and structural_matches and site_name in structural_matches:
            matches = structural_matches[site_name]
            logger.info(f"[DEBUG process_matches] Found {len(matches)} structural matches for {site_name}")
            
            # Convert RMSD to similarity score (1/RMSD) and sort
            weighted_matches = []
            for match in matches:
                if not isinstance(match, dict):
                    logger.warning(f"[DEBUG process_matches] Skipping non-dict match: {match}")
                    continue
                    
                # Skip if missing essential data
                if 'target_uniprot' not in match or 'target_site' not in match:
                    logger.warning(f"[DEBUG process_matches] Skipping match without target info: {match}")
                    continue
                    
                # Get RMSD with safety checks
                rmsd = match.get('rmsd')
                if rmsd is None or not isinstance(rmsd, (int, float)) or rmsd <= 0:
                    logger.warning(f"[DEBUG process_matches] Skipping match with invalid RMSD: {match}")
                    continue
                
                # Skip if above threshold
                if rmsd > rm_threshold:
                    continue
                
                # Skip tyrosine sites
                if match['target_site'] and match['target_site'][0] == 'Y':
                    continue
                
                # Convert RMSD to similarity (higher is better)
                similarity = min(1.0, 1.0 / rmsd)
                
                # Skip if below threshold
                if similarity < min_structural_similarity:
                    continue
                
                # Create target site ID - extract residue number
                target_site = match['target_site']
                target_resno = None
                
                # Try different methods to extract the residue number
                if target_site[0] in ['S', 'T', 'Y'] and target_site[1:].isdigit():
                    # Format like S123
                    target_resno = int(target_site[1:])
                else:
                    # Try to extract digits
                    digits = ''.join(c for c in target_site if c.isdigit())
                    if digits:
                        target_resno = int(digits)
                    # If still not found, check if target_id contains it
                    elif 'target_id' in match and '_' in match['target_id']:
                        # Format like P04637_123
                        parts = match['target_id'].split('_')
                        if len(parts) > 1 and parts[1].isdigit():
                            target_resno = int(parts[1])
                
                if target_resno is None:
                    logger.warning(f"[DEBUG process_matches] Couldn't extract residue number from: {target_site}")
                    continue
                
                target_id = f"{match['target_uniprot']}_{target_resno}"
                
                # Skip self-references
                if target_id == site_id:
                    continue
                    
                # Skip if we've already included this match
                if target_id in used_match_ids:
                    continue
                    
                # Add to weighted matches
                weighted_matches.append((target_id, similarity, match))
                used_match_ids.add(target_id)
            
            # Sort by similarity (highest first) and take top matches
            weighted_matches.sort(key=lambda x: x[1], reverse=True)
            for target_id, similarity, match in weighted_matches[:max_matches_per_site]:
                site_match_ids.append(target_id)
                site_match_weights[target_id] = similarity
            
            # Update match count
            total_structural_matches += len(site_match_ids)
            logger.info(f"[DEBUG process_matches] Added {len(site_match_ids)} structural matches for {site_name}")
        
        # Process sequence matches if requested
        if match_type in ["sequence", "combined"] and sequence_matches and site_name in sequence_matches:
            matches = sequence_matches[site_name]
            logger.info(f"[DEBUG process_matches] Found {len(matches)} sequence matches for {site_name}")
            
            # Filter and sort by similarity
            weighted_matches = []
            for match in matches:
                if not isinstance(match, dict):
                    continue
                    
                # Skip if missing essential data
                if 'target_uniprot' not in match or 'target_site' not in match:
                    continue
                    
                # Get similarity with safety checks
                similarity = match.get('similarity')
                if similarity is None or not isinstance(similarity, (int, float)) or similarity <= 0:
                    continue
                
                # Skip if below threshold
                if similarity < min_sequence_similarity:
                    continue
                
                # Skip tyrosine sites
                if match['target_site'] and match['target_site'][0] == 'Y':
                    continue
                
                # Create target site ID - extract residue number
                target_site = match['target_site']
                target_resno = None
                
                # Try different methods to extract the residue number
                if target_site[0] in ['S', 'T', 'Y'] and target_site[1:].isdigit():
                    # Format like S123
                    target_resno = int(target_site[1:])
                else:
                    # Try to extract digits
                    digits = ''.join(c for c in target_site if c.isdigit())
                    if digits:
                        target_resno = int(digits)
                    # If still not found, check if target_id contains it
                    elif 'target_id' in match and '_' in match['target_id']:
                        # Format like P04637_123
                        parts = match['target_id'].split('_')
                        if len(parts) > 1 and parts[1].isdigit():
                            target_resno = int(parts[1])
                
                if target_resno is None:
                    logger.warning(f"[DEBUG process_matches] Couldn't extract residue number from: {target_site}")
                    continue
                
                target_id = f"{match['target_uniprot']}_{target_resno}"
                
                # Skip self-references
                if target_id == site_id:
                    continue
                    
                # Skip if we've already included this match
                if target_id in used_match_ids:
                    continue
                    
                # Add to weighted matches
                weighted_matches.append((target_id, similarity, match))
                used_match_ids.add(target_id)
            
            # Sort by similarity (highest first) and take top matches
            weighted_matches.sort(key=lambda x: x[1], reverse=True)
            for target_id, similarity, match in weighted_matches[:max_matches_per_site]:
                site_match_ids.append(target_id)
                site_match_weights[target_id] = similarity
            
            # Update match count
            total_sequence_matches += len(site_match_ids)
            logger.info(f"[DEBUG process_matches] Added {len(site_match_ids)} sequence matches for {site_name}")
        
        # Only add sites that have matches
        if site_match_ids:
            # Add this site to the result
            result["sites"].append({
                "site_id": site_id,
                "site_name": site_name,
                "residue_type": site_name[0] if site_name else "S",
                "match_count": len(site_match_ids)
            })
            
            # Add matches and weights
            result["site_matches"][site_id] = site_match_ids
            result["site_weights"][site_id] = site_match_weights
            logger.info(f"[DEBUG process_matches] Site {site_name} has {len(site_match_ids)} total matches")
    
    # Add match count statistics
    result["match_counts"] = {
        "total_sites": len(result["sites"]),
        "total_structural_matches": total_structural_matches,
        "total_sequence_matches": total_sequence_matches,
        "unique_matches": len(set().union(*[set(matches) for matches in result["site_matches"].values()])) 
              if result["site_matches"] else 0
    }
    
    logger.info(f"[DEBUG process_matches] Processed match data: {result['match_counts']}")
    logger.info(f"[DEBUG process_matches] PROCESSED KINASE SCORING DATA")
    logger.info(f"{result}")
    
    return result

def create_protein_kinase_report(
    protein_uniprot_id: str,
    phosphosites: List[Dict],
    structural_matches: Dict[str, List[Dict]] = None,
    sequence_matches: Dict[str, List[Dict]] = None,
    match_type: str = "structural",
    include_heatmap: bool = True,
    include_family_analysis: bool = True,
    top_n_kinases: int = 10
) -> Dict:
    """
    Create a comprehensive kinase analysis report for a protein's phosphosites.
    
    Args:
        protein_uniprot_id: UniProt ID of the protein
        phosphosites: List of phosphosite dicts
        structural_matches: Dict mapping site names to lists of structural match dicts
        sequence_matches: Dict mapping site names to lists of sequence match dicts
        match_type: Which type of matches to process - "structural", "sequence", or "combined"
        include_heatmap: Whether to include heatmap visualization data
        include_family_analysis: Whether to include kinase family analysis
        top_n_kinases: Number of top kinases to include per site
        
    Returns:
        Dict with comprehensive kinase analysis data
    """
    # Process matches to prepare for kinase analysis
    processed_data = process_matches_for_kinase_scoring(
        protein_uniprot_id,
        phosphosites,
        structural_matches,
        sequence_matches,
        match_type=match_type
    )
    print("PROCESSED KINASE SCORING DATA")
    print(processed_data)
    # Initialize result
    result = {
        "protein_id": protein_uniprot_id,
        "sites": processed_data["sites"],
        "match_stats": processed_data["match_counts"],
        "site_predictions": {},
        "match_type": match_type
    }
    
    # Analyze each site with matches
    for site_info in processed_data["sites"]:
        site_id = site_info["site_id"]
        
        # Get match data for this site
        match_ids = processed_data["site_matches"].get(site_id, [])
        match_weights = processed_data["site_weights"].get(site_id, {})
        
        # Get top kinases for this site
        top_kinases = predict_top_kinases(
            site_id,
            match_ids,
            match_weights,
            top_n=top_n_kinases
        )
        
        # Store in result
        result["site_predictions"][site_id] = {
            "site_name": site_info["site_name"],
            "top_kinases": top_kinases,
            "match_count": len(match_ids)
        }
        
        # Add family analysis if requested
        if include_family_analysis:
            family_scores = get_kinase_family_distribution(
                site_id,
                match_ids,
                match_weights,
                top_n_families=8
            )
            
            result["site_predictions"][site_id]["family_scores"] = family_scores
    
    # Add heatmap data if requested
    if include_heatmap and processed_data["sites"]:
        # Get all site IDs
        all_site_ids = []
        for site_info in processed_data["sites"]:
            all_site_ids.append(site_info["site_id"])
            # Include a subset of matches for each site
            site_id = site_info["site_id"]
            match_ids = processed_data["site_matches"].get(site_id, [])
            match_weights = processed_data["site_weights"].get(site_id, {})
            
            # Sort matches by weight and take top 5
            top_matches = sorted(
                [(match_id, match_weights.get(match_id, 0.5)) for match_id in match_ids],
                key=lambda x: x[1],
                reverse=True
            )[:5]
            
            all_site_ids.extend([match_id for match_id, _ in top_matches])
        
        # Generate heatmap data
        heatmap_data = create_kinase_heatmap(
            all_site_ids,
            top_n_kinases=20,
            score_type="auto"
        )
        
        # Add to result
        result["heatmap_data"] = heatmap_data
    
    return result

def get_html_protein_kinase_report(report_data: Dict) -> str:
    """
    Generate an HTML representation of a protein kinase report.
    
    Args:
        report_data: Dict with kinase analysis data from create_protein_kinase_report
        
    Returns:
        HTML string with formatted report
    """
    if not report_data or "sites" not in report_data or not report_data["sites"]:
        return '<div class="alert alert-warning">No phosphosite data available for kinase analysis.</div>'
    
    # Extract data
    protein_id = report_data["protein_id"]
    sites = report_data["sites"]
    match_stats = report_data["match_stats"]
    site_predictions = report_data["site_predictions"]
    match_type = report_data["match_type"]
    
    # Create HTML sections
    sections = []
    
    # Summary section
    summary_html = f"""
    <div class="card mb-4">
        <div class="card-header">
            <h5 class="mb-0">Kinase Prediction Summary</h5>
        </div>
        <div class="card-body">
            <p>Analysis for protein <strong>{protein_id}</strong> using <strong>{match_type}</strong> similarity.</p>
            <div class="mb-3">
                <h6>Match Statistics:</h6>
                <ul>
                    <li><strong>Phosphosites analyzed:</strong> {match_stats["total_sites"]}</li>
                    <li><strong>Total structural matches:</strong> {match_stats["total_structural_matches"]}</li>
                    <li><strong>Total sequence matches:</strong> {match_stats["total_sequence_matches"]}</li>
                    <li><strong>Unique match sites:</strong> {match_stats["unique_matches"]}</li>
                </ul>
            </div>
            <p class="text-muted">
                This analysis predicts kinases for each phosphosite based on similar sites in other proteins,
                aggregating scores from the Cantley lab kinase prediction models.
            </p>
        </div>
    </div>
    """
    sections.append(summary_html)
    
    # Add heatmap if available
    if "heatmap_data" in report_data and "error" not in report_data["heatmap_data"]:
        heatmap_html = convert_to_html_heatmap(report_data["heatmap_data"])
        sections.append(heatmap_html)
    
    # Site predictions section
    for site_info in sites:
        site_id = site_info["site_id"]
        site_name = site_info["site_name"]
        
        # Skip if no predictions
        if site_id not in site_predictions:
            continue
            
        prediction = site_predictions[site_id]
        top_kinases = prediction["top_kinases"]
        match_count = prediction["match_count"]
        
        # Skip if no top kinases
        if not top_kinases:
            continue
        
        # Create kinase table
        kinase_rows = ""
        for i, kinase_info in enumerate(top_kinases):
            kinase = kinase_info["kinase"]
            score = kinase_info["score"]
            family = kinase_info["family"]
            
            # Format score as percentage
            score_percent = f"{score * 100:.1f}%"
            
            # Create progress bar with color based on rank
            color_class = "bg-success" if i < 3 else "bg-primary" if i < 5 else "bg-secondary"
            progress_width = f"{min(100, score * 100):.1f}%"
            
            # Create table row
            kinase_rows += f"""
            <tr>
                <td>{i+1}</td>
                <td>{kinase}</td>
                <td>{family}</td>
                <td>
                    <div class="d-flex align-items-center">
                        <div class="progress flex-grow-1" style="height: 12px;">
                            <div class="progress-bar {color_class}" role="progressbar" 
                                 style="width: {progress_width};" 
                                 aria-valuenow="{score * 100}" aria-valuemin="0" aria-valuemax="100">
                            </div>
                        </div>
                        <span class="ms-2">{score_percent}</span>
                    </div>
                </td>
            </tr>
            """
        
        # Create family visualization if available
        family_html = ""
        if "family_scores" in prediction and prediction["family_scores"]:
            family_scores = prediction["family_scores"]
            
            # Create family bars
            family_bars = ""
            for family, score in family_scores.items():
                # Format score and width
                score_percent = f"{score * 100:.1f}%"
                width_percent = f"{min(100, score * 100):.1f}%"
                
                family_bars += f"""
                <div class="mb-2">
                    <div class="d-flex justify-content-between mb-1">
                        <span>{family}</span>
                        <span>{score_percent}</span>
                    </div>
                    <div class="progress" style="height: 10px;">
                        <div class="progress-bar bg-info" role="progressbar" 
                             style="width: {width_percent};" 
                             aria-valuenow="{score * 100}" aria-valuemin="0" aria-valuemax="100">
                        </div>
                    </div>
                </div>
                """
            
            family_html = f"""
            <div class="col-md-5">
                <div class="card h-100">
                    <div class="card-header">
                        <h6 class="mb-0">Kinase Family Distribution</h6>
                    </div>
                    <div class="card-body">
                        {family_bars}
                    </div>
                </div>
            </div>
            """
        
        # Create site card
        site_html = f"""
        <div class="card mb-4">
            <div class="card-header bg-primary text-white">
                <h5 class="mb-0">Phosphosite {site_name} - Top Predicted Kinases</h5>
            </div>
            <div class="card-body">
                <p class="mb-3">
                    Based on analysis of <strong>{match_count}</strong> similar phosphosites, 
                    the following kinases are predicted to phosphorylate this site:
                </p>
                
                <div class="row">
                    <div class="col-md-{7 if family_html else 12}">
                        <div class="table-responsive">
                            <table class="table table-sm">
                                <thead>
                                    <tr>
                                        <th>#</th>
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
                    {family_html}
                </div>
                
                <div class="mt-3">
                    <a href="/site/{protein_id}/{site_name}" class="btn btn-outline-primary">
                        View Site Details
                    </a>
                </div>
            </div>
        </div>
        """
        sections.append(site_html)
    
    # Combine sections
    full_html = "\n".join(sections)
    
    return full_html

# ---- Utility functions for use in the web application ----

def get_site_prediction_data(site_id: str, match_site_ids: List[str] = None, 
                          similarity_weights: Dict[str, float] = None) -> Dict:
    """
    Get comprehensive kinase prediction data for a specific site.
    Suitable for use in web application routes.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        match_site_ids: Optional list of similar site IDs
        similarity_weights: Optional dict mapping site IDs to similarity weights
        
    Returns:
        Dict with prediction data, top kinases, and family analysis
    """
    # Get top kinases
    top_kinases = predict_top_kinases(
        site_id,
        match_site_ids,
        similarity_weights,
        top_n=20
    )
    
    # Get family distribution
    family_scores = get_kinase_family_distribution(
        site_id,
        match_site_ids,
        similarity_weights,
        top_n_families=10
    )
    
    # If matches were provided, create heatmap data
    heatmap_data = None
    if match_site_ids and len(match_site_ids) > 0:
        all_site_ids = [site_id] + match_site_ids[:20]  # Limit to top 20 matches
        heatmap_data = create_kinase_heatmap(
            all_site_ids,
            top_n_kinases=15,
            score_type="auto",
            similarity_weights=similarity_weights
        )
    
    # Determine residue type
    residue_type = get_residue_type(site_id)
    
    # Return comprehensive result
    return {
        "site_id": site_id,
        "residue_type": residue_type, 
        "top_kinases": top_kinases,
        "family_scores": family_scores,
        "heatmap_data": heatmap_data,
        "match_count": len(match_site_ids) if match_site_ids else 0
    }

def get_proteins_kinase_heatmap(protein_uniprot_id: str, phosphosites: List[Dict], 
                              structural_matches: Dict[str, List[Dict]] = None,
                              sequence_matches: Dict[str, List[Dict]] = None,
                              match_type: str = "structural") -> str:
    """
    Generate an HTML heatmap visualization for a protein's phosphosites.
    This is a convenience function for direct use in web application routes.
    
    Args:
        protein_uniprot_id: UniProt ID of the protein
        phosphosites: List of phosphosite dicts
        structural_matches: Dict mapping site names to lists of structural match dicts
        sequence_matches: Dict mapping site names to lists of sequence match dicts
        match_type: Which type of matches to process - "structural", "sequence", or "combined"
        
    Returns:
        HTML string with heatmap visualization
    """
    # Process matches
    processed_data = process_matches_for_kinase_scoring(
        protein_uniprot_id,
        phosphosites,
        structural_matches,
        sequence_matches,
        match_type=match_type
    )
    
    if not processed_data["sites"]:
        return '<div class="alert alert-warning">No phosphosite data available for kinase analysis.</div>'
    
    # Get all site IDs
    all_site_ids = []
    for site_info in processed_data["sites"]:
        all_site_ids.append(site_info["site_id"])
        # Include a subset of matches for each site
        site_id = site_info["site_id"]
        match_ids = processed_data["site_matches"].get(site_id, [])
        match_weights = processed_data["site_weights"].get(site_id, {})
        
        # Sort matches by weight and take top 5
        top_matches = sorted(
            [(match_id, match_weights.get(match_id, 0.5)) for match_id in match_ids],
            key=lambda x: x[1],
            reverse=True
        )[:5]
        
        all_site_ids.extend([match_id for match_id, _ in top_matches])
    
    # Generate heatmap data
    heatmap_data = create_kinase_heatmap(
        all_site_ids,
        top_n_kinases=20,
        score_type="auto"
    )
    
    # Convert to HTML
    return convert_to_html_heatmap(heatmap_data)