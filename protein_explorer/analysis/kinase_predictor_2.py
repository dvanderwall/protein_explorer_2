"""
Functions for kinase prediction based on structural and sequence similarity.
Uses database queries instead of file loading.
"""

import os
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Union
import logging

# Import database module
from protein_explorer.db import get_kinase_scores, get_kinase_scores_batch

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_kinase_scores(file_path: str = None, score_type: str = 'structure') -> bool:
    """
    Load kinase scores from the database.
    
    Args:
        file_path: Parameter kept for backward compatibility but not used
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Boolean indicating if scores can be accessed
    """
    # This function is kept for backward compatibility but now uses the database
    logger.info(f"Using database for {score_type} kinase scores")
    return True

def get_site_kinase_scores(site_id: str, score_type: str = 'structure') -> Dict:
    """
    Get kinase scores for a specific site.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary with kinase names as keys and scores as values
    """
    scores_data = get_kinase_scores(site_id, score_type)
    
    if not scores_data:
        logger.warning(f"No {score_type} kinase scores available for {site_id}")
        return {}
    
    return scores_data

def predict_kinases(site_id: str, top_n: int = 5, score_type: str = 'structure') -> List[Dict]:
    """
    Get top N predicted kinases for a site.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        top_n: Number of top kinases to return
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        List of dictionaries with kinase names and scores
    """
    # Get all kinase scores for the site
    site_data = get_site_kinase_scores(site_id, score_type)
    
    if not site_data or 'scores' not in site_data:
        logger.warning(f"No {score_type} kinase scores available for {site_id}")
        return []
    
    # Get all kinase scores
    scores = site_data['scores']
    
    # Sort by score (descending)
    sorted_scores = sorted(scores.items(), key=lambda x: float(x[1]), reverse=True)
    
    # Return top N kinases
    top_kinases = []
    for kinase, score in sorted_scores[:top_n]:
        top_kinases.append({
            'kinase': kinase,
            'score': float(score)  # Ensure score is a float
        })
    
    return top_kinases

def compare_kinase_scores(site_ids: List[str], top_n: int = 5, score_type: str = 'structure') -> Dict:
    """
    Compare kinase scores across multiple sites.
    
    Args:
        site_ids: List of site IDs
        top_n: Number of top kinases to consider for each site
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary with comparison data
    """
    if not site_ids:
        return {}
    
    # Get top kinases for each site
    site_kinases = {}
    all_kinases = set()
    
    for site_id in site_ids:
        top_kinases = predict_kinases(site_id, top_n, score_type)
        site_kinases[site_id] = top_kinases
        
        # Collect all unique kinases
        all_kinases.update([k['kinase'] for k in top_kinases])
    
    # Prepare comparison data
    comparison = {
        'sites': site_ids,
        'kinases': list(all_kinases),
        'data': {}
    }
    
    # Get all site scores in a batch query
    batch_scores = get_kinase_scores_batch(site_ids, score_type)
    
    # Add scores for each site and kinase
    for site_id in site_ids:
        comparison['data'][site_id] = {}
        site_data = batch_scores.get(site_id, {})
        
        if site_data and 'scores' in site_data:
            scores = site_data['scores']
            for kinase in all_kinases:
                comparison['data'][site_id][kinase] = scores.get(kinase, 0)
    
    return comparison

def get_heatmap_data(site_ids: List[str], top_n: int = 10, score_type: str = 'structure') -> Dict:
    """
    Get data for heatmap visualization of kinase scores.
    
    Args:
        site_ids: List of site IDs
        top_n: Number of top kinases to include
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary with heatmap data
    """
    if not site_ids:
        return {}
    
    # Get all site data in a single batch query
    batch_scores = get_kinase_scores_batch(site_ids, score_type)
    
    if not batch_scores:
        logger.warning(f"No {score_type} kinase scores available")
        return {}
    
    try:
        # Get data for the specified sites
        sites_data = []
        for site_id in site_ids:
            if site_id in batch_scores and 'scores' in batch_scores[site_id]:
                # Add this site's scores to our dataset
                site_scores = batch_scores[site_id]['scores']
                sites_data.append((site_id, site_scores))
        
        if not sites_data:
            logger.warning(f"None of the specified sites found in {score_type} kinase scores")
            return {}
        
        # Convert to a combined dictionary
        combined_scores = {}
        for site_id, scores in sites_data:
            for kinase, score in scores.items():
                if kinase not in combined_scores:
                    combined_scores[kinase] = []
                combined_scores[kinase].append(score)
        
        # Calculate mean score for each kinase across all sites
        mean_scores = {}
        for kinase, scores in combined_scores.items():
            mean_scores[kinase] = sum(scores) / len(scores)
        
        # Get top N kinases by mean score
        top_kinases = sorted(mean_scores.items(), key=lambda x: x[1], reverse=True)[:top_n]
        top_kinase_names = [k[0] for k in top_kinases]
        
        # Prepare heatmap data
        heatmap_data = {
            'sites': site_ids,
            'kinases': top_kinase_names,
            'scores': []
        }
        
        # Add scores for each site and kinase
        for site_id in site_ids:
            if site_id in batch_scores and 'scores' in batch_scores[site_id]:
                site_scores = batch_scores[site_id]['scores']
                for kinase in top_kinase_names:
                    score = site_scores.get(kinase, 0)
                    heatmap_data['scores'].append({
                        'site': site_id,
                        'kinase': kinase,
                        'score': float(score)
                    })
        
        return heatmap_data
    except Exception as e:
        logger.error(f"Error generating heatmap data: {e}")
        return {}

def get_kinase_radar_data(site_id: str, top_n: int = 5, score_type: str = 'structure') -> Dict:
    """
    Get data for radar chart visualization of kinase scores.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        top_n: Number of top kinases to include
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary with radar chart data
    """
    # Get top kinases for the site
    top_kinases = predict_kinases(site_id, top_n, score_type)
    
    if not top_kinases:
        return {}
    
    # Prepare radar chart data
    radar_data = {
        'labels': [k['kinase'] for k in top_kinases],
        'datasets': [{
            'label': f"{score_type.capitalize()} Kinase Scores",
            'data': [k['score'] for k in top_kinases]
        }]
    }
    
    return radar_data

def get_kinase_comparison_data(site_id: str, score_types: List[str] = ['structure', 'sequence'], top_n: int = 5) -> Dict:
    """
    Get comparison data between structure and sequence kinase scores.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        score_types: List of score types to compare
        top_n: Number of top kinases to include
        
    Returns:
        Dictionary with comparison data
    """
    # Get all unique kinases from both score types
    all_kinases = set()
    
    for score_type in score_types:
        top_kinases = predict_kinases(site_id, top_n, score_type)
        all_kinases.update([k['kinase'] for k in top_kinases])
    
    # Prepare comparison data
    comparison_data = {
        'kinases': list(all_kinases),
        'datasets': []
    }
    
    # Add scores for each score type
    for score_type in score_types:
        site_data = get_site_kinase_scores(site_id, score_type)
        
        if site_data and 'scores' in site_data:
            scores = site_data['scores']
            
            dataset = {
                'label': f"{score_type.capitalize()} Score",
                'data': [scores.get(kinase, 0) for kinase in all_kinases]
            }
            
            comparison_data['datasets'].append(dataset)
    
    return comparison_data

def get_known_kinase_info(site_id: str, score_type: str = 'structure') -> Dict:
    """
    Get information about the known kinase for a site, if available.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        score_type: Type of scores to check
        
    Returns:
        Dictionary with known kinase information
    """
    # Get site data
    site_data = get_site_kinase_scores(site_id, score_type)
    
    if not site_data:
        return {'has_known_kinase': False}
    
    known_kinase = site_data.get('known_kinase', 'unlabeled')
    
    if known_kinase == 'unlabeled':
        return {'has_known_kinase': False}
    
    return {
        'has_known_kinase': True,
        'kinase': known_kinase
    }

def categorize_kinases_by_family(kinases: List[Dict]) -> Dict:
    """
    Categorize kinases by family.
    
    Args:
        kinases: List of dictionaries with kinase names and scores
        
    Returns:
        Dictionary with kinase families and their scores
    """
    # Define kinase families
    kinase_families = {
        'CDK': ['CDK1', 'CDK2', 'CDK4', 'CDK5', 'CDK6', 'CDK7', 'CDK8', 'CDK9'],
        'MAPK': ['ERK1', 'ERK2', 'p38', 'JNK1', 'JNK2', 'JNK3'],
        'GSK': ['GSK3', 'GSK3A', 'GSK3B'],
        'CK': ['CK1', 'CK2', 'CSNK1', 'CSNK2'],
        'PKC': ['PKC', 'PKCALPHA', 'PKCBETA', 'PKCDELTA', 'PKCEPSILON', 'PKCGAMMA', 'PKCZETA'],
        'PKA': ['PKA', 'PKACA', 'PKACB', 'PKACG'],
        'AKT': ['AKT', 'AKT1', 'AKT2', 'AKT3'],
        'SRC': ['SRC', 'FYN', 'LCK', 'LYN', 'HCK', 'FGR', 'BLK', 'YES'],
        'CAMK': ['CAMK', 'CAMK1', 'CAMK2', 'CAMK4'],
        'ATM/ATR': ['ATM', 'ATR', 'DNAPK'],
        'PLK': ['PLK1', 'PLK2', 'PLK3', 'PLK4'],
        'AURORA': ['AURKA', 'AURKB', 'AURKC'],
        'Other': []
    }
    
    # Categorize kinases
    family_scores = {}
    
    for kinase_data in kinases:
        kinase_name = kinase_data['kinase']
        kinase_score = kinase_data['score']
        
        # Find the family for this kinase
        assigned = False
        for family, members in kinase_families.items():
            if any(member in kinase_name.upper() for member in members):
                if family not in family_scores:
                    family_scores[family] = 0
                family_scores[family] += kinase_score
                assigned = True
                break
        
        # If not assigned to any family, put in 'Other'
        if not assigned:
            if 'Other' not in family_scores:
                family_scores['Other'] = 0
            family_scores['Other'] += kinase_score
    
    # Return sorted by score (descending)
    return dict(sorted(family_scores.items(), key=lambda x: x[1], reverse=True))