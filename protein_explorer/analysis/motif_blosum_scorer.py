"""
Module for calculating BLOSUM62 similarity scores between phosphosite motifs.

This module provides functions to calculate similarity scores between 
phosphosite motifs using BLOSUM62 matrix, particularly focusing on 
the -4 to +4 region around the phosphorylation site.
"""

import numpy as np
import pandas as pd
import logging
import time
from functools import lru_cache
from typing import Dict, List, Tuple, Set, Optional
from protein_explorer.db.db import get_phosphosites_batch, get_phosphosite_data, execute_query, init_db

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# BLOSUM62 scoring matrix
BLOSUM62 = {
    'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0, 'B': -2, 'Z': -1, 'X': -1, '*': -4},
    'R': {'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3, 'B': -1, 'Z': 0, 'X': -1, '*': -4},
    'N': {'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3, 'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3, 'B': 3, 'Z': 0, 'X': -1, '*': -4},
    'D': {'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'Z': 1, 'X': -1, '*': -4},
    'C': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1, 'B': -3, 'Z': -3, 'X': -1, '*': -4},
    'Q': {'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2, 'B': 0, 'Z': 3, 'X': -1, '*': -4},
    'E': {'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'Z': 4, 'X': -1, '*': -4},
    'G': {'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3, 'B': -1, 'Z': -2, 'X': -1, '*': -4},
    'H': {'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3, 'B': 0, 'Z': 0, 'X': -1, '*': -4},
    'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4, 'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3, 'B': -3, 'Z': -3, 'X': -1, '*': -4},
    'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2, 'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1, 'B': -4, 'Z': -3, 'X': -1, '*': -4},
    'K': {'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'Z': 1, 'X': -1, '*': -4},
    'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1, 'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1, 'B': -3, 'Z': -1, 'X': -1, '*': -4},
    'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0, 'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1, 'B': -3, 'Z': -3, 'X': -1, '*': -4},
    'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2, 'B': -2, 'Z': -1, 'X': -1, '*': -4},
    'S': {'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2, 'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'Z': 0, 'X': -1, '*': -4},
    'T': {'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0, 'B': -1, 'Z': -1, 'X': -1, '*': -4},
    'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3, 'B': -4, 'Z': -3, 'X': -1, '*': -4},
    'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1, 'B': -3, 'Z': -2, 'X': -1, '*': -4},
    'V': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3, 'L': 1, 'K': -2, 'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4, 'B': -3, 'Z': -2, 'X': -1, '*': -4},
    'B': {'A': -2, 'R': -1, 'N': 3, 'D': 4, 'C': -3, 'Q': 0, 'E': 1, 'G': -1, 'H': 0, 'I': -3, 'L': -4, 'K': 0, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'Z': 1, 'X': -1, '*': -4},
    'Z': {'A': -1, 'R': 0, 'N': 0, 'D': 1, 'C': -3, 'Q': 3, 'E': 4, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'Z': 4, 'X': -1, '*': -4},
    'X': {'A': -1, 'R': -1, 'N': -1, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -1, 'H': -1, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -1, 'P': -1, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': -1, 'B': -1, 'Z': -1, 'X': -1, '*': -4},
    '*': {'A': -4, 'R': -4, 'N': -4, 'D': -4, 'C': -4, 'Q': -4, 'E': -4, 'G': -4, 'H': -4, 'I': -4, 'L': -4, 'K': -4, 'M': -4, 'F': -4, 'P': -4, 'S': -4, 'T': -4, 'W': -4, 'Y': -4, 'V': -4, 'B': -4, 'Z': -4, 'X': -4, '*': 1}
}


def extract_central_motif(motif: str, window_size: int = 4) -> str:
    """
    Extract the central region of a motif (-window_size to +window_size).
    
    Args:
        motif: The full motif sequence
        window_size: Size of the window on each side (default: 4)
        
    Returns:
        Central region of the motif
    """
    if not motif:
        return ""
    
    # Find center position (phosphosite position)
    center_pos = len(motif) // 2
    
    # Extract window around center
    start_pos = max(0, center_pos - window_size)
    end_pos = min(len(motif), center_pos + window_size + 1)
    central_motif = motif[start_pos:end_pos]
    
    return central_motif

def standardize_motif(motif: str, window_size: int = 4) -> str:
    """
    Standardize a motif to have exactly window_size positions before and after the center.
    Handles varying length motifs and pads with 'X' where needed.
    
    Args:
        motif: The motif sequence
        window_size: Size of the window on each side (default: 4)
        
    Returns:
        Standardized motif with exactly 2*window_size+1 characters
    """
    if not motif:
        return "X" * (2 * window_size + 1)
    
    # Find center position (phosphosite)
    center_pos = len(motif) // 2
    
    # Get the phosphosite and parts before/after
    site_char = motif[center_pos]
    before_site = motif[:center_pos]
    after_site = motif[center_pos + 1:]
    
    # Ensure exactly window_size positions before
    if len(before_site) < window_size:
        # Pad with X if shorter
        before_site = "X" * (window_size - len(before_site)) + before_site
    else:
        # Truncate if longer
        before_site = before_site[-window_size:]
    
    # Ensure exactly window_size positions after
    if len(after_site) < window_size:
        # Pad with X if shorter
        after_site = after_site + "X" * (window_size - len(after_site))
    else:
        # Truncate if longer
        after_site = after_site[:window_size]
    
    # Combine standardized parts
    return before_site + site_char + after_site

@lru_cache(maxsize=100000)
def calculate_blosum62_score_cached(motif1: str, motif2: str) -> float:
    """
    Calculate BLOSUM62 similarity score between two motifs with caching.
    
    Args:
        motif1: First motif sequence
        motif2: Second motif sequence
        
    Returns:
        BLOSUM62 similarity score
    """
    if not motif1 or not motif2 or len(motif1) != len(motif2):
        return 0.0
    
    # Calculate raw score
    raw_score = 0
    for i in range(len(motif1)):
        aa1 = motif1[i].upper()
        aa2 = motif2[i].upper()
        
        # Replace any unusual characters with X
        if aa1 not in BLOSUM62:
            aa1 = 'X'
        if aa2 not in BLOSUM62:
            aa2 = 'X'
        
        raw_score += BLOSUM62[aa1][aa2]
    
    # Normalize by motif length
    return raw_score / len(motif1)

def normalize_similarity_score(score: float) -> float:
    """
    Normalize BLOSUM62 score to a 0-1 range for compatibility with your similarity network.
    This is a heuristic normalization based on typical BLOSUM62 scores.
    
    Args:
        score: Raw BLOSUM62 score
        
    Returns:
        Normalized similarity score (0-1)
    """
    # Typical range of BLOSUM62 scores for -4:+4 motifs is roughly -20 to +20
    # We'll map this range to 0-1 for compatibility with your existing similarity scores
    min_score = -20.0
    max_score = 20.0
    
    # Normalize to 0-1 range
    normalized = (score - min_score) / (max_score - min_score)
    
    # Clamp to 0-1 range
    return max(0.0, min(1.0, normalized))

def find_blosum_matches_batch_optimized(unmatched_sites: List[Dict], min_similarity: float = 0.4) -> Dict[str, List[Dict]]:
    """
    Optimized version: Find BLOSUM matches for multiple unmatched sites in batch.
    This version uses trigram-based prefiltering and efficient database queries.
    
    Args:
        unmatched_sites: List of site dictionaries that don't have sequence matches
        min_similarity: Minimum normalized similarity score (0-1)
        
    Returns:
        Dictionary mapping site names to lists of match dictionaries
    """
    if not unmatched_sites:
        return {}
    
    logger.info(f"Finding optimized BLOSUM matches for {len(unmatched_sites)} unmatched sites")
    
    # Initialize results dictionary
    blosum_matches = {}
    
    # Keep track of timing for performance monitoring
    total_start_time = time.time()
    
    try:
        # Process in smaller batches for better memory management
        batch_size = 100
        for i in range(0, len(unmatched_sites), batch_size):
            batch = unmatched_sites[i:i+batch_size]
            logger.info(f"Processing batch {i//batch_size + 1}/{(len(unmatched_sites) + batch_size - 1)//batch_size} with {len(batch)} sites")
            batch_start_time = time.time()
            
            # Process each site in this batch
            for site in batch:
                
                if 'motif' not in site or not site['motif'] or '_' in site['motif']:
                    continue
                    
                site_name = site.get('site', '')
                if not site_name:
                    continue
                
                if site['is_known'] or (not site['is_known'] and site['surface_accessibility'] > 50):

                    site_type = site_name[0] if site_name and site_name[0] in 'STY' else 'S'
                    site_number = site_name[1:] if site_name else ''
                    
                    # Extract central motif (-4:+4) from the full motif
                    # Only do this for the query motif since we don't have central_motif stored for it yet
                    query_motif = site['motif']
                    central_query_motif = extract_central_motif(query_motif, window_size=4)
                    std_query_motif = standardize_motif(central_query_motif)
                    

                    # Skip if not exactly 9 characters (our expected motif length)
                    if len(std_query_motif) != 9:
                        logger.warning(f"Skipping site {site_name} with invalid motif length: {len(std_query_motif)}")
                        continue
                    
                    # Extract trigrams for pre-filtering
                    start_trigram = std_query_motif[0:3]
                    middle_trigram = std_query_motif[3:6]
                    end_trigram = std_query_motif[6:9]
                    
                    # Pre-filter potential matches using trigrams for speed
                    # This query directly uses the central_motif column we added
                    potential_matches_query = """
                        SELECT PhosphositeID, uniprot_id, Residue_Number, Residue, central_motif
                        FROM Phosphosite_Supplementary_Data
                        WHERE central_motif IS NOT NULL
                        AND LENGTH(central_motif) = 9
                        AND (
                            LEFT(central_motif, 3) = :start_trigram OR
                            SUBSTRING(central_motif, 4, 3) = :middle_trigram OR
                            RIGHT(central_motif, 3) = :end_trigram
                        )
                        AND `SITE_+/-7_AA` IS NOT NULL
                        AND `SITE_+/-7_AA` NOT LIKE '%\\_%'
                        AND `SITE_+/-7_AA` NOT LIKE '%?%'
                        LIMIT 200000
                    """
                    
                    # Execute query with trigram parameters
                    df = execute_query(potential_matches_query, {
                        "start_trigram": start_trigram,
                        "middle_trigram": middle_trigram,
                        "end_trigram": end_trigram
                    })
                    
                    if df.empty:
                        logger.info(f"No potential matches found for site {site_name}")
                        continue
                    
                    logger.info(f"Found {len(df)} potential matches for site {site_name}")
                    
                    # Process potential matches
                    site_matches = []
                    # Process each potential match
                    for _, row in df.iterrows():
                        # Skip self-matches
                        if (site.get('uniprot', '') == row['uniprot_id'] and 
                            str(site.get('resno', '')) == str(row['Residue_Number'])):
                            continue
                        
                        # Use the central_motif directly from the database - no need to extract again
                        match_motif = row['central_motif'].upper()
                        
                        # Calculate BLOSUM score using cached function
                        blosum_score = calculate_blosum62_score_cached(std_query_motif, match_motif)
                        
                        # Normalize to 0-1 range
                        similarity = normalize_similarity_score(blosum_score)
                        
                        # Filter by minimum similarity
                        if similarity >= min_similarity:
                            # Determine site type from residue column
                            residue_type = row['Residue'] if row['Residue'] in 'STY' else 'S'
                            # Build match dictionary
                            match_dict = {
                                'target_id': row['PhosphositeID'], 'target_uniprot': row['uniprot_id'],
                                'target_site': f"{residue_type}{row['Residue_Number']}", 'site_type': residue_type,
                                'similarity': similarity, 'blosum_score': blosum_score, 'match_type': 'blosum'
                            }
                            site_matches.append(match_dict)
                    
                    # Sort by similarity (highest first)
                    site_matches.sort(key=lambda x: x['similarity'], reverse=True)
                    
                    # Limit to top 50 matches
                    site_matches = site_matches[:20]
                    
                    # Add to results if we found matches
                    if site_matches:
                        blosum_matches[site_name] = site_matches
                        logger.info(f"Found {len(site_matches)} BLOSUM matches for site {site_name}")
                    else:
                        logger.info(f"No matching sites found for {site_name} after filtering")
                else:
                    logger.info(f"Skipping site {site_name} due to known status or low surface accessibility")
            
            batch_end_time = time.time()
            logger.info(f"Batch processed in {batch_end_time - batch_start_time:.2f} seconds")
        
        total_end_time = time.time()
        logger.info(f"All BLOSUM matches found in {total_end_time - total_start_time:.2f} seconds")
        return blosum_matches
    
    except Exception as e:
        logger.error(f"Error in batch BLOSUM matching: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return {}

def find_blosum_matches(query_motif: str, min_similarity: float = 0.4) -> List[Dict]:
    """
    Find matches based on BLOSUM62 scores with existing motifs in the database.
    This is a legacy function kept for backward compatibility.
    
    Args:
        query_motif: Motif sequence to find matches for
        min_similarity: Minimum normalized similarity score (0-1)
        
    Returns:
        List of matching site dictionaries with similarity scores
    """
    if not query_motif:
        return []
    
    logger.info(f"Finding BLOSUM matches for motif: {query_motif}")
    
    # Create a dummy site dictionary
    dummy_site = {
        'site': 'S0',
        'motif': query_motif
    }
    
    # Use the optimized batch function with a single site
    results = find_blosum_matches_batch_optimized([dummy_site], min_similarity)
    
    # Return the matches for the dummy site if found
    if 'S0' in results and results['S0']:
        return results['S0']
    
    return []