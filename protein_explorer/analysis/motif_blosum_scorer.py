"""
Module for calculating BLOSUM62 similarity scores between phosphosite motifs.

This module provides functions to calculate similarity scores between 
phosphosite motifs using BLOSUM62 matrix, particularly focusing on 
the -4 to +4 region around the phosphorylation site.
"""

import numpy as np
import pandas as pd
import logging
from typing import Dict, List, Tuple, Set, Optional
from protein_explorer.db.db import get_phosphosites_batch, get_phosphosite_data

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

def calculate_blosum62_score(motif1: str, motif2: str) -> float:
    """
    Calculate BLOSUM62 similarity score between two motifs.
    
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
    
    # Normalize by motif length to get scores in a more consistent range
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
    # Adjust these bounds if you find different score ranges in your actual data
    min_score = -20.0
    max_score = 20.0
    
    # Normalize to 0-1 range
    normalized = (score - min_score) / (max_score - min_score)
    
    # Clamp to 0-1 range
    return max(0.0, min(1.0, normalized))

def find_blosum_matches(query_motif: str, min_similarity: float = 0.4) -> List[Dict]:
    """
    Find matches based on BLOSUM62 scores with existing motifs in the database.
    
    Args:
        query_motif: Motif sequence to find matches for
        min_similarity: Minimum normalized similarity score (0-1)
        
    Returns:
        List of matching site dictionaries with similarity scores
    """
    if not query_motif:
        return []
    
    logger.info(f"Finding BLOSUM matches for motif: {query_motif}")
    
    # Extract central -4:+4 region
    central_query_motif = extract_central_motif(query_motif, window_size=4)
    
    # Standardize query motif
    std_query_motif = standardize_motif(central_query_motif)
    
    try:
        # Fetch a batch of phosphosites from the database to compare against
        # Limiting to a reasonable number to avoid excessive processing
        # Need to implement a smart batching strategy here based on your database structure
        from protein_explorer.db.db import execute_query
        
        # Query to get sites with valid motifs
        query = """
            SELECT PhosphositeID, uniprot_id, Residue_Number, `SITE_+/-7_AA` as motif
            FROM Phosphosite_Supplementary_Data
            WHERE `SITE_+/-7_AA` IS NOT NULL 
            AND `SITE_+/-7_AA` NOT LIKE '%\\_%'
            AND `SITE_+/-7_AA` NOT LIKE '%_%'
            LIMIT 200000
        """
        
        df = execute_query(query)
        
        if df.empty:
            logger.warning("No phosphosites found in database for BLOSUM comparison")
            return []
        
        logger.info(f"Found {len(df)} phosphosites in database for BLOSUM comparison")
        
        # Calculate BLOSUM scores for each potential match
        matches = []
        
        for _, row in df.iterrows():
            # Extract data
            site_id = row['PhosphositeID']
            uniprot_id = row['uniprot_id']
            residue_number = row['Residue_Number']
            motif = row['motif']
            
            # Skip if missing essential data
            if not site_id or not uniprot_id or not residue_number or not motif:
                continue
            
            # Extract central -4:+4 region
            central_motif = extract_central_motif(motif, window_size=4)
            
            # Standardize comparison motif
            std_motif = standardize_motif(central_motif)
            
            # Calculate BLOSUM62 score
            blosum_score = calculate_blosum62_score(std_query_motif, std_motif)
            
            # Normalize to 0-1 range for compatibility with your similarity network
            similarity = normalize_similarity_score(blosum_score)
            
            # Filter by minimum similarity
            if similarity >= min_similarity:
                # Determine site type from site_id or motif
                site_type = 'S'  # Default to S if we can't determine
                
                # Try to extract from site_id
                if '_' in site_id:
                    residue_part = site_id.split('_')[1]
                    if residue_part and len(residue_part) > 0:
                        if residue_part[0] in 'STY':
                            site_type = residue_part[0]
                
                # Build match dictionary with similar structure to your existing matches
                match = {
                    'target_id': site_id,
                    'target_uniprot': uniprot_id,
                    'target_site': f"{site_type}{residue_number}",
                    'site_type': site_type,
                    'similarity': similarity,
                    'blosum_score': blosum_score,  # Include raw score for reference
                    'motif': motif,
                    'match_type': 'blosum'  # Flag to indicate this is a BLOSUM-based match
                }
                
                matches.append(match)
        
        # Sort by similarity (highest first)
        matches.sort(key=lambda x: x['similarity'], reverse=True)
        
        # Limit the number of matches to a reasonable number (e.g., top 50)
        # Adjust this limit based on your application needs
        matches = matches[:50]
        
        logger.info(f"Found {len(matches)} BLOSUM matches with similarity >= {min_similarity}")
        
        return matches
    
    except Exception as e:
        logger.error(f"Error finding BLOSUM matches: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return []

def find_blosum_matches_batch(unmatched_sites: List[Dict], min_similarity: float = 0.4) -> Dict[str, List[Dict]]:
    """
    Find BLOSUM matches for multiple unmatched sites in batch.
    
    Args:
        unmatched_sites: List of site dictionaries that don't have sequence matches
        min_similarity: Minimum normalized similarity score (0-1)
        
    Returns:
        Dictionary mapping site names to lists of match dictionaries
    """
    if not unmatched_sites:
        return {}
    
    logger.info(f"Finding BLOSUM matches for {len(unmatched_sites)} unmatched sites")
    
    # Initialize results dictionary
    blosum_matches = {}
    
    try:
        # Fetch a batch of phosphosites from the database to compare against
        from protein_explorer.db.db import execute_query
        
        # Simplified query - SQL syntax might be causing issues
        # Adjusted to avoid potential SQL errors with escape sequences and wildcards
        query = """
            SELECT PhosphositeID, uniprot_id, Residue_Number, Residue, `SITE_+/-7_AA` as motif
            FROM Phosphosite_Supplementary_Data
            WHERE `SITE_+/-7_AA` IS NOT NULL 
            LIMIT 200000
        """
        
        df = execute_query(query)
        
        if df.empty:
            logger.warning("No phosphosites found in database for BLOSUM comparison")
            return {}
        
        logger.info(f"Found {len(df)} phosphosites in database for BLOSUM comparison")
        
        # Prepare database motifs for comparison - manually filter out ones with underscores
        db_motifs = []
        for _, row in df.iterrows():
            site_id = row['PhosphositeID']
            uniprot_id = row['uniprot_id']
            residue_number = row['Residue_Number']
            motif = row['motif'] if 'motif' in row else None
            
            # Get residue type, defaulting to S if not available
            residue_type = 'S'
            if 'Residue' in row and row['Residue'] and row['Residue'] in 'STY':
                residue_type = row['Residue']
            
            # Skip if missing essential data or if motif contains underscore
            if not site_id or not residue_number or not motif:
                continue
                
            # Skip motifs with underscores
            if '_' in motif:
                continue

            if not uniprot_id:
                if '_' in site_id:
                    uniprot_id = site_id.split('_')[0]
            
            # First extract just the central -4:+4 region (for +/-7 motifs)
            central_motif = extract_central_motif(motif, window_size=4)
            
            # Then standardize for comparison
            std_motif = standardize_motif(central_motif)
            
            db_motifs.append({
                'site_id': site_id,
                'uniprot_id': uniprot_id,
                'residue_number': residue_number,
                'residue_type': residue_type,
                'motif': motif,
                'std_motif': std_motif
            })
        
        logger.info(f"Prepared {len(db_motifs)} database motifs for comparison")
        
        # Process each unmatched site
        for site in unmatched_sites:
            # Skip sites without motifs
            if 'motif' not in site or not site['motif']:
                continue
            
            # Get site info
            site_name = site.get('site', '')
            if not site_name:
                continue
            
            site_type = site_name[0] if site_name and site_name[0] in 'STY' else 'S'
            site_number = site_name[1:] if site_name else ''
            
            # Extract central -4:+4 region of query motif
            query_motif = site['motif']
            central_query_motif = extract_central_motif(query_motif, window_size=4)
            
            # Standardize for comparison
            std_query_motif = standardize_motif(central_query_motif)
            
            # Calculate BLOSUM scores and find matches
            site_matches = []
            
            for db_site in db_motifs:
                # Skip self-matches (same uniprot and residue number)
                if (site.get('uniprot', '') == db_site['uniprot_id'] and 
                    str(site.get('resno', '')) == str(db_site['residue_number'])):
                    continue
                
                # Calculate BLOSUM62 score
                blosum_score = calculate_blosum62_score(std_query_motif, db_site['std_motif'])
                
                # Normalize to 0-1 range
                similarity = normalize_similarity_score(blosum_score)
                
                # Filter by minimum similarity
                if similarity >= min_similarity:
                    # Build match dictionary
                    match = {
                        'target_id': db_site['site_id'],
                        'target_uniprot': db_site['uniprot_id'],
                        'target_site': f"{db_site['residue_type']}{db_site['residue_number']}",
                        'site_type': db_site['residue_type'],
                        'similarity': similarity,
                        'blosum_score': blosum_score,
                        'motif': db_site['motif'],
                        'match_type': 'blosum'
                    }
                    
                    site_matches.append(match)
            
            # Sort by similarity (highest first)
            site_matches.sort(key=lambda x: x['similarity'], reverse=True)
            
            # Limit to top 50 matches
            site_matches = site_matches[:50]
            
            # Add to results if we found matches
            if site_matches:
                blosum_matches[site_name] = site_matches
                logger.info(f"Found {len(site_matches)} BLOSUM matches for site {site_name}")
        
        return blosum_matches
    
    except Exception as e:
        logger.error(f"Error in batch BLOSUM matching: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return {}