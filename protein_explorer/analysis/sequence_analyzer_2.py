"""
Functions for analyzing sequence similarity between phosphorylation sites.

This module handles sequence-based analysis using database queries and provides
visualization-ready outputs. It integrates with supplementary data to include
motif information.
"""

import os
import pandas as pd
import numpy as np
import logging
from typing import Dict, List, Tuple, Optional, Union
from collections import Counter, defaultdict
import re

# Import database functions
from protein_explorer.db import (
    get_phosphosite_data, get_phosphosites_batch,
    find_sequence_matches, find_sequence_matches_batch
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def preload_sequence_data(file_path: str = None) -> None:
    """
    Function kept for backward compatibility.
    Sequence data is now loaded from the database as needed.
    
    Args:
        file_path: Path to the sequence similarity data file (not used)
    """
    logger.info("Using database for sequence similarity data")

def find_sequence_matches(site_id: str, 
                         top_n: int = 200, 
                         min_similarity: float = 0.4) -> List[Dict]:
    """
    Find sequence similarity matches for a site.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        top_n: Maximum number of results to return
        min_similarity: Minimum similarity score to include (0-1)
        
    Returns:
        List of dictionaries with match information
    """
    try:
        # Get sequence matches from database
        matches = from protein_explorer.db import find_sequence_matches
        raw_matches = find_sequence_matches(site_id, min_similarity)
        
        if not raw_matches:
            logger.warning(f"No sequence matches found for {site_id}")
            return []
            
        # Sort by similarity (highest first) and take top N
        matches = sorted(raw_matches, key=lambda x: x['similarity'], reverse=True)
        if top_n:
            matches = matches[:top_n]
            
        # Enhance matches with motif data
        enhanced_matches = enhance_sequence_matches(matches)
        
        return enhanced_matches
    except Exception as e:
        logger.error(f"Error finding sequence matches: {e}")
        return []

def extract_motif_from_site_id(site_id: str) -> Optional[str]:
    """
    Extract motif sequence for a given site ID using the database.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        
    Returns:
        Motif sequence or None if not found
    """
    # Get from supplementary data
    try:
        supp_data = get_phosphosite_data(site_id)
        if supp_data:
            if 'SITE_+/-7_AA' in supp_data and supp_data['SITE_+/-7_AA']:
                return supp_data['SITE_+/-7_AA']
            elif 'motif' in supp_data and supp_data['motif']:
                return supp_data['motif']
    except Exception as e:
        logger.error(f"Error getting motif from supplementary data: {e}")
    
    return None

def enhance_sequence_matches(matches: List[Dict]) -> List[Dict]:
    """
    Enhance sequence matches with additional information from the database.
    
    Args:
        matches: List of match dictionaries from find_sequence_matches
        
    Returns:
        Enhanced list of match dictionaries
    """
    # If no matches, return empty list
    if not matches:
        return []
    
    enhanced_matches = []
    
    # Extract target IDs for batch query
    target_ids = [match['target_id'] for match in matches]
    
    # Batch query for supplementary data
    supp_data_batch = get_phosphosites_batch(target_ids)
    
    for match in matches:
        # Create enhanced match from original
        enhanced = match.copy()
        
        # Add motif if available
        target_id = match['target_id']
        
        # First check if motif is already in the match
        if 'motif' not in enhanced or not enhanced['motif']:
            # If not, try to get from the batch query results
            if target_id in supp_data_batch:
                site_data = supp_data_batch[target_id]
                if 'SITE_+/-7_AA' in site_data and site_data['SITE_+/-7_AA']:
                    enhanced['motif'] = site_data['SITE_+/-7_AA']
                elif 'motif' in site_data and site_data['motif']:
                    enhanced['motif'] = site_data['motif']
            
            # If still not found, try direct query
            if 'motif' not in enhanced or not enhanced['motif']:
                motif = extract_motif_from_site_id(target_id)
                if motif:
                    enhanced['motif'] = motif
        
        enhanced_matches.append(enhanced)
    
    return enhanced_matches

def analyze_motif_conservation(matches: List[Dict], 
                             query_motif: str = None) -> Dict:
    """
    Analyze conservation patterns in motifs of sequence-similar sites.
    
    Args:
        matches: List of match dictionaries from find_sequence_matches
        query_motif: Motif of the query site
        
    Returns:
        Dictionary with conservation analysis results
    """
    logger.info(f"Analyzing conservation with {len(matches)} matches")
    # If no matches, return empty results
    if not matches:
        return {
            'motif_count': 0,
            'conserved_positions': [],
            'n_term_analysis': {},
            'c_term_analysis': {},
            'position_frequencies': {}
        }
    
    # Standardize motif length by assuming phosphosite is in the middle
    # and padding with X if necessary
    std_motifs = []
    for motif in motifs:
        center_pos = len(motif) // 2
        site_char = motif[center_pos]
        
        before_site = motif[:center_pos]
        after_site = motif[center_pos+1:]
        
        # Ensure we have 7 positions before and after
        if len(before_site) < 7:
            before_site = 'X' * (7 - len(before_site)) + before_site
        else:
            before_site = before_site[-7:]
            
        if len(after_site) < 7:
            after_site = after_site + 'X' * (7 - len(after_site))
        else:
            after_site = after_site[:7]
            
        std_motifs.append(before_site + site_char + after_site)
    
    # Analyze conservation at each position
    position_counts = []
    for i in range(15):  # -7 to +7 positions
        counts = Counter()
        for motif in std_motifs:
            if i < len(motif):
                counts[motif[i]] += 1
        position_counts.append(counts)
    
    # Calculate frequency of each amino acid at each position
    position_frequencies = {}
    motif_count = len(std_motifs)
    
    for i, counts in enumerate(position_counts):
        position = i - 7  # Convert to -7 to +7 positions
        position_frequencies[position] = {}
        
        for aa, count in counts.items():
            position_frequencies[position][aa] = count / motif_count
    
    # Identify positions with strong conservation (>50% same AA)
    conserved_positions = []
    for pos, freqs in position_frequencies.items():
        if pos == 0:  # Skip the phosphosite itself
            continue
            
        most_common = max(freqs.items(), key=lambda x: x[1], default=(None, 0))
        if most_common[1] >= 0.5:  # 50% or more conservation
            conserved_positions.append({
                'position': pos,
                'amino_acid': most_common[0],
                'frequency': most_common[1] * 100
            })
    
    # Analyze N-terminal and C-terminal regions separately
    n_term_motifs = [m[:7] for m in std_motifs]  # -7 to -1 positions
    c_term_motifs = [m[8:] for m in std_motifs]  # +1 to +7 positions
    
    # Amino acid group classification
    aa_groups = {
        'polar': 'STYCNQ',
        'nonpolar': 'AVILMFWPG',
        'acidic': 'DE',
        'basic': 'KRH',
        'other': 'X'
    }
    
    # Function to analyze region
    def analyze_region(motifs):
        aa_composition = defaultdict(int)
        aa_group_composition = defaultdict(int)
        total_aa = len(motifs) * 7  # 7 positions per motif
        
        for motif in motifs:
            for aa in motif:
                aa_composition[aa] += 1
                
                # Classify by group
                for group, members in aa_groups.items():
                    if aa in members:
                        aa_group_composition[group] += 1
                        break
        
        # Convert to percentages
        aa_percentages = {aa: count/total_aa*100 for aa, count in aa_composition.items()}
        group_percentages = {group: count/total_aa*100 for group, count in aa_group_composition.items()}
        
        return {
            'aa_composition': dict(sorted(aa_percentages.items(), key=lambda x: x[1], reverse=True)),
            'group_composition': dict(sorted(group_percentages.items(), key=lambda x: x[1], reverse=True))
        }
    
    n_term_analysis = analyze_region(n_term_motifs)
    c_term_analysis = analyze_region(c_term_motifs)
    
    # Calculate a consensus motif
    consensus_motif = ""
    for i in range(15):  # -7 to +7 positions
        if i == 7:  # Phosphosite position
            consensus_motif += std_motifs[0][7] if std_motifs else "X"
            continue
            
        counts = position_counts[i]
        if counts:
            # Get the most common AA, but exclude X unless it's the only one
            filtered_counts = {aa: count for aa, count in counts.items() if aa != 'X'}
            if filtered_counts:
                most_common = max(filtered_counts.items(), key=lambda x: x[1])[0]
            else:
                most_common = 'X'
            consensus_motif += most_common
        else:
            consensus_motif += "X"
    
    return {
        'motif_count': len(motifs),
        'consensus_motif': consensus_motif,
        'conserved_positions': conserved_positions,
        'n_term_analysis': n_term_analysis,
        'c_term_analysis': c_term_analysis,
        'position_frequencies': position_frequencies
    }

def create_sequence_network_data(query_site_id: str, 
                               matches: List[Dict],
                               query_motif: str = None) -> Dict:
    """
    Create data for sequence similarity network visualization.
    
    Args:
        query_site_id: Site ID of the query
        matches: List of match dictionaries from find_sequence_matches
        query_motif: Motif of the query site
        
    Returns:
        Dictionary with nodes and links for network visualization
    """
    if not matches:
        return {'nodes': [], 'links': []}
    
    # Extract parts of query_site_id
    query_parts = query_site_id.split('_')
    query_uniprot = query_parts[0] if len(query_parts) > 0 else ""
    query_site = query_parts[1] if len(query_parts) > 1 else ""
    
    # Extract site type if possible
    query_site_type = None
    site_match = re.match(r'([STY])(\d+)', query_site)
    if site_match:
        query_site_type = site_match.group(1)
        query_site_number = site_match.group(2)
        query_site = f"{query_site_type}{query_site_number}"
    
    # Create nodes list starting with query node
    nodes = [{
        'id': query_site_id,
        'name': query_site,
        'display_name': query_site_id,
        'uniprot': query_uniprot,
        'type': 'query',
        'site_type': query_site_type,
        'size': 12,
        'motif': query_motif
    }]
    
    # Create links list
    links = []
    
    # Process each match
    seen_nodes = {query_site_id}  # Track nodes we've already added
    
    for match in matches:
        target_id = match['target_id']
        
        # Skip self-matches
        if target_id == query_site_id:
            continue
            
        # Skip duplicates
        if target_id in seen_nodes:
            continue
            
        seen_nodes.add(target_id)
        
        # Create node for this match
        nodes.append({
            'id': target_id,
            'name': match['target_site'],
            'display_name': target_id,
            'uniprot': match['target_uniprot'],
            'type': 'target',
            'site_type': match.get('site_type'),
            'similarity': match['similarity'],
            'size': 8,
            'motif': match.get('motif')
        })
        
        # Create link between query and target
        links.append({
            'source': query_site_id,
            'target': target_id,
            'similarity': match['similarity']
        })
    
    return {
        'nodes': nodes,
        'links': links
    }

def get_motif_enrichment(matches: List[Dict], 
                       query_motif: str = None,
                       background_frequencies: Dict = None) -> Dict:
    """
    Calculate amino acid enrichment in motifs compared to background frequencies.
    
    Args:
        matches: List of match dictionaries
        query_motif: Motif of the query site
        background_frequencies: Background AA frequencies (defaults to UniProt averages)
        
    Returns:
        Dictionary with enrichment analysis results
    """
    # Default background frequencies from UniProt
    if background_frequencies is None:
        background_frequencies = {
            'A': 0.0825, 'R': 0.0553, 'N': 0.0406, 'D': 0.0545,
            'C': 0.0137, 'Q': 0.0393, 'E': 0.0675, 'G': 0.0707,
            'H': 0.0227, 'I': 0.0595, 'L': 0.0965, 'K': 0.0584,
            'M': 0.0241, 'F': 0.0386, 'P': 0.0470, 'S': 0.0656,
            'T': 0.0534, 'W': 0.0108, 'Y': 0.0292, 'V': 0.0687
        }
    
    # Collect motifs
    motifs = []
    for match in matches:
        if 'motif' in match and match['motif']:
            motifs.append(match['motif'])
    
    # Add query motif if provided
    if query_motif:
        motifs = [query_motif] + motifs
    
    # If no motifs, return empty results
    if not motifs:
        return {
            'position_enrichment': {},
            'overall_enrichment': {}
        }
    
    # Standardize motifs
    std_motifs = []
    for motif in motifs:
        center_pos = len(motif) // 2
        site_char = motif[center_pos]
        
        before_site = motif[:center_pos]
        after_site = motif[center_pos+1:]
        
        # Ensure we have 7 positions before and after
        if len(before_site) < 7:
            before_site = 'X' * (7 - len(before_site)) + before_site
        else:
            before_site = before_site[-7:]
            
        if len(after_site) < 7:
            after_site = after_site + 'X' * (7 - len(after_site))
        else:
            after_site = after_site[:7]
            
        std_motifs.append(before_site + site_char + after_site)
    
    # Count AAs at each position
    position_counts = []
    for i in range(15):  # -7 to +7 positions
        counts = Counter()
        for motif in std_motifs:
            if i < len(motif) and motif[i] != 'X':  # Skip placeholder X
                counts[motif[i]] += 1
        position_counts.append(counts)
    
    # Count overall AA frequencies (excluding phosphosite and X)
    overall_counts = Counter()
    for motif in std_motifs:
        for i, aa in enumerate(motif):
            if i != 7 and aa != 'X':  # Skip phosphosite and placeholder X
                overall_counts[aa] += 1
    
    # Calculate position-specific enrichment
    position_enrichment = {}
    for i, counts in enumerate(position_counts):
        position = i - 7  # Convert to -7 to +7 positions
        
        if position == 0:  # Skip phosphosite position
            continue
            
        position_enrichment[position] = {}
        total_aas = sum(counts.values())
        
        if total_aas == 0:
            continue
            
        for aa, count in counts.items():
            if aa == 'X' or aa not in background_frequencies:
                continue
                
            observed_freq = count / total_aas
            expected_freq = background_frequencies[aa]
            enrichment = observed_freq / expected_freq if expected_freq > 0 else 0
            
            position_enrichment[position][aa] = {
                'observed_freq': observed_freq,
                'expected_freq': expected_freq,
                'enrichment': enrichment,
                'count': count
            }
    
    # Calculate overall enrichment
    overall_enrichment = {}
    total_aas = sum(overall_counts.values())
    
    if total_aas > 0:
        for aa, count in overall_counts.items():
            if aa not in background_frequencies:
                continue
                
            observed_freq = count / total_aas
            expected_freq = background_frequencies[aa]
            enrichment = observed_freq / expected_freq if expected_freq > 0 else 0
            
            overall_enrichment[aa] = {
                'observed_freq': observed_freq,
                'expected_freq': expected_freq,
                'enrichment': enrichment,
                'count': count
            }
    
    return {
        'position_enrichment': position_enrichment,
        'overall_enrichment': overall_enrichment
    }

def create_sequence_motif_visualization(query_site_id: str, 
                                         query_motif: str,
                                         matches: List[Dict],
                                         max_matches: int = 10) -> str:
    """
    Create HTML for a comparative visualization of sequence motifs.
    
    Args:
        query_site_id: ID of the query site
        query_motif: Motif of the query site
        matches: List of match dictionaries
        max_matches: Maximum number of matches to display
        
    Returns:
        HTML code for the visualization
    """
    # This function generates HTML and doesn't depend on database access,
    # so it doesn't need to be refactored.': [],
            'n_term_analysis': {},
            'c_term_analysis': {},
            'position_frequencies': {}
        }
    
    # Collect motifs from matches
    motifs = []
    for match in matches:
        motif = match.get('motif')
        if motif:
            motifs.append(motif)
    
    # Add query motif if provided
    if query_motif:
        motifs = [query_motif] + motifs
    
    # If no motifs found, return empty results
    if not motifs:
        return {
            'motif_count': 0,
            'conserved_positions