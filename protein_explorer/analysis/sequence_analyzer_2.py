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
from protein_explorer.db.db import (
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
        from protein_explorer.db import find_sequence_matches
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
    Makes sure site types and motif sequences are consistent.
    
    Args:
        matches: List of match dictionaries from find_sequence_matches
        
    Returns:
        Enhanced list of match dictionaries
    """
    # Configure logging
    import logging
    logger = logging.getLogger(__name__)
    
    # If no matches, return empty list
    if not matches:
        logger.info("No matches provided to enhance_sequence_matches")
        return []
    
    logger.info(f"Enhancing {len(matches)} sequence matches")
    
    # Debug: Check first match structure
    if matches and len(matches) > 0:
        logger.info(f"Sample match keys: {list(matches[0].keys())}")
    
    enhanced_matches = []
    
    # Extract target IDs for batch query
    target_ids = []
    for match in matches:
        if 'target_id' in match and match['target_id']:
            target_ids.append(match['target_id'])
    
    logger.info(f"Found {len(target_ids)} target IDs for batch query")
    
    # Batch query for supplementary data if we have any target IDs
    supp_data_batch = {}
    if target_ids:
        try:
            supp_data_batch = get_phosphosites_batch(target_ids)
            logger.info(f"Retrieved supplementary data for {len(supp_data_batch)} sites")
        except Exception as e:
            logger.error(f"Error getting phosphosite batch data: {e}")
    
    # Process each match individually with careful error handling
    for match_index, match in enumerate(matches):
        try:
            # Create enhanced match from original
            enhanced = match.copy()
            
            # Add motif if available
            target_id = match.get('target_id')
            
            # First check if motif is already in the match
            if ('motif' not in enhanced or enhanced.get('motif') is None or 
                    (isinstance(enhanced.get('motif'), str) and not enhanced['motif'].strip())):
                # If not, try to get from the batch query results
                if target_id and target_id in supp_data_batch:
                    site_data = supp_data_batch[target_id]
                    
                    if 'SITE_+/-7_AA' in site_data and site_data.get('SITE_+/-7_AA'):
                        enhanced['motif'] = str(site_data['SITE_+/-7_AA']).upper()  # Ensure uppercase
                    elif 'motif' in site_data and site_data.get('motif'):
                        enhanced['motif'] = str(site_data['motif']).upper()  # Ensure uppercase
            # If motif exists but isn't uppercase, make sure it is
            elif enhanced.get('motif') and isinstance(enhanced['motif'], str):
                enhanced['motif'] = enhanced['motif'].upper()
                
            # Check for site type in Residue field from supplementary data
            if target_id and target_id in supp_data_batch:
                site_data = supp_data_batch[target_id]
                print("HERE IS THE SITE DATA, IS RESIDUE COLUMN PRESENT?")
                print(site_data)
                # Get site type from Residue/residue column if available
                if 'Residue' in site_data and site_data.get('Residue') in 'STY':
                    enhanced['site_type'] = site_data['Residue']
                    logger.debug(f"Using site type {enhanced['site_type']} from Residue column for {target_id}")
                elif 'residue' in site_data and site_data.get('residue') in 'STY':
                    enhanced['site_type'] = site_data['residue']
                    logger.debug(f"Using site type {enhanced['site_type']} from residue column for {target_id}")
                    
                # Add is_known_phosphosite if available
                if 'is_known_phosphosite' in site_data:
                    try:
                        enhanced['is_known_phosphosite'] = float(site_data['is_known_phosphosite'])
                    except (ValueError, TypeError):
                        enhanced['is_known_phosphosite'] = 0.0  # Default to not known
                
                # Add other useful fields
                for field in ['site_plddt', 'mean_plddt', 'nearby_count', 
                             'surface_accessibility', 'acidic_percentage', 
                             'basic_percentage', 'aromatic_percentage']:
                    if field in site_data and site_data[field] is not None:
                        try:
                            enhanced[field] = site_data[field]
                        except Exception as e:
                            logger.warning(f"Error adding field {field}: {e}")
                
            # IMPORTANT: Leave target_site as-is - don't add 'S' prefix
            # Just make sure it exists
            if 'target_site' not in enhanced or enhanced.get('target_site') is None:
                # Extract site number from target_id if possible
                if target_id and '_' in target_id:
                    site_number = target_id.split('_')[1]
                    enhanced['target_site'] = site_number
                    logger.debug(f"Set target_site to {site_number} from target_id for {target_id}")
            
            enhanced_matches.append(enhanced)
            
        except Exception as e:
            logger.error(f"Error enhancing match {match_index}: {e}")
            # Add the original match without enhancement rather than skipping it
            enhanced_matches.append(match)
    
    logger.info(f"Finished enhancing {len(enhanced_matches)} sequence matches")
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
    
    # Collect motifs from matches
    motifs = []
    for match in matches:
        if 'motif' in match and match['motif']:
            motifs.append(match['motif'])
    
    # Add query motif if provided
    if query_motif:
        motifs = [query_motif] + motifs
    
    # If no motifs found, return empty results
    if not motifs:
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
    # so it doesn't need to be significantly refactored.
    
    logger.info(f"Creating motif comparison with {len(matches)} matches")
    logger.info(f"Query motif: {query_motif}")

    # Extract uniprot and site from query_site_id
    query_parts = query_site_id.split('_')
    query_uniprot = query_parts[0] if len(query_parts) > 0 else ""
    query_site = query_parts[1] if len(query_parts) > 1 else ""
    
    site_match = re.match(r'([STY])(\d+)', query_site)
    if site_match:
        query_site = site_match.group(0)
    
    # If the query motif is empty or None, try to get it from supplementary data
    if not query_motif:
        try:
            supp_data = get_phosphosite_data(query_site_id)
            if supp_data:
                if 'SITE_+/-7_AA' in supp_data and supp_data['SITE_+/-7_AA']:
                    query_motif = supp_data['SITE_+/-7_AA']
                elif 'motif' in supp_data and supp_data['motif']:
                    query_motif = supp_data['motif']
                logger.info(f"Retrieved query motif from supplementary data: {query_motif}")
        except Exception as e:
            logger.error(f"Error retrieving query motif from supplementary data: {e}")
    
    # Filter matches with motifs
    valid_matches = []
    match_ids_without_motifs = []
    
    for match in matches:
        # Check if match already has a motif
        if 'motif' in match and match['motif']:
            valid_matches.append(match)
        else:
            match_ids_without_motifs.append(match['target_id'])
    
    # Batch query for motifs if needed
    if match_ids_without_motifs:
        try:
            # Get supplementary data for matches without motifs
            supp_data_batch = get_phosphosites_batch(match_ids_without_motifs)
            
            # Add motifs to matches
            for match in matches:
                if 'motif' not in match or not match['motif']:
                    target_id = match['target_id']
                    if target_id in supp_data_batch:
                        site_data = supp_data_batch[target_id]
                        
                        # Create a new match with the motif added
                        enhanced_match = match.copy()
                        
                        if 'SITE_+/-7_AA' in site_data and site_data['SITE_+/-7_AA']:
                            enhanced_match['motif'] = site_data['SITE_+/-7_AA']
                            valid_matches.append(enhanced_match)
                        elif 'motif' in site_data and site_data['motif']:
                            enhanced_match['motif'] = site_data['motif']
                            valid_matches.append(enhanced_match)
        except Exception as e:
            logger.error(f"Error retrieving motifs from supplementary data: {e}")
    
    logger.info(f"Found {len(valid_matches)} matches with motifs after enhancement")

    # If no valid matches, return simple message
    if not valid_matches:
        return f"""
        <div class="alert alert-info">
            No sequence motif data available for comparison with {query_site_id}.
        </div>
        """
    
    # Sort by similarity (highest first)
    sorted_matches = sorted(valid_matches, key=lambda x: x.get('similarity', 0), reverse=True)
    
    # Take top N matches
    top_matches = sorted_matches[:max_matches]
    
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
    
    # Standardization function for motifs
    def standardize_motif(motif):
        # Find the center position (phosphosite)
        center_pos = len(motif) // 2
        
        # Get the phosphosite and parts before/after
        site_char = motif[center_pos]
        before_site = motif[:center_pos]
        after_site = motif[center_pos + 1:]
        
        # Ensure we have exactly 7 characters before and after
        if len(before_site) < 7:
            before_site = "X" * (7 - len(before_site)) + before_site
        else:
            before_site = before_site[-7:]
            
        if len(after_site) < 7:
            after_site = after_site + "X" * (7 - len(after_site))
        else:
            after_site = after_site[:7]
            
        return before_site + site_char + after_site
    
    # Function to trim to just -5 to +5 range for display
    def trim_to_central_range(motif_str):
        # We want to keep positions 2-12 (0-indexed) from a 15-char motif
        # which represent positions -5 to +5 around the phosphosite
        return motif_str[2:13]
    
    # Modified helper function to create HTML for a motif
    def create_motif_html(motif):
        # First standardize to full 15 chars
        std_motif = standardize_motif(motif)
        
        # Then trim to -5 to +5 range
        trimmed_motif = trim_to_central_range(std_motif)
        
        # Create HTML
        html = '<div class="motif-sequence" style="display: flex; flex-wrap: nowrap;">'
        for i, aa in enumerate(trimmed_motif):
            aa_class = get_aa_class(aa)
            highlight_class = "highlighted" if i == 5 else aa_class  # Position 5 is the phosphosite in trimmed motif
            html += f'<div class="motif-aa {highlight_class}" style="width: 24px; height: 24px; display: flex; align-items: center; justify-content: center; margin: 0 1px; border-radius: 3px;">{aa}</div>'
        html += '</div>'
        return html
    
    # Create HTML for the visualization
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
        <h5 class="mb-3">Sequence Similarity Motif Comparison</h5>
        
        <!-- Position markers -->
        <div class="motif-position">
    """
    
    # Add position markers from -5 to +5
    for i in range(-5, 6):
        html += f'<span>{i}</span>'
    
    html += """
        </div>
    """
    
    # Add query motif row
    html += f"""
        <div class="motif-row">
            <div class="motif-label">{query_uniprot}_{query_site}:</div>
            {create_motif_html(query_motif)}
            <div class="match-info">
                Query
            </div>
        </div>
    """
    
    # Add match motifs
    for match in top_matches:
        motif = match.get('motif', '')
        target_site = match.get('target_site', 'Unknown')
        target_uniprot = match.get('target_uniprot', 'Unknown')
        similarity = match.get('similarity', 0.0)
        
        html += f"""
        <div class="motif-row">
            <div class="motif-label">{target_uniprot}_{target_site}:</div>
            {create_motif_html(motif)}
            <div class="match-info">
                Similarity: {similarity:.2f} | <a href="/site/{target_uniprot}/{target_site}" class="text-decoration-none">View site</a>
            </div>
        </div>
        """
    
    html += """
    </div>
    """
    
    return html

# Function to batch find sequence matches for multiple sites
def find_sequence_matches_batch(site_ids: List[str], 
                              min_similarity: float = 0.4) -> Dict[str, List[Dict]]:
    """
    Find sequence similarity matches for multiple sites in a batch.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        min_similarity: Minimum similarity score to include (0-1)
        
    Returns:
        Dictionary mapping site IDs to lists of match dictionaries
    """
    try:
        # Get sequence matches from database in batch
        from protein_explorer.db import find_sequence_matches_batch
        raw_matches_dict = find_sequence_matches_batch(site_ids, min_similarity)
        
        if not raw_matches_dict:
            logger.warning(f"No sequence matches found for any of the {len(site_ids)} sites")
            return {}
            
        # Collect all target IDs for batch motif lookup
        all_target_ids = set()
        for site_matches in raw_matches_dict.values():
            for match in site_matches:
                all_target_ids.add(match['target_id'])
        
        # Batch query for motifs
        supp_data_batch = get_phosphosites_batch(list(all_target_ids))
        
        # Process each site's matches
        result = {}
        for site_id, matches in raw_matches_dict.items():
            # Sort by similarity (highest first)
            sorted_matches = sorted(matches, key=lambda x: x['similarity'], reverse=True)
            
            # Enhance with motifs from batch data
            enhanced_matches = []
            for match in sorted_matches:
                enhanced = match.copy()
                target_id = match['target_id']
                
                # Add motif if available in batch data
                if target_id in supp_data_batch:
                    site_data = supp_data_batch[target_id]
                    if 'SITE_+/-7_AA' in site_data and site_data['SITE_+/-7_AA']:
                        enhanced['motif'] = site_data['SITE_+/-7_AA']
                    elif 'motif' in site_data and site_data['motif']:
                        enhanced['motif'] = site_data['motif']
                
                enhanced_matches.append(enhanced)
            
            result[site_id] = enhanced_matches
        
        return result
    except Exception as e:
        logger.error(f"Error in batch sequence match finding: {e}")
        return {}
    


# This function should be added to protein_explorer/analysis/sequence_analyzer_2.py
def find_sequence_matches_with_connections(site_ids: List[str], min_similarity: float = 0.4) -> Dict:
    """
    Find sequence similarity matches for multiple sites in a batch,
    including connections between the similar sites.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        min_similarity: Minimum similarity score to include (0-1)
        
    Returns:
        Dictionary with nodes, edges and required data for network visualization
    """
    if not site_ids:
        return {'nodes': [], 'links': [], 'site_matches': {}}
    
    # Get batch match results for the sites
    site_matches_dict = find_sequence_matches_batch(site_ids, min_similarity)
    
    # Set to track all unique match IDs for the next query
    all_match_ids = set()
    
    # Collect all matches for each site
    for site_id, matches in site_matches_dict.items():
        for match in matches:
            all_match_ids.add(match['target_id'])
    
    # Find interconnections between the matches themselves
    # This requires an additional batch query
    match_interconnections = {}
    if all_match_ids:
        # Convert to list for the query
        match_id_list = list(all_match_ids)
        
        # Find matches among the matches themselves
        match_matches = find_sequence_matches_batch(match_id_list, min_similarity)
        
        # Store valid interconnections
        for match_id, connections in match_matches.items():
            valid_connections = []
            for connection in connections:
                # Only include connections to other matches we already know about
                # Filter by our threshold and avoid self-connections
                if (connection['target_id'] in all_match_ids and 
                    connection['target_id'] != match_id and
                    connection['similarity'] >= min_similarity):
                    valid_connections.append(connection)
            
            if valid_connections:
                match_interconnections[match_id] = valid_connections
    
    # Build network with nodes and links
    nodes = []
    links = []
    node_map = {}  # Track nodes we've added
    link_map = {}  # Track links we've added to avoid duplicates
    
    # Get all supplementary data in one batch query
    # Collect all site IDs that need supplementary data (original sites + matches)
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
            target_id = match['target_id']
            
            # Skip if already added
            if target_id in node_map:
                continue
                
            # Extract target info
            target_uniprot = match['target_uniprot']
            target_site = match['target_site']
            similarity = match['similarity']
            
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
                'siteType': match.get('site_type', 'S'),
                'isKnown': False,  # Assume false for matches
                'similarity': similarity,
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
                'similarity': similarity
            })
            
            link_map[link_id] = True
            # Also add reverse direction to avoid duplicate links
            link_map[f"{target_id}-{site_id}"] = True
    
    # Finally add links between matches
    for match_id, connections in match_interconnections.items():
        for connection in connections:
            target_id = connection['target_id']
            similarity = connection['similarity']
            
            # Skip if link already added
            link_id = f"{match_id}-{target_id}"
            if link_id in link_map:
                continue
                
            links.append({
                'source': match_id,
                'target': target_id,
                'similarity': similarity
            })
            
            link_map[link_id] = True
            # Also add reverse direction to avoid duplicate links
            link_map[f"{target_id}-{match_id}"] = True
    
    return {
        'nodes': nodes,
        'links': links,
        'site_matches': site_matches_dict  # Include original match data for reference
    }