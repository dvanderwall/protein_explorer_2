
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, Selection, Vector, MMCIFParser, MMCIF2Dict
import argparse
import warnings
from math import degrees
import re
import time
import gzip
from Bio.PDB.Polypeptide import three_to_index, is_aa
from scipy.spatial import cKDTree  # Fast neighbor search
import os

# Add this function to replace three_to_one
def three_to_one(res_name):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    try:
        index = three_to_index(res_name)
        return amino_acids[index]
    except:
        return "X"

# Suppress warnings
warnings.filterwarnings("ignore")

# Define amino acid properties
aa_properties = {
    'A': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'small'},
    'R': {'hydrophobic': False, 'polar': True, 'charged': True, 'charge': 1, 'size': 'large'},
    'N': {'hydrophobic': False, 'polar': True, 'charged': False, 'charge': 0, 'size': 'medium'},
    'D': {'hydrophobic': False, 'polar': True, 'charged': True, 'charge': -1, 'size': 'medium'},
    'C': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'medium'},
    'Q': {'hydrophobic': False, 'polar': True, 'charged': False, 'charge': 0, 'size': 'medium'},
    'E': {'hydrophobic': False, 'polar': True, 'charged': True, 'charge': -1, 'size': 'medium'},
    'G': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'small'},
    'H': {'hydrophobic': False, 'polar': True, 'charged': True, 'charge': 0.5, 'size': 'medium'},
    'I': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'large'},
    'L': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'large'},
    'K': {'hydrophobic': False, 'polar': True, 'charged': True, 'charge': 1, 'size': 'large'},
    'M': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'large'},
    'F': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'large'},
    'P': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'medium'},
    'S': {'hydrophobic': False, 'polar': True, 'charged': False, 'charge': 0, 'size': 'small'},
    'T': {'hydrophobic': False, 'polar': True, 'charged': False, 'charge': 0, 'size': 'medium'},
    'W': {'hydrophobic': True, 'polar': True, 'charged': False, 'charge': 0, 'size': 'large'},
    'Y': {'hydrophobic': False, 'polar': True, 'charged': False, 'charge': 0, 'size': 'large'},
    'V': {'hydrophobic': True, 'polar': False, 'charged': False, 'charge': 0, 'size': 'medium'},
}

# Precompute some common values
STY_RESIDUES = set(['SER', 'THR', 'TYR'])
STY_SINGLE_LETTER = set(['S', 'T', 'Y'])
MAX_DISTANCE_SQUARED = 16  # 4Å squared


def parse_plddt_from_structure_file(structure_file):
    """Extract pLDDT scores from AlphaFold structure files (B-factor column)
    Works with both .pdb, .cif and their gzipped versions"""
    plddt_dict = {}
    
    # Check if file is gzipped
    is_gzipped = structure_file.endswith('.gz')
    
    # Check if file is CIF
    is_cif = structure_file.endswith('.cif') or structure_file.endswith('.cif.gz')
    
    # For CIF files, try to extract B-factors using MMCIF2Dict
    if is_cif:
        try:
            if is_gzipped:
                with gzip.open(structure_file, 'rt') as f:
                    mmcif_dict = MMCIF2Dict.MMCIF2Dict(f)
            else:
                mmcif_dict = MMCIF2Dict.MMCIF2Dict(structure_file)
            
            # Check if the necessary data is in the CIF dictionary
            if '_atom_site.B_iso_or_equiv' in mmcif_dict and '_atom_site.label_seq_id' in mmcif_dict:
                b_factors = mmcif_dict['_atom_site.B_iso_or_equiv']
                seq_ids = mmcif_dict['_atom_site.label_seq_id']
                
                # Convert values to appropriate types
                b_factors = [float(b) for b in b_factors]
                seq_ids = [int(id) if id != '?' else None for id in seq_ids]
                
                # Create a dictionary mapping residue IDs to average B-factors
                b_factor_sums = {}
                b_factor_counts = {}
                
                for i in range(len(seq_ids)):
                    if seq_ids[i] is not None:
                        if seq_ids[i] not in b_factor_sums:
                            b_factor_sums[seq_ids[i]] = 0
                            b_factor_counts[seq_ids[i]] = 0
                        b_factor_sums[seq_ids[i]] += b_factors[i]
                        b_factor_counts[seq_ids[i]] += 1
                
                # Calculate average B-factor for each residue
                for res_id in b_factor_sums:
                    plddt_dict[res_id] = b_factor_sums[res_id] / b_factor_counts[res_id]
                
                return plddt_dict
        except Exception as e:
            print(f"  Warning: Could not extract B-factors from CIF dictionary: {e}")
            print("  Falling back to text parsing...")
    
    # Fall back to text parsing for both PDB and CIF
    # Open file accordingly
    if is_gzipped:
        opener = gzip.open(structure_file, 'rt')
    else:
        opener = open(structure_file, 'r')
        
    with opener as f:
        for line in f:
            if line.startswith('ATOM'):
                res_id = int(line[22:26].strip())
                b_factor = float(line[60:66].strip())
                plddt_dict[res_id] = b_factor
    
    return plddt_dict

def find_hydrogen_bonds_vectorized(hydroxyl_atom, neighbor_residues, max_distance=3.5):
    """Vectorized calculation of hydrogen bonds"""
    if not hydroxyl_atom:
        return 0
    
    max_distance_squared = max_distance * max_distance
    
    # Extract all oxygen and nitrogen atoms from neighbor residues
    potential_partners = []
    
    for res in neighbor_residues:
        # Skip the residue itself
        if res is hydroxyl_atom.get_parent():
            continue
            
        for atom in res:
            # Only consider oxygen and nitrogen atoms as potential partners
            if atom.element in ['O', 'N']:
                potential_partners.append(atom.coord)
    
    if not potential_partners:
        return 0
    
    # Convert to numpy arrays
    potential_partners = np.array(potential_partners)
    hydroxyl_coord = np.array(hydroxyl_atom.coord)
    
    # Calculate squared distances
    vectors = potential_partners - hydroxyl_coord
    distances_sq = np.sum(vectors * vectors, axis=1)
    
    # Count potential hydrogen bonds
    h_bonds = np.sum(distances_sq <= max_distance_squared)
    
    return h_bonds

def calculate_hydroxyl_exposure_vectorized(hydroxyl_atom, neighbor_residues):
    """Vectorized calculation of hydroxyl exposure"""
    if not hydroxyl_atom or not neighbor_residues:
        return None, 0
    
    # Extract all atom coordinates from neighbor residues
    neighbor_atoms = []
    backbone_atoms = []
    
    for res in neighbor_residues:
        for atom in res:
            # Skip if this is from the same residue as hydroxyl_atom
            if atom.get_parent() is hydroxyl_atom.get_parent():
                continue
                
            neighbor_atoms.append(atom.coord)
            
            # Check if it's a backbone atom
            if atom.name in ['N', 'CA', 'C', 'O']:
                backbone_atoms.append(atom.coord)
    
    if not neighbor_atoms:
        return 1.0, 0  # Fully exposed if no neighbors
    
    # Convert to numpy arrays
    neighbor_atoms = np.array(neighbor_atoms)
    hydroxyl_coord = np.array(hydroxyl_atom.coord)
    
    # Calculate squared distances to hydroxyl
    # Use broadcasting to calculate all distances at once
    vectors = neighbor_atoms - hydroxyl_coord
    distances_sq = np.sum(vectors * vectors, axis=1)
    
    # Count atoms within cutoff distance
    nearby_atoms = np.sum(distances_sq <= MAX_DISTANCE_SQUARED)
    
    # Calculate backbone contacts if we have backbone atoms
    if backbone_atoms:
        backbone_atoms = np.array(backbone_atoms)
        backbone_vectors = backbone_atoms - hydroxyl_coord
        backbone_distances_sq = np.sum(backbone_vectors * backbone_vectors, axis=1)
        nearby_backbone = np.sum(backbone_distances_sq <= MAX_DISTANCE_SQUARED)
    else:
        nearby_backbone = 0
    
    # Calculate exposure
    exposure = 1.0 / (1.0 + nearby_atoms/10.0)
    
    return exposure, nearby_backbone


def find_neighbors_fast(chain_data, residue, radius=10.0):
    """Use KD-tree for fast neighbor finding"""
    # Get CA atom coordinates for the central residue
    try:
        ca_atom = residue['CA']
    except KeyError:
        return []
    
    ca_coord = ca_atom.coord
    
    # Query KD-tree for neighbors
    if chain_data['kd_tree'] is None:
        return []
    
    # Find all atoms within radius
    neighbor_indices = chain_data['kd_tree'].query_ball_point(ca_coord, radius)
    
    # Get unique residues
    neighbor_residues = set()
    for idx in neighbor_indices:
        if idx < len(chain_data['atom_to_residue']):
            neighbor_residues.add(chain_data['atom_to_residue'][idx])
    
    # Remove the center residue from neighbors
    if residue in neighbor_residues:
        neighbor_residues.remove(residue)
    
    return list(neighbor_residues)



def build_chain_data(chain):
    """Pre-compute and cache chain-level data for reuse"""
    chain_data = {}
    
    # Build complete chain sequence once
    chain_seq = ''
    seq_ids = []
    residue_map = {}  # Map residue ID to position in chain
    
    for i, res in enumerate(chain):
        if is_aa(res):
            chain_seq += three_to_one(res.resname)
            seq_ids.append(res.id[1])
            residue_map[res.id[1]] = i
    
    chain_data['sequence'] = chain_seq
    chain_data['seq_ids'] = seq_ids
    chain_data['residue_map'] = residue_map
    
    # Collect atom coordinates for KD-tree
    atom_coords = []
    atom_to_residue = []  # Map atom index to residue
    
    # Also collect STY residues
    sty_residues = []
    
    for res in chain:
        if is_aa(res):
            if res.resname in STY_RESIDUES:
                sty_residues.append(res)
            
            for atom in res:
                atom_coords.append(atom.coord)
                atom_to_residue.append(res)
    
    chain_data['atom_coords'] = np.array(atom_coords)
    chain_data['atom_to_residue'] = atom_to_residue
    chain_data['sty_residues'] = sty_residues
    
    # Build KD-tree for fast neighbor lookup if we have atoms
    if atom_coords:
        chain_data['kd_tree'] = cKDTree(chain_data['atom_coords'])
    else:
        chain_data['kd_tree'] = None
    
    return chain_data


def calc_dihedral(v1, v2, v3, v4):
    """Calculate dihedral angle between 4 vectors"""
    try:
        ab = v1 - v2
        cb = v3 - v2
        db = v4 - v3
        
        u = cb.cross(ab)
        v = db.cross(cb)
        
        u.normalize()
        v.normalize()
        w = u.cross(v)
        
        angle = np.arctan2(w * cb, u * v)
        return angle
    except:
        return None


def get_sequence_motif(chain, residue, window=7):
    """Get the -7:+7 sequence motif around a residue"""
    res_id = residue.id[1]
    chain_seq = ''
    seq_ids = []
    
    # Build ordered sequence and track residue IDs
    for res in chain:
        if is_aa(res):
            chain_seq += three_to_one(res.resname)
            seq_ids.append(res.id[1])
    
    # Find the position in the sequence
    if res_id not in seq_ids:
        return "X" * (2 * window + 1)
    
    seq_pos = seq_ids.index(res_id)
    
    # Extract motif
    start = max(0, seq_pos - window)
    end = min(len(chain_seq), seq_pos + window + 1)
    
    # Pad with X if needed
    prefix = "X" * max(0, window - seq_pos)
    suffix = "X" * max(0, window - (len(chain_seq) - seq_pos - 1))
    
    motif = prefix + chain_seq[start:end] + suffix
    
    # Ensure motif is the correct length
    if len(motif) > 2 * window + 1:
        motif = motif[:2 * window + 1]
    elif len(motif) < 2 * window + 1:
        motif = motif + "X" * (2 * window + 1 - len(motif))
    
    return motif


def calculate_hydroxyl_exposure_vectorized(hydroxyl_atom, neighbors):
    """Vectorized calculation of hydroxyl exposure"""
    if not hydroxyl_atom or not neighbors:
        return None, 0
    
    # Extract all atom coordinates from neighbor residues
    neighbor_atoms = []
    backbone_atoms = []
    
    for res in neighbors:
        for atom in res:
            # Skip if this is from the same residue as hydroxyl_atom
            if atom.get_parent() is hydroxyl_atom.get_parent():
                continue
                
            neighbor_atoms.append(atom.coord)
            
            # Check if it's a backbone atom
            if atom.name in ['N', 'CA', 'C', 'O']:
                backbone_atoms.append(atom.coord)
    
    if not neighbor_atoms:
        return 1.0, 0  # Fully exposed if no neighbors
    
    # Convert to numpy arrays
    neighbor_atoms = np.array(neighbor_atoms)
    hydroxyl_coord = np.array(hydroxyl_atom.coord)
    
    # Calculate squared distances to hydroxyl
    # Use broadcasting to calculate all distances at once
    vectors = neighbor_atoms - hydroxyl_coord
    distances_sq = np.sum(vectors * vectors, axis=1)
    
    # Count atoms within cutoff distance
    nearby_atoms = np.sum(distances_sq <= MAX_DISTANCE_SQUARED)
    
    # Calculate backbone contacts if we have backbone atoms
    if backbone_atoms:
        backbone_atoms = np.array(backbone_atoms)
        backbone_vectors = backbone_atoms - hydroxyl_coord
        backbone_distances_sq = np.sum(backbone_vectors * backbone_vectors, axis=1)
        nearby_backbone = np.sum(backbone_distances_sq <= MAX_DISTANCE_SQUARED)
    else:
        nearby_backbone = 0
    
    # Calculate exposure
    exposure = 1.0 / (1.0 + nearby_atoms/10.0)
    
    return exposure, nearby_backbone

def analyze_pocket_properties(neighbors):
    """Analyze biochemical properties of residues in the pocket"""
    if not neighbors:
        return {
            'hydrophobic_count': 0,
            'polar_count': 0,
            'positive_count': 0,
            'negative_count': 0,
            'net_charge': 0,
            'small_count': 0,
            'medium_count': 0,
            'large_count': 0,
            'hydrophobic_ratio': 0,
            'charged_ratio': 0
        }
    
    pocket_aa = [three_to_one(res.resname) for res in neighbors if is_aa(res)]
    total_count = len(pocket_aa)
    
    if total_count == 0:
        return {
            'hydrophobic_count': 0,
            'polar_count': 0,
            'positive_count': 0,
            'negative_count': 0,
            'net_charge': 0,
            'small_count': 0,
            'medium_count': 0,
            'large_count': 0,
            'hydrophobic_ratio': 0,
            'charged_ratio': 0
        }
    
    hydrophobic_count = sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['hydrophobic'])
    polar_count = sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['polar'])
    positive_count = sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['charged'] and aa_properties[aa]['charge'] > 0)
    negative_count = sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['charged'] and aa_properties[aa]['charge'] < 0)
    
    properties = {
        'hydrophobic_count': hydrophobic_count,
        'polar_count': polar_count,
        'positive_count': positive_count,
        'negative_count': negative_count,
        'net_charge': sum(aa_properties[aa]['charge'] for aa in pocket_aa if aa in aa_properties),
        'small_count': sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['size'] == 'small'),
        'medium_count': sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['size'] == 'medium'),
        'large_count': sum(1 for aa in pocket_aa if aa in aa_properties and aa_properties[aa]['size'] == 'large'),
        'hydrophobic_ratio': hydrophobic_count / total_count if total_count > 0 else 0,
        'charged_ratio': (positive_count + negative_count) / total_count if total_count > 0 else 0
    }
    
    return properties


def identify_domain_and_motifs(sequence, uniprot_id, res_pos):
    """Identify known domains and motifs around the STY site"""
    # Pre-compile regex patterns for better performance
    motifs = {
        # Generic kinase motifs
        'PKA_consensus': re.compile(r'R[RK].S'),
        'PKC_consensus': re.compile(r'[ST].[RK]'),
        'CK2_consensus': re.compile(r'S..E'),
        'CK1_consensus': re.compile(r'S..[ST]'),
        'proline_directed': re.compile(r'[ST]P'),
        'GSK3_consensus': re.compile(r'S...S'),
        'acidic_directed': re.compile(r'S.[DE]'),
        'basic_directed': re.compile(r'S.[RK]'),
        'CDK_consensus': re.compile(r'[ST]P.[RK]'),
        'MAPK_consensus': re.compile(r'P.[ST]P'),
        'ATM_ATR_consensus': re.compile(r'[ST]Q'),
    }
    
    found_motifs = []
    window_size = 10
    start = max(0, res_pos - window_size)
    end = min(len(sequence), res_pos + window_size + 1)
    
    window_seq = sequence[start:end]
    centered_pos = min(res_pos, window_size)
    
    for motif_name, pattern in motifs.items():
        matches = pattern.finditer(window_seq)
        for match in matches:
            # Check if the match includes our STY site
            match_start = match.start()
            match_end = match.end()
            if match_start <= centered_pos < match_end:
                found_motifs.append(motif_name)
    
    return ', '.join(found_motifs) if found_motifs else 'None'


def find_hydrogen_bonds(residue, neighbors, max_distance=3.5):
    """Find potential hydrogen bonds involving the STY hydroxyl group"""
    # Get the hydroxyl oxygen atom
    hydroxyl_atom = None
    if residue.resname == 'SER' and 'OG' in residue:
        hydroxyl_atom = residue['OG']
    elif residue.resname == 'THR' and 'OG1' in residue:
        hydroxyl_atom = residue['OG1']
    elif residue.resname == 'TYR' and 'OH' in residue:
        hydroxyl_atom = residue['OH']
    
    if not hydroxyl_atom:
        return 0
    
    # Look for potential hydrogen bond acceptors/donors
    h_bonds = 0
    max_distance_squared = max_distance * max_distance
    
    for neighbor_res in neighbors:
        # Skip the residue itself
        if neighbor_res is residue:
            continue
            
        for atom in neighbor_res:
            # Only consider oxygen and nitrogen atoms as potential partners
            if atom.element not in ['O', 'N']:
                continue
                
            # Calculate distance squared (faster than taking square root)
            dist_squared = sum((hydroxyl_atom.coord - atom.coord) ** 2)
            if dist_squared <= max_distance_squared:
                h_bonds += 1
    
    return h_bonds


def calculate_b_factor_stats(residue, neighbors):
    """Calculate B-factor statistics for the site and surroundings"""
    # Get B-factors for the STY residue atoms
    sty_b_factors = [atom.bfactor for atom in residue]
    sty_avg_b = np.mean(sty_b_factors) if sty_b_factors else None
    
    # Get B-factors for the neighbor residues
    neighbor_b_factors = []
    for neighbor in neighbors:
        if neighbor is not residue:  # Skip the STY residue itself
            neighbor_b_factors.extend([atom.bfactor for atom in neighbor])
    
    neighbor_avg_b = np.mean(neighbor_b_factors) if neighbor_b_factors else None
    b_factor_ratio = sty_avg_b / neighbor_avg_b if (sty_avg_b is not None and neighbor_avg_b is not None and neighbor_avg_b > 0) else None
    
    return sty_avg_b, neighbor_avg_b, b_factor_ratio


def extract_secondary_structure(structure_file):
    """Extract secondary structure information from a CIF file with optimized performance"""
    ss_dict = {}  # Map of residue_id -> secondary structure type
    
    # Early return for non-CIF files (unchanged, this is already efficient)
    if not (structure_file.endswith('.cif') or structure_file.endswith('.cif.gz')):
        return ss_dict
    
    try:
        # Handle gzipped files (no change needed here, this is standard practice)
        if structure_file.endswith('.gz'):
            with gzip.open(structure_file, 'rt') as f:
                mmcif_dict = MMCIF2Dict.MMCIF2Dict(f)
        else:
            mmcif_dict = MMCIF2Dict.MMCIF2Dict(structure_file)
        
        # Check if struct_conf is in the CIF file
        if '_struct_conf.conf_type_id' not in mmcif_dict:
            return ss_dict  # Early return if no secondary structure data
        
        # Get the secondary structure data
        # Use dict.get() with default empty list to avoid potential KeyErrors
        begin_chain_ids = mmcif_dict.get('_struct_conf.beg_auth_asym_id', [])
        begin_res_ids = mmcif_dict.get('_struct_conf.beg_auth_seq_id', [])
        end_chain_ids = mmcif_dict.get('_struct_conf.end_auth_asym_id', [])
        end_res_ids = mmcif_dict.get('_struct_conf.end_auth_seq_id', [])
        ss_types = mmcif_dict.get('_struct_conf.conf_type_id', [])
        
        # Skip if any required data is missing
        if not (begin_chain_ids and begin_res_ids and end_chain_ids and end_res_ids and ss_types):
            return ss_dict
        
        # Pre-calculate length for performance
        num_entries = len(ss_types)
        
        # Pre-convert all residue IDs to integers at once rather than repeatedly
        try:
            begin_res_ids = [int(x) for x in begin_res_ids]
            end_res_ids = [int(x) for x in end_res_ids]
        except ValueError:
            # Handle potential non-integer values
            begin_res_ids_clean = []
            end_res_ids_clean = []
            for i in range(num_entries):
                try:
                    begin_res_ids_clean.append(int(begin_res_ids[i]))
                    end_res_ids_clean.append(int(end_res_ids[i]))
                except ValueError:
                    # Skip entries with non-integer residue IDs
                    continue
            begin_res_ids = begin_res_ids_clean
            end_res_ids = end_res_ids_clean
            
            # Update num_entries in case we filtered out some values
            num_entries = min(len(begin_res_ids), len(end_res_ids), len(ss_types))
        
        # Map each residue to its secondary structure
        # Use zip to iterate over all lists simultaneously
        for i, (chain_id, start_res, end_res, ss_type) in enumerate(
                zip(begin_chain_ids[:num_entries], begin_res_ids[:num_entries], 
                    end_res_ids[:num_entries], ss_types[:num_entries])):
            
            # Add each residue in the range to the dictionary
            # Use a dictionary comprehension for better performance
            ss_dict.update({(chain_id, res_id): ss_type for res_id in range(start_res, end_res + 1)})
    
    except Exception as e:
        print(f"Error extracting secondary structure: {e}")
    
    return ss_dict


def calculate_sasa_vectorized(residue, neighbors, probe_radius=1.4, n_points=100):
    """Calculate approximate SASA using a vectorized ray-casting approach
    
    Args:
        residue: The residue to calculate SASA for
        neighbors: List of neighboring residues
        probe_radius: Radius of water probe (default 1.4 Å)
        n_points: Number of points to use on unit sphere (higher = more accurate, slower)
        
    Returns:
        total_sasa: Approximate SASA in Å²
        hydroxyl_sasa: SASA specific to hydroxyl oxygen (for STY residues)
    """
    # Standard VDW radii
    vdw_radii = {
        'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8,
        'P': 1.8, 'H': 1.2, 'F': 1.47, 'CL': 1.75,
        'BR': 1.85, 'I': 1.98
    }
    default_radius = 1.7  # Default for unknown atoms
    
    # Generate points on unit sphere (evenly distributed)
    # Using Fibonacci sphere algorithm for uniform distribution
    phi = (1 + np.sqrt(5)) / 2  # Golden ratio
    indices = np.arange(n_points)
    y = 1 - (indices / (n_points - 1)) * 2  # y goes from 1 to -1
    radius = np.sqrt(1 - y * y)  # Radius at each height
    
    theta = 2 * np.pi * indices / phi  # Golden angle increment
    
    x = np.cos(theta) * radius
    z = np.sin(theta) * radius
    
    # Combine into unit vectors
    unit_sphere = np.column_stack((x, y, z))
    
    # Get hydroxyl oxygen for STY residues
    hydroxyl_atom = None
    if residue.resname == 'SER' and 'OG' in residue:
        hydroxyl_atom = residue['OG']
    elif residue.resname == 'THR' and 'OG1' in residue:
        hydroxyl_atom = residue['OG1']
    elif residue.resname == 'TYR' and 'OH' in residue:
        hydroxyl_atom = residue['OH']
    
    # Extract neighbor atom positions and radii
    neighbor_positions = []
    neighbor_radii = []
    
    # Add atoms from neighboring residues
    for res in neighbors:
        for atom in res:
            # Skip if this is from the central residue
            if atom.get_parent() is residue:
                continue
                
            neighbor_positions.append(atom.coord)
            element = atom.element.upper()
            neighbor_radii.append(vdw_radii.get(element, default_radius) + probe_radius)
    
    # Convert to numpy arrays
    neighbor_positions = np.array(neighbor_positions)
    neighbor_radii = np.array(neighbor_radii)
    
    # Process each atom in the residue
    total_area = 0
    hydroxyl_area = 0
    
    for atom in residue:
        # Get VDW radius
        element = atom.element.upper()
        atom_radius = vdw_radii.get(element, default_radius)
        atom_coord = atom.coord
        
        # Expand sphere points to actual position
        expanded_radius = atom_radius + probe_radius
        test_points = unit_sphere * expanded_radius + atom_coord
        
        # Area of each point on expanded surface
        point_area = 4 * np.pi * expanded_radius**2 / n_points
        
        # Test each ray for overlap with neighbors
        if len(neighbor_positions) > 0:
            # For each test point, compute distance to all neighbor centers
            # Shape: (n_points, n_neighbors)
            diff = test_points[:, np.newaxis, :] - neighbor_positions[np.newaxis, :, :]
            sq_dists = np.sum(diff**2, axis=2)
            
            # Compare squared distances to squared neighbor radii
            # If any squared distance is less than squared radius, point is blocked
            sq_radii = neighbor_radii**2
            blocked_mask = np.any(sq_dists <= sq_radii, axis=1)
            
            # Count accessible points
            accessible_points = n_points - np.sum(blocked_mask)
            accessible_area = accessible_points * point_area
        else:
            # No neighbors, all points are accessible
            accessible_area = 4 * np.pi * expanded_radius**2
            
        total_area += accessible_area
        
        # Track hydroxyl SASA separately
        if atom is hydroxyl_atom:
            hydroxyl_area = accessible_area
    
    return total_area, hydroxyl_area


def distance_to_nearest_vectorized(residue, chain_data):
    """Calculate distance to the nearest residue of each type and property category
    
    Args:
        residue: The target residue
        chain_data: Pre-computed chain data with KD-tree
        
    Returns:
        dict: Distances to nearest residue of each type and property category
    """
    # Get CA atom coordinates for the central residue
    try:
        ca_atom = residue['CA']
    except KeyError:
        # Return default values if CA atom is missing
        return {f"dist_to_{aa}": None for aa in "ACDEFGHIKLMNPQRSTVWY"}, {
            "dist_to_hydrophobic": None,
            "dist_to_polar": None,
            "dist_to_charged": None,
            "dist_to_acidic": None,
            "dist_to_basic": None,
            "dist_to_small": None,
            "dist_to_medium": None,
            "dist_to_large": None
        }
    
    ca_coord = np.array(ca_atom.coord)
    
    # Get distances to all residues in the chain
    distances = {}
    property_distances = {
        "hydrophobic": [],
        "polar": [],
        "charged": [],
        "acidic": [],
        "basic": [],
        "small": [],
        "medium": [],
        "large": []
    }
    
    # Collect all residue coordinates by type
    residue_coords = {aa: [] for aa in "ACDEFGHIKLMNPQRSTVWY"}
    
    # Iterate through all residues in the chain
    for res in chain_data.get('all_residues', []):
        # Skip the target residue
        if res.id == residue.id and res.get_parent() == residue.get_parent():
            continue
        
        # Skip non-amino acids
        if not is_aa(res):
            continue
            
        # Get residue type
        res_type = three_to_one(res.resname)
        if res_type == 'X':
            continue
            
        # Get CA coordinates
        try:
            res_ca = res['CA']
            coords = res_ca.coord
            
            # Add to appropriate lists
            residue_coords[res_type].append(coords)
            
            # Add to property-based lists
            props = aa_properties.get(res_type, {})
            if props.get('hydrophobic', False):
                property_distances["hydrophobic"].append(coords)
            if props.get('polar', False):
                property_distances["polar"].append(coords)
            if props.get('charged', False):
                property_distances["charged"].append(coords)
                
                # Acidic (negative) or basic (positive)
                charge = props.get('charge', 0)
                if charge < 0:
                    property_distances["acidic"].append(coords)
                elif charge > 0:
                    property_distances["basic"].append(coords)
            
            # Size properties
            size = props.get('size', '')
            if size == 'small':
                property_distances["small"].append(coords)
            elif size == 'medium':
                property_distances["medium"].append(coords)
            elif size == 'large':
                property_distances["large"].append(coords)
                
        except KeyError:
            # Skip residues without CA atoms
            continue
    
    # Calculate distances for each amino acid type
    aa_distances = {}
    for aa, coords_list in residue_coords.items():
        if coords_list:
            # Convert to numpy array
            coords_array = np.array(coords_list)
            
            # Calculate Euclidean distances
            vectors = coords_array - ca_coord
            dist_sq = np.sum(vectors * vectors, axis=1)
            min_dist = np.sqrt(np.min(dist_sq)) if len(dist_sq) > 0 else None
            
            aa_distances[f"dist_to_{aa}"] = min_dist
        else:
            aa_distances[f"dist_to_{aa}"] = None
    
    # Calculate distances for each property category
    property_results = {}
    for prop, coords_list in property_distances.items():
        if coords_list:
            # Convert to numpy array
            coords_array = np.array(coords_list)
            
            # Calculate Euclidean distances
            vectors = coords_array - ca_coord
            dist_sq = np.sum(vectors * vectors, axis=1)
            min_dist = np.sqrt(np.min(dist_sq)) if len(dist_sq) > 0 else None
            
            property_results[f"dist_to_{prop}"] = min_dist
        else:
            property_results[f"dist_to_{prop}"] = None
    
    return aa_distances, property_results


def num_nearest_type_vectorized(residue, neighbors):
    """Count the number of each amino acid type and property within 10Å
    
    Args:
        residue: The target residue
        neighbors: List of neighboring residues within 10Å
        
    Returns:
        dict: Counts of each amino acid type and property category
    """
    # Initialize counts for each amino acid
    aa_counts = {f"num_{aa}_10A": 0 for aa in "ACDEFGHIKLMNPQRSTVWY"}
    
    # Initialize counts for property categories
    property_counts = {
        "num_hydrophobic_10A": 0,
        "num_polar_10A": 0,
        "num_charged_10A": 0,
        "num_acidic_10A": 0,
        "num_basic_10A": 0,
        "num_small_10A": 0,
        "num_medium_10A": 0,
        "num_large_10A": 0
    }
    
    # Count residues by type and property
    for res in neighbors:
        # Skip non-amino acids
        if not is_aa(res):
            continue
            
        # Get residue type
        res_type = three_to_one(res.resname)
        if res_type == 'X':
            continue
            
        # Increment type count
        aa_counts[f"num_{res_type}_10A"] += 1
        
        # Increment property counts
        props = aa_properties.get(res_type, {})
        if props.get('hydrophobic', False):
            property_counts["num_hydrophobic_10A"] += 1
        if props.get('polar', False):
            property_counts["num_polar_10A"] += 1
        if props.get('charged', False):
            property_counts["num_charged_10A"] += 1
            
            # Acidic (negative) or basic (positive)
            charge = props.get('charge', 0)
            if charge < 0:
                property_counts["num_acidic_10A"] += 1
            elif charge > 0:
                property_counts["num_basic_10A"] += 1
        
        # Size properties
        size = props.get('size', '')
        if size == 'small':
            property_counts["num_small_10A"] += 1
        elif size == 'medium':
            property_counts["num_medium_10A"] += 1
        elif size == 'large':
            property_counts["num_large_10A"] += 1
    
    return aa_counts, property_counts

def calculate_cavity_volume_vectorized(hydroxyl_atom, neighbors, grid_spacing=0.7, radius=6.0):
    """Calculate approximate cavity volume around STY hydroxyl group using a grid-based approach
    
    Args:
        hydroxyl_atom: The hydroxyl oxygen atom of the STY residue
        neighbors: List of neighboring residues
        grid_spacing: Grid resolution in Ångstroms (lower = more accurate but slower)
        radius: Maximum distance from hydroxyl to consider for cavity in Ångstroms
        
    Returns:
        float: Approximate cavity volume in Å³
    """
    if not hydroxyl_atom:
        return 0.0
    
    # Standard VDW radii
    vdw_radii = {
        'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8,
        'P': 1.8, 'H': 1.2, 'F': 1.47, 'CL': 1.75,
        'BR': 1.85, 'I': 1.98
    }
    default_radius = 1.7  # Default for unknown atoms
    water_radius = 1.4  # Radius of water probe
    
    # Get hydroxyl center
    hydroxyl_coord = np.array(hydroxyl_atom.coord)
    
    # Create 3D grid around hydroxyl
    grid_size = int(2 * radius / grid_spacing)
    
    # Performance optimization: Use a sparse approach instead of a full grid
    # Generate randomly distributed points within the sphere
    num_points = int((4/3) * np.pi * (radius**3) / (grid_spacing**3))
    # Limit maximum points for performance
    num_points = min(num_points, 10000)
    
    # Generate random points within a cube
    random_points = np.random.uniform(-radius, radius, (num_points, 3))
    
    # Filter to keep only points within the sphere
    distances_sq = np.sum(random_points**2, axis=1)
    mask = distances_sq <= radius**2
    grid_points = random_points[mask]
    
    # Calculate actual coordinates by adding hydroxyl position
    grid_coords = grid_points + hydroxyl_coord
    
    # Extract neighbor atom coordinates and VDW radii
    atom_coords = []
    atom_radii = []
    
    for res in neighbors:
        for atom in res:
            # Skip if this is the hydroxyl atom itself
            if atom == hydroxyl_atom:
                continue
                
            atom_coords.append(atom.coord)
            element = atom.element.upper()
            atom_radii.append(vdw_radii.get(element, default_radius) + water_radius)
    
    # Handle the case with no neighbors
    if not atom_coords:
        return (4/3) * np.pi * radius**3
    
    # Convert to numpy arrays
    atom_coords = np.array(atom_coords)
    atom_radii = np.array(atom_radii)
    
    # Calculate squared distances from each grid point to each atom center
    # Use a chunked approach for large number of atoms to avoid memory issues
    chunk_size = 10000
    is_occupied = np.zeros(len(grid_coords), dtype=bool)
    
    for i in range(0, len(grid_coords), chunk_size):
        chunk_coords = grid_coords[i:i+chunk_size]
        
        # Reshape for broadcasting: (n_grid_points, 1, 3) - (1, n_atoms, 3)
        diff = chunk_coords.reshape(-1, 1, 3) - atom_coords.reshape(1, -1, 3)
        dist_sq = np.sum(diff * diff, axis=2)  # (n_grid_points, n_atoms)
        
        # Compare squared distances to squared radii
        radii_sq = atom_radii * atom_radii  # Squared radii
        chunk_occupied = (dist_sq <= radii_sq).any(axis=1)
        
        is_occupied[i:i+len(chunk_coords)] = chunk_occupied
    
    # Grid points that are not occupied are accessible (part of cavity)
    is_accessible = ~is_occupied
    
    # Calculate cavity volume
    accessible_count = np.sum(is_accessible)
    total_count = len(grid_coords)
    
    # Calculate volume based on proportion of accessible points
    sphere_volume = (4/3) * np.pi * radius**3
    cavity_volume = (accessible_count / total_count) * sphere_volume
    
    return cavity_volume


def calculate_motif_plddt(chain_data, residue, plddt_dict, window=7):
    """Calculate mean pLDDT for the sequence motif around a residue
    
    Args:
        chain_data: Pre-computed chain data
        residue: The target residue
        plddt_dict: Dictionary mapping residue IDs to pLDDT values
        window: Size of the sequence window on each side (default: 7)
        
    Returns:
        float: Mean pLDDT of residues in the motif
    """
    res_id = residue.id[1]
    
    # Get the sequence IDs from chain data
    seq_ids = chain_data['seq_ids']
    
    # Find the position in the sequence
    if res_id not in seq_ids:
        return None
    
    seq_pos = seq_ids.index(res_id)
    
    # Get residue IDs in the motif window
    start = max(0, seq_pos - window)
    end = min(len(seq_ids), seq_pos + window + 1)
    
    motif_res_ids = seq_ids[start:end]
    
    # Get pLDDT values for each residue in the motif
    motif_plddt = [plddt_dict.get(res_id, 0) for res_id in motif_res_ids]
    
    # Remove missing values (0)
    motif_plddt = [p for p in motif_plddt if p > 0]
    
    if not motif_plddt:
        return None
    
    return np.mean(motif_plddt)


def calculate_neighbor_plddt(neighbors, plddt_dict):
    """Calculate mean pLDDT for the neighboring residues
    
    Args:
        neighbors: List of neighboring residues
        plddt_dict: Dictionary mapping residue IDs to pLDDT values
        
    Returns:
        float: Mean pLDDT of neighboring residues
    """
    if not neighbors:
        return None
    
    # Get pLDDT values for each neighbor
    neighbor_plddt = []
    
    for res in neighbors:
        # Get residue ID
        res_id = res.id[1]
        
        # Get pLDDT value
        plddt = plddt_dict.get(res_id, 0)
        
        if plddt > 0:
            neighbor_plddt.append(plddt)
    
    if not neighbor_plddt:
        return None
    
    return np.mean(neighbor_plddt)


def cache_turns_and_helices_fast(model, ss_dict):
    """Cache only turn and helix residues in the model for faster lookup"""
    if not ss_dict:
        return {'turn': {'residues': [], 'ca_coords': {}}, 
                'helix': {'residues': [], 'ca_coords': {}}}

    # Define secondary structure keywords (exact match instead of substring)
    turn_keywords = {
        'TURN_TY1_P', 'TURN_TY2_P', 'TURN_TY1_PM', 'TURN_TY2_PM', 'TURN_TY3_P', 'BEND'
    }
    helix_keywords = {
        'HELX_LH_PP_P', 'HELX_RH_PP_P', 'HELX_RH_3T_P', 'HELX_RH_PI_P', 'HELX_LH_P'
    }

    turn_data = {'residues': [], 'ca_coords': {}}
    helix_data = {'residues': [], 'ca_coords': {}}

    for (chain_id, res_id), ss_type in ss_dict.items():
        try:
            res = model[chain_id][(' ', res_id, ' ')]
            if 'CA' not in res:
                continue
            coord = res['CA'].coord

            if ss_type in turn_keywords:
                turn_data['residues'].append((chain_id, res_id))
                turn_data['ca_coords'][(chain_id, res_id)] = coord
            elif ss_type in helix_keywords:
                helix_data['residues'].append((chain_id, res_id))
                helix_data['ca_coords'][(chain_id, res_id)] = coord
        except KeyError:
            continue

    return {
        'turn': turn_data,
        'helix': helix_data
    }


def calculate_turn_helix_distances_fast(residue, seq_ids, ss_data):
    """Vectorized distance calculation to nearest turn and helix"""
    chain_id = residue.get_parent().id
    res_id = residue.id[1]

    try:
        ca_coord = np.array(residue['CA'].coord)
    except KeyError:
        return None, None, None, None

    target_seq_pos = seq_ids.index(res_id) if res_id in seq_ids else None

    def get_min_distances(ss_type):
        if ss_type not in ss_data:
            return None, None

        ss_residues = ss_data[ss_type]['residues']
        ss_coords = ss_data[ss_type]['ca_coords']

        # Sequence distance (same chain only)
        if target_seq_pos is not None:
            same_chain_ids = [r for (c, r) in ss_residues if c == chain_id and r in seq_ids]
            if same_chain_ids:
                ss_seq_pos = np.array([seq_ids.index(r) for r in same_chain_ids])
                seq_dists = np.abs(ss_seq_pos - target_seq_pos)
                min_seq = int(np.min(seq_dists))
            else:
                min_seq = None
        else:
            min_seq = None

        # Spatial distance (across all chains)
        all_coords = [coord for key, coord in ss_coords.items() if key != (chain_id, res_id)]
        if all_coords:
            all_coords = np.array(all_coords)
            dists = np.sqrt(np.sum((all_coords - ca_coord) ** 2, axis=1))
            min_spatial = float(np.min(dists))
        else:
            min_spatial = None

        return min_seq, min_spatial

    seq_turn, spatial_turn = get_min_distances('turn')
    seq_helix, spatial_helix = get_min_distances('helix')

    return seq_turn, spatial_turn, seq_helix, spatial_helix



from Bio.PDB.HSExposure import HSExposureCA, HSExposureCB, ExposureCN

from Bio.PDB.HSExposure import HSExposureCA, HSExposureCB, ExposureCN


def calculate_hse_for_residue(model, target_residue, radius=13.0):
    """
    Calculate HSE and coordination number for a specific residue using BioPython.

    Args:
        model: Bio.PDB Model object
        target_residue: Residue object to analyze
        radius: Sphere radius in Å (default: 13.0)

    Returns:
        dict: Dictionary of HSE values and coordination number for the residue
    """
    try:
        # These update all residues in the model with xtra annotations
        HSExposureCA(model, radius)
        HSExposureCB(model, radius)
        ExposureCN(model, radius)

        # Extract the info for the specific residue
        res_id = target_residue.id[1]
        result = {}

        # CA-based HSE
        if "EXP_HSE_A_U" in target_residue.xtra:
            result["HSE_CA_U"] = target_residue.xtra["EXP_HSE_A_U"]
            result["HSE_CA_D"] = target_residue.xtra["EXP_HSE_A_D"]
            u = result["HSE_CA_U"]
            d = result["HSE_CA_D"]
            result["HSE_CA_RATIO"] = u / d if d > 0 else None

        # CB-based HSE
        if "EXP_HSE_B_U" in target_residue.xtra:
            result["HSE_CB_U"] = target_residue.xtra["EXP_HSE_B_U"]
            result["HSE_CB_D"] = target_residue.xtra["EXP_HSE_B_D"]
            u = result["HSE_CB_U"]
            d = result["HSE_CB_D"]
            result["HSE_CB_RATIO"] = u / d if d > 0 else None

        # Coordination number
        if "EXP_CN" in target_residue.xtra:
            result["EXP_CN"] = target_residue.xtra["EXP_CN"]
            result["EXP_CN_COUNT"] = target_residue.xtra.get("EXP_CN_COUNT", None)
        #print("RESULT")
        #print(result)
        return result

    except Exception as e:
        print(f"Error calculating HSE for residue: {e}")
        return {}
    



def calculate_hse_for_residues(model, residues, radius=10.0):
    """
    Calculate HSE efficiently for a list of residues using pre-computation.
    
    Args:
        model: Bio.PDB Model object
        residues: List of Bio.PDB Residue objects
        radius: Sphere radius in Å (default: 13.0)
        
    Returns:
        dict: Mapping from residue IDs to HSE dictionaries
    """
    # Get residue IDs for faster lookup later
    residue_ids = {residue.id[1]: residue for residue in residues}
    
    try:
        # Pre-calculate all HSE values at once for the model
        # This is the costly part, but we only do it once
        hse_ca = HSExposureCA(model, radius)
        hse_cb = HSExposureCB(model, radius)
        
        # Create results dictionary directly using dictionary comprehension
        # This avoids the overhead of repeatedly building dictionaries in a loop
        results = {}
        
        # Process all residues at once
        for res_id, residue in residue_ids.items():
            hse = {}
            
            # HSE-CA (most important values)
            if "EXP_HSE_A_U" in residue.xtra:
                u = residue.xtra["EXP_HSE_A_U"]
                d = residue.xtra["EXP_HSE_A_D"]
                hse["HSE_CA_U"] = u
                hse["HSE_CA_D"] = d
                hse["HSE_CA_RATIO"] = u / d if d > 0 else None
            
            # HSE-CB
            if "EXP_HSE_B_U" in residue.xtra:
                u = residue.xtra["EXP_HSE_B_U"]
                d = residue.xtra["EXP_HSE_B_D"]
                hse["HSE_CB_U"] = u
                hse["HSE_CB_D"] = d
                hse["HSE_CB_RATIO"] = u / d if d > 0 else None
            
            results[res_id] = hse
            
        return results
    
    except Exception as e:
        print(f"Warning: HSE calculation failed: {e}")
        # Return empty results rather than failing completely
        return {res_id: {} for res_id in residue_ids}
    

import numpy as np

def calculate_hse_fast(model, residues, chain_data, radius=13.0):
    """
    Fast HSE calculation that only processes specific residues using the existing
    find_neighbors_fast function for better performance.
    
    Args:
        model: BioPython Model object
        residues: List of Bio.PDB Residue objects to calculate HSE for
        chain_data: Pre-computed chain data with KD-tree
        radius: Sphere radius in Å (default: 13.0)
        
    Returns:
        dict: Mapping from residue IDs to HSE dictionaries
    """
    results = {}
    
    # Process each requested residue
    for residue in residues:
        res_id = residue.id[1]
        
        # Skip non-amino acids
        if not is_aa(residue):
            results[res_id] = {}
            continue
            
        # Get CA and CB atoms (if they exist)
        ca_atom = residue['CA'] if 'CA' in residue else None
        cb_atom = residue['CB'] if 'CB' in residue else None
        
        # Skip residues without required atoms
        if not ca_atom:
            results[res_id] = {}
            continue
            
        # Find neighbors using our fast function
        neighbors = find_neighbors_fast(chain_data, residue, radius)
        
        # Initialize HSE data
        hse = {}
        
        # Define up and down regions based on CA-CB vector (or pseudo-CB)
        if cb_atom:
            # Use actual CB atom if available
            ca_cb_vector = cb_atom.coord - ca_atom.coord
        else:
            # Calculate pseudo-CB vector based on backbone
            n_atom = residue['N'] if 'N' in residue else None
            c_atom = residue['C'] if 'C' in residue else None
            
            if n_atom and c_atom:
                # Pseudo-CB direction for Glycine (perpendicular to peptide plane)
                ca_n_vector = n_atom.coord - ca_atom.coord
                ca_c_vector = c_atom.coord - ca_atom.coord
                ca_cb_vector = np.cross(ca_n_vector, ca_c_vector)
                # Normalize vector
                norm = np.sqrt(np.sum(ca_cb_vector * ca_cb_vector))
                if norm > 0:
                    ca_cb_vector = ca_cb_vector / norm
                else:
                    ca_cb_vector = np.array([0, 0, 0])
            else:
                ca_cb_vector = np.array([0, 0, 0])
        
        # Count neighbors in up and down half-spheres for CA-based HSE
        up_count = 0
        down_count = 0
        
        for neighbor_res in neighbors:
            # Skip if it's the same residue
            if neighbor_res is residue:
                continue
                
            # Count contributions from each atom in the neighbor residue
            for atom in neighbor_res:
                # Calculate vector from CA to neighbor atom
                ca_neighbor_vector = atom.coord - ca_atom.coord
                
                # Check if this atom is within the radius
                dist_sq = np.sum(ca_neighbor_vector * ca_neighbor_vector)
                if dist_sq <= radius * radius:
                    # Dot product determines if neighbor is in up or down half-sphere
                    if np.dot(ca_cb_vector, ca_neighbor_vector) >= 0:
                        up_count += 1
                    else:
                        down_count += 1
        
        # Store CA-based HSE results
        hse["HSE_CA_U"] = up_count
        hse["HSE_CA_D"] = down_count
        hse["HSE_CA_RATIO"] = up_count / down_count if down_count > 0 else None
        
        # Calculate CB-based HSE if CB exists
        if cb_atom:
            cb_up_count = 0
            cb_down_count = 0
            
            for neighbor_res in neighbors:
                # Skip if it's the same residue
                if neighbor_res is residue:
                    continue
                    
                # Count contributions from each atom in the neighbor residue
                for atom in neighbor_res:
                    # Calculate vector from CB to neighbor atom
                    cb_neighbor_vector = atom.coord - cb_atom.coord
                    
                    # Check if this atom is within the radius
                    dist_sq = np.sum(cb_neighbor_vector * cb_neighbor_vector)
                    if dist_sq <= radius * radius:
                        # Dot product determines if neighbor is in up or down half-sphere
                        if np.dot(ca_cb_vector, cb_neighbor_vector) >= 0:
                            cb_up_count += 1
                        else:
                            cb_down_count += 1
            
            # Store CB-based HSE results
            hse["HSE_CB_U"] = cb_up_count
            hse["HSE_CB_D"] = cb_down_count
            hse["HSE_CB_RATIO"] = cb_up_count / cb_down_count if cb_down_count > 0 else None
        
        # Store results for this residue
        results[res_id] = hse
    
    return results

def cache_all_secondary_structures(model, ss_dict):
    """Cache all secondary structure residues in the model for faster lookup
    
    Args:
        model: BioPython Model object
        ss_dict: Dictionary mapping (chain_id, res_id) to secondary structure
        
    Returns:
        dict: Dictionary with keys for each secondary structure type, each containing:
            - 'residues': List of (chain_id, res_id) tuples
            - 'ca_coords': Dictionary mapping (chain_id, res_id) to CA coordinates
    """
    if not ss_dict:
        return {}
    
    # Define secondary structure categories
    ss_categories = {
        'helix': ['HELX_LH_PP_P', 'HELX_RH_PP_P', 'HELX_RH_3T_P', 'HELX_RH_PI_P', 'HELX_LH_P'],
        'strand': ['STRN', 'SHEET'],
        'turn': ['TURN_TY1_P', 'TURN_TY2_P', 'TURN_TY1_PM', 'TURN_TY2_PM', 'TURN_TY3_P', 'BEND'],
        'coil': ['COIL']
    }
    
    # Initialize result dictionary
    ss_data = {ss_type: {'residues': [], 'ca_coords': {}} for ss_type in ss_categories}
    
    # Create a mapping from specific SS codes to category
    ss_type_to_category = {}
    for category, ss_types in ss_categories.items():
        for ss_type in ss_types:
            ss_type_to_category[ss_type] = category
    
    # Process all residues in the SS dictionary
    for (c_id, r_id), ss_type in ss_dict.items():
        # Determine the secondary structure category
        category = None
        for cat, patterns in ss_categories.items():
            if any(pattern in ss_type for pattern in patterns):
                category = cat
                break
        
        # If no category matched, use the raw SS type
        if category is None:
            category = 'other'
            if 'other' not in ss_data:
                ss_data['other'] = {'residues': [], 'ca_coords': {}}
        
        # Add the residue to the appropriate category
        ss_data[category]['residues'].append((c_id, r_id))
        
        # Try to get CA coordinates
        if c_id in model:
            chain = model[c_id]
            try:
                # BioPython residue selector
                res = chain[(' ', r_id, ' ')]
                if 'CA' in res:
                    ss_data[category]['ca_coords'][(c_id, r_id)] = res['CA'].coord
            except (KeyError, Exception):
                continue
    
    return ss_data


def calculate_secondary_structure_distances(residue, seq_ids, ss_data):
    """Vectorized calculation of distances to nearest secondary structures of each type
    
    Args:
        residue: The target residue
        seq_ids: List of residue IDs in sequence order
        ss_data: Dictionary of secondary structure data from cache_all_secondary_structures
        
    Returns:
        dict: Dictionary mapping secondary structure types to tuples of 
              (sequence distance, spatial distance in Ångstroms)
    """
    # Get residue info
    chain_id = residue.get_parent().id
    res_id = residue.id[1]
    
    # Get CA atom of the residue
    try:
        ca_atom = residue['CA']
        ca_coord = np.array(ca_atom.coord)
    except KeyError:
        return {ss_type: (None, None) for ss_type in ss_data}
    
    # Get target sequence position
    target_seq_pos = None
    if res_id in seq_ids:
        target_seq_pos = seq_ids.index(res_id)
    
    # Initialize results dictionary
    distances = {}
    
    # Calculate distances for each secondary structure type
    for ss_type, data in ss_data.items():
        # Get residues and coordinates for this SS type
        ss_residues = data['residues']
        ss_ca_coords = data['ca_coords']
        
        if not ss_residues:
            distances[ss_type] = (None, None)
            continue
        
        # Initialize distances
        min_seq_distance = None
        min_spatial_distance = None
        
        # Calculate sequence distance (if in same chain)
        if target_seq_pos is not None:
            # Find SS elements in the same chain
            same_chain_ss = [(c, r) for (c, r) in ss_residues if c == chain_id]
            
            if same_chain_ss:
                # Extract residue IDs and find those that are in seq_ids
                ss_res_ids = [r for (_, r) in same_chain_ss]
                valid_ss_res_ids = [r for r in ss_res_ids if r in seq_ids]
                
                if valid_ss_res_ids:
                    # Calculate sequence distances vectorized
                    ss_seq_pos = [seq_ids.index(r) for r in valid_ss_res_ids]
                    seq_distances = np.abs(np.array(ss_seq_pos) - target_seq_pos)
                    min_seq_distance = np.min(seq_distances)
        
        # Calculate spatial distance
        all_ss_coords = []
        
        for key in ss_ca_coords:
            if key != (chain_id, res_id):  # Skip the residue itself
                all_ss_coords.append(ss_ca_coords[key])
        
        if all_ss_coords:
            # Convert to numpy array for vectorized operations
            all_ss_coords = np.array(all_ss_coords)
            
            # Calculate distances vectorized
            vectors = all_ss_coords - ca_coord
            distances_sq = np.sum(vectors * vectors, axis=1)
            min_spatial_distance = np.sqrt(np.min(distances_sq))
        
        distances[ss_type] = (min_seq_distance, min_spatial_distance)
    
    return distances

def process_structure_file(structure_file, total_files, file_index):
    """Process a single structure file (PDB or CIF, potentially gzipped) and extract features for STY residues"""
    results = []
    start_time = time.time()
    
    # Extract UniProt ID from filename (assuming format like AF-P12345-F1-model_v4.pdb.gz or .cif.gz)
    filename = os.path.basename(structure_file)
    if 'AF-' in filename and '-F' in filename:
        uniprot_id = filename.split('AF-')[1].split('-F')[0]
    else:
        uniprot_id = filename.split('.')[0]  # Fallback
    
    print(f"[{file_index}/{total_files}] Processing {filename}...")
    
    # Determine file type and select appropriate parser
    if filename.endswith('.cif.gz') or filename.endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
        is_cif = True
    else:
        parser = PDBParser(QUIET=True)
        is_cif = False
    
    total_sty_found = 0
    
    try:
        # Extract secondary structure information if it's a CIF file
        ss_dict = {}
        if is_cif:
            ss_dict = extract_secondary_structure(structure_file)
        
        # Handle gzipped files
        if structure_file.endswith('.gz'):
            with gzip.open(structure_file, 'rt') as f:
                temp_file = f"{structure_file}.temp"
                with open(temp_file, 'w') as temp:
                    temp.write(f.read())
                structure = parser.get_structure(uniprot_id, temp_file)
                os.remove(temp_file)  # Clean up
        else:
            structure = parser.get_structure(uniprot_id, structure_file)
        
        model = structure[0]
        
        # Get pLDDT scores - only do this once per file
        plddt_dict = parse_plddt_from_structure_file(structure_file)
        
        # Cache all turn residues for faster lookup
        ss_data = cache_all_secondary_structures(model, ss_dict)
        
        
        
        # Process each chain
        for chain_idx, chain in enumerate(model):
            chain_id = chain.id
            print(f"  Processing chain {chain_id} ({chain_idx+1}/{len(model)} chains)...")
            
            # Build chain data once - contains sequence, KD-tree, etc.
            chain_data = build_chain_data(chain)
            
            # Add all residues to chain_data for distance calculations
            chain_data['all_residues'] = list(chain.get_residues())
            
            # Get STY residues from pre-computed data
            sty_residues = chain_data['sty_residues']
            
            if not sty_residues:
                print(f"  No S/T/Y residues found in chain {chain_id}")
                continue
                
            print(f"  Found {len(sty_residues)} S/T/Y residues in chain {chain_id}")
            
            # Get HSE data for this chain
            hse_total_result = calculate_hse_fast(model,sty_residues, chain_data)
            #print("HSE_TOTAL_RESULT")
            #print(hse_total_result)
            # Process each STY residue
            for i, residue in enumerate(sty_residues):
                if (i + 1) % 100 == 0:
                    print(f"    Processed {i+1}/{len(sty_residues)} STY residues in chain {chain_id} (current: {residue.id[1]})")
                #print(residue)
                res_id = residue.id[1]
                res_code = three_to_one(residue.resname)
                
                
                # Get pLDDT score (faster lookup)
                plddt = plddt_dict.get(res_id, None)
                
                # Find neighbors using KD-tree (much faster)
                neighbors = find_neighbors_fast(chain_data, residue)
                
                # Get sequence motif
                motif = get_sequence_motif(chain, residue)
                
                # Analyze pocket properties
                pocket_props = analyze_pocket_properties(neighbors)
                
                # Get the hydroxyl oxygen atom
                hydroxyl_atom = None
                if residue.resname == 'SER' and 'OG' in residue:
                    hydroxyl_atom = residue['OG']
                elif residue.resname == 'THR' and 'OG1' in residue:
                    hydroxyl_atom = residue['OG1']
                elif residue.resname == 'TYR' and 'OH' in residue:
                    hydroxyl_atom = residue['OH']
                
                # Calculate hydroxyl exposure using vectorized function
                hydroxyl_exposure, backbone_contacts = calculate_hydroxyl_exposure_vectorized(
                    hydroxyl_atom, neighbors)
                
                # Calculate distance to nearest residue types (vectorized)
                nearest_aa_distances, nearest_property_distances = distance_to_nearest_vectorized(
                    residue, chain_data)
                
                # Count nearby residues by type (vectorized)
                aa_counts, property_counts = num_nearest_type_vectorized(
                    residue, neighbors)
                
                # Find hydrogen bonds using vectorized function
                h_bonds = find_hydrogen_bonds_vectorized(hydroxyl_atom, neighbors)
                
                # Get secondary structure (if available)
                ss_key = (chain_id, res_id)
                secondary_structure = ss_dict.get(ss_key, 'Unknown')
                
                # Map secondary structure codes to more descriptive names
                ss_descriptions = {
                    'HELX_LH_PP_P': 'Alpha Helix',
                    'HELX_RH_PP_P': 'Right-handed Helix',
                    'HELX_RH_3T_P': '3/10 Helix',
                    'HELX_RH_PI_P': 'Pi Helix',
                    'HELX_LH_P': 'Left-handed Helix',
                    'STRN': 'Beta Strand',
                    'TURN_TY1_P': 'Type I Turn',
                    'TURN_TY2_P': 'Type II Turn',
                    'TURN_TY1_PM': 'Type I Prime Turn',
                    'TURN_TY2_PM': 'Type II Prime Turn',
                    'TURN_TY3_P': 'Type III Turn',
                    'BEND': 'Bend',
                    'SHEET': 'Beta Sheet',
                    'COIL': 'Coil'
                }
                ss_description = ss_descriptions.get(secondary_structure, secondary_structure)
                
                # Calculate additional metrics
                
                # 1. Cavity volume around hydroxyl
                #cavity_volume = calculate_cavity_volume_vectorized(hydroxyl_atom, neighbors)
                
                # 2. Half Sphere Exposure from pre-calculated data
                #residue_hse = chain_hse.get(res_id, {})
                
                # 3. Mean pLDDT for the motif
                motif_plddt = calculate_motif_plddt(chain_data, residue, plddt_dict)
                
                # 4. Mean pLDDT for neighbors
                #neighbor_plddt = calculate_neighbor_plddt(neighbors, plddt_dict)
                
                # 5. Distance to nearest turn
                ss_distances = calculate_secondary_structure_distances(
                    residue, chain_data['seq_ids'], ss_data)
                #print("SS__DISTANCES")
                #print(ss_distances)
                residue_hse = hse_total_result.get(res_id, {})
                #print(residue_hse)
                # Compile result
                result = {
                    'UniProtID': uniprot_id,
                    'ResidueNumber': res_id,
                    'Site': f"{uniprot_id}_{res_id}",
                    'ResidueType': res_code,
                    'Motif': motif,
                    'pLDDT': plddt,
                    'NeighborCount': len(neighbors),
                    'ChainID': chain_id,
                    'SecondaryStructure': secondary_structure,
                    'SecondaryStructureDesc': ss_description,
                    'HydroxylExposure': hydroxyl_exposure,
                    'BackboneContacts': backbone_contacts,
                    'HydrogenBonds': h_bonds,
                    
                    # New metrics
                    #'CavityVolume': cavity_volume,
                    
                    
                    
                    'MotifPLDDT': motif_plddt,
                    # Secondary structure distances
                    'SeqDistToHelix': ss_distances.get('helix', (None, None))[0],
                    'SpatialDistToHelix': ss_distances.get('helix', (None, None))[1],
                    'SeqDistToStrand': ss_distances.get('strand', (None, None))[0],
                    'SpatialDistToStrand': ss_distances.get('strand', (None, None))[1],
                    'SeqDistToTurn': ss_distances.get('turn', (None, None))[0],
                    'SpatialDistToTurn': ss_distances.get('turn', (None, None))[1],
                    'SeqDistToCoil': ss_distances.get('coil', (None, None))[0],
                    'SpatialDistToCoil': ss_distances.get('coil', (None, None))[1],
                }
                
                # Add pocket properties
                result.update(pocket_props)
                
                # Add nearest residue distances
                result.update(nearest_aa_distances)
                result.update(nearest_property_distances)
                
                # Add residue counts
                result.update(aa_counts)
                result.update(property_counts)
                
                result.update({
                    'HSE_CA_U': residue_hse.get('HSE_CA_U'),
                    'HSE_CA_D': residue_hse.get('HSE_CA_D'),
                    'HSE_CA_RATIO': residue_hse.get('HSE_CA_RATIO'),
                    'HSE_CB_U': residue_hse.get('HSE_CB_U'),
                    'HSE_CB_D': residue_hse.get('HSE_CB_D'),
                    'HSE_CB_RATIO': residue_hse.get('HSE_CB_RATIO'),
                    #'EXP_CN': residue_hse.get('EXP_CN'),
                    #'EXP_CN_COUNT': residue_hse.get('EXP_CN_COUNT'),
                })

                result.update(residue_hse)

                results.append(result)
                total_sty_found += 1
    
    except Exception as e:
        print(f"Error processing {structure_file}: {e}")
    
    elapsed_time = time.time() - start_time
    print(f"✓ [{file_index}/{total_files}] Completed {filename} in {elapsed_time:.1f}s, found {total_sty_found} STY sites")
    
    return results


def main():
    parser = argparse.ArgumentParser(description='Analyze STY residues in AlphaFold structure files')
    parser.add_argument('structure_dir', nargs='?', default='F:/Kinome/Complete_AF_Proteome',
                        help='Directory containing AlphaFold structure files (.pdb, .cif, or gzipped versions)')
    parser.add_argument('--output', '-o', default='sty_analysis.csv', help='Output CSV file')
    parser.add_argument('--temp-dir', '-d', default='./temp_results', help='Directory to store temporary result files')
    parser.add_argument('--file-types', '-t', default='cif', help='Comma-separated list of file types to process (default: cif)')
    parser.add_argument('--max-files', '-m', type=int, default=0, help='Maximum number of files to process (0 for all)')
    parser.add_argument('--summary', '-s', action='store_true', help='Generate summary statistics after analysis')
    parser.add_argument('--fallback-pdb', '-f', action='store_true', help='Fallback to PDB files if CIF not available')
    args = parser.parse_args()
    
    # Create temp directory if it doesn't exist
    if not os.path.exists(args.temp_dir):
        print(f"Creating temporary directory: {args.temp_dir}")
        os.makedirs(args.temp_dir)
    
    # Parse file types to look for
    file_extensions = args.file_types.split(',')
    
    # Get list of structure files, prioritizing CIF over PDB for the same UniProt ID
    structure_files = []
    uniprot_files = {}  # Maps UniProt ID to file path
    
    # First pass: collect all files by UniProt ID
    for f in os.listdir(args.structure_dir):
        # Extract UniProt ID from filename
        if 'AF-' in f and '-F' in f:
            uniprot_id = f.split('AF-')[1].split('-F')[0]
            file_path = os.path.join(args.structure_dir, f)
            
            # Check if this is a supported file type
            is_supported = any(f.endswith(ext) or f.endswith(f"{ext}.gz") for ext in file_extensions)
            
            # If fallback is enabled, also check for PDB files
            if args.fallback_pdb and not is_supported:
                is_supported = any(f.endswith(ext) or f.endswith(f"{ext}.gz") for ext in ['pdb'])
            
            if is_supported:
                if uniprot_id not in uniprot_files:
                    uniprot_files[uniprot_id] = []
                uniprot_files[uniprot_id].append(file_path)
    
    # Second pass: prioritize CIF files over PDB
    for uniprot_id, files in uniprot_files.items():
        # Sort files by extension preference
        sorted_files = sorted(files, 
                            key=lambda x: 0 if x.endswith('.cif') or x.endswith('.cif.gz') else 1)
        
        # Add the highest priority file for this UniProt ID
        if sorted_files:
            structure_files.append(sorted_files[0])
    
    # Limit number of files if requested
    if args.max_files > 0 and args.max_files < len(structure_files):
        print(f"Limiting analysis to first {args.max_files} files out of {len(structure_files)} total")
        structure_files = structure_files[:args.max_files]
    
    total_files = len(structure_files)
    print(f"Found {total_files} unique structure files. Starting analysis...")
    
    # Track overall progress
    start_time = time.time()
    total_sty_sites = 0
    processed_file_count = 0
    
    # Process files sequentially, saving results to individual files
    for file_index, structure_file in enumerate(structure_files, 1):
        results = process_structure_file(structure_file, total_files, file_index)
        total_sty_sites += len(results)
        
        # Save this file's results to a temporary CSV
        if results:
            temp_file = os.path.join(args.temp_dir, f"temp_results_{file_index}.csv")
            temp_df = pd.DataFrame(results)
            temp_df.to_csv(temp_file, index=False)
            processed_file_count += 1
        
        # Print progress summary
        elapsed = time.time() - start_time
        files_left = total_files - file_index
        if files_left > 0 and file_index > 0:
            avg_time_per_file = elapsed / file_index
            est_time_left = avg_time_per_file * files_left
            print(f"\nPROGRESS: {file_index}/{total_files} files processed, {total_sty_sites} STY sites found")
            print(f"ESTIMATED TIME REMAINING: {est_time_left/60:.1f} minutes ({est_time_left/3600:.2f} hours)")
    
    # Combine all temp files into final result
    print(f"\nCombining {processed_file_count} temporary files into final output...")
    
    temp_files = [os.path.join(args.temp_dir, f) for f in os.listdir(args.temp_dir) if f.startswith("temp_results_") and f.endswith(".csv")]
    
    if temp_files:
        # Read and combine all temporary CSV files
        dfs = []
        for temp_file in temp_files:
            try:
                df = pd.read_csv(temp_file)
                dfs.append(df)
            except Exception as e:
                print(f"Error reading {temp_file}: {e}")
        
        if dfs:
            # Combine all dataframes
            combined_df = pd.concat(dfs, ignore_index=True)
            combined_df.to_csv(args.output, index=False)
            
            total_time = time.time() - start_time
            print(f"\nAnalysis complete. Combined {len(combined_df)} STY sites from {processed_file_count} structures.")
            print(f"Results saved to {args.output}")
            print(f"Total execution time: {total_time/60:.1f} minutes ({total_time/3600:.2f} hours)")
            
            # Generate summary statistics if requested
            if args.summary and not combined_df.empty:
                print("\nGenerating summary statistics...")
                
                # Count of each residue type
                residue_counts = combined_df['ResidueType'].value_counts()
                print("\nResidue Type Counts:")
                for res, count in residue_counts.items():
                    print(f"  {res}: {count} ({count/len(combined_df)*100:.1f}%)")
                
                # Secondary structure distribution
                if 'SecondaryStructure' in combined_df.columns:
                    ss_counts = combined_df['SecondaryStructure'].value_counts()
                    print("\nSecondary Structure Distribution:")
                    for ss, count in ss_counts.items():
                        print(f"  {ss}: {count} ({count/len(combined_df)*100:.1f}%)")
                
                # HSE statistics (report for CB which is most relevant for phosphosites)
                if 'HSE_CB_RATIO' in combined_df.columns:
                    print("\nHSE-CB Ratio Statistics (higher ratio = more exposed):")
                    for res in ['S', 'T', 'Y']:
                        res_df = combined_df[combined_df['ResidueType'] == res]
                        if not res_df.empty and not res_df['HSE_CB_RATIO'].isna().all():
                            median_ratio = res_df['HSE_CB_RATIO'].median()
                            mean_ratio = res_df['HSE_CB_RATIO'].mean()
                            print(f"  {res}: Median={median_ratio:.2f}, Mean={mean_ratio:.2f}")
                
                # Cavity volume statistics
                if 'CavityVolume' in combined_df.columns:
                    print("\nCavity Volume Statistics (Å³):")
                    for res in ['S', 'T', 'Y']:
                        res_df = combined_df[combined_df['ResidueType'] == res]
                        if not res_df.empty and not res_df['CavityVolume'].isna().all():
                            median_vol = res_df['CavityVolume'].median()
                            mean_vol = res_df['CavityVolume'].mean()
                            print(f"  {res}: Median={median_vol:.1f}Å³, Mean={mean_vol:.1f}Å³")
                
                # pLDDT statistics
                if 'pLDDT' in combined_df.columns:
                    mean_plddt = combined_df['pLDDT'].mean()
                    median_plddt = combined_df['pLDDT'].median()
                    print(f"\npLDDT Statistics:")
                    print(f"  Mean: {mean_plddt:.1f}")
                    print(f"  Median: {median_plddt:.1f}")
                    
                    # pLDDT distribution by residue type
                    print("\npLDDT by Residue Type:")
                    for res in ['S', 'T', 'Y']:
                        res_df = combined_df[combined_df['ResidueType'] == res]
                        if not res_df.empty:
                            print(f"  {res}: Mean={res_df['pLDDT'].mean():.1f}, Median={res_df['pLDDT'].median():.1f}")
                
                # Summary stats file
                summary_file = f"{os.path.splitext(args.output)[0]}_summary.txt"
                with open(summary_file, 'w') as f:
                    f.write(f"STY Analysis Summary\n")
                    f.write(f"===================\n\n")
                    f.write(f"Total STY sites analyzed: {len(combined_df)}\n")
                    f.write(f"Structures processed: {processed_file_count}\n\n")
                    
                    f.write("Residue Type Counts:\n")
                    for res, count in residue_counts.items():
                        f.write(f"  {res}: {count} ({count/len(combined_df)*100:.1f}%)\n")
                    
                    if 'SecondaryStructure' in combined_df.columns:
                        f.write("\nSecondary Structure Distribution:\n")
                        for ss, count in ss_counts.items():
                            f.write(f"  {ss}: {count} ({count/len(combined_df)*100:.1f}%)\n")
                    
                    # Add more summary statistics to the file
                    if 'CavityVolume' in combined_df.columns:
                        f.write("\nCavity Volume Statistics (Å³):\n")
                        for res in ['S', 'T', 'Y']:
                            res_df = combined_df[combined_df['ResidueType'] == res]
                            if not res_df.empty and not res_df['CavityVolume'].isna().all():
                                f.write(f"  {res}: Mean={res_df['CavityVolume'].mean():.1f}Å³, Median={res_df['CavityVolume'].median():.1f}Å³\n")
                    
                    # HSE statistics
                    if 'HSE_CB_RATIO' in combined_df.columns:
                        f.write("\nHSE-CB Ratio Statistics:\n")
                        for res in ['S', 'T', 'Y']:
                            res_df = combined_df[combined_df['ResidueType'] == res]
                            if not res_df.empty and not res_df['HSE_CB_RATIO'].isna().all():
                                f.write(f"  {res}: Mean={res_df['HSE_CB_RATIO'].mean():.2f}, Median={res_df['HSE_CB_RATIO'].median():.2f}\n")
                    
                print(f"Summary statistics saved to {summary_file}")
        else:
            print("\nNo valid temporary files found to combine.")
    else:
        print("\nNo results were found.")


if __name__ == "__main__":
    main()