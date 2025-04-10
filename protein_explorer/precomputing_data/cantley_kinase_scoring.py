import pandas as pd
import numpy as np
from pathlib import Path
import re
from tqdm import tqdm
import os
import shutil

# Function to load the kinase data
def load_kinase_data(file_path):
    """
    Load kinase position-specific scoring data from Excel or CSV file.
    
    Args:
        file_path: Path to the Excel or CSV file
        
    Returns:
        DataFrame with kinase data
    """
    try:
        # Check file extension
        if file_path.lower().endswith('.xlsx') or file_path.lower().endswith('.xls'):
            try:
                # Try to load the Excel file
                print("Loading Excel file...")
                df = pd.read_excel(file_path, index_col=0)
                return df
            except ImportError:
                print("Excel loading failed. Do you have openpyxl installed?")
                print("Try: python -m pip install openpyxl")
                print("Or convert your Excel file to CSV and use that instead.")
                
                # Ask if there's a CSV version available
                csv_path = input("Enter path to CSV version of the file (or press Enter to exit): ")
                if not csv_path:
                    return None
                file_path = csv_path
        
        # If we get here, either it was a CSV file to begin with, or we're falling back to CSV
        if file_path.lower().endswith('.csv'):
            print("Loading CSV file...")
            df = pd.read_csv(file_path, index_col=0)
            return df
        else:
            print(f"Unsupported file format: {file_path}")
            return None
    except Exception as e:
        print(f"Error loading kinase data: {e}")
        return None

# Function to parse the column names into position and amino acid
def parse_column_names(df):
    """
    Parse column names to extract position and amino acid information.
    
    Args:
        df: DataFrame with kinase data
        
    Returns:
        Dictionary mapping (position, amino_acid) to column name
    """
    position_aa_map = {}
    
    for col in df.columns:
        # Extract position and amino acid using regex
        match = re.match(r'(-?\d+)([A-Za-z])', col)
        if match:
            position = int(match.group(1))
            amino_acid = match.group(2)
            position_aa_map[(position, amino_acid)] = col
    
    return position_aa_map

# Function to create a matrix representation of kinase data for faster processing
def create_kinase_matrix(kinase_data, position_aa_map):
    """
    Create a 3D matrix for efficient scoring:
    - First dimension: kinases (corresponding to rows in kinase_data)
    - Second dimension: positions (-5 to +4, excluding 0)
    - Third dimension: amino acids (A-Z)
    
    Returns:
        - 3D numpy array
        - Mapping of positions to indices
        - Mapping of amino acids to indices
    """
    # Define positions (excluding 0) and amino acids
    positions = list(range(-5, 0)) + list(range(1, 5))
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # 20 standard amino acids
    
    # Create mappings for positions and amino acids
    pos_to_idx = {pos: idx for idx, pos in enumerate(positions)}
    aa_to_idx = {aa: idx for idx, aa in enumerate(amino_acids)}
    
    # Print some debug info
    print(f"Processing {len(positions)} positions and {len(amino_acids)} amino acids")
    print(f"Matrix shape will be: ({len(kinase_data)}, {len(positions)}, {len(amino_acids)})")
    
    # Check and report on missing positions or amino acids in the data
    found_pos_aa = set()
    for pos, aa in position_aa_map.keys():
        found_pos_aa.add((pos, aa.upper()))
    
    # Initialize the matrix with ones (neutral value for multiplication)
    # Shape: (num_kinases, num_positions, num_amino_acids)
    matrix = np.ones((len(kinase_data), len(positions), len(amino_acids)))
    
    # Fill the matrix with values from kinase_data
    filled_count = 0
    for (pos, aa), col in position_aa_map.items():
        if pos in pos_to_idx and aa.upper() in aa_to_idx:
            pos_idx = pos_to_idx[pos]
            aa_idx = aa_to_idx[aa.upper()]
            matrix[:, pos_idx, aa_idx] = kinase_data[col].values
            filled_count += 1
    
    print(f"Filled {filled_count} position-amino acid combinations in the matrix")
    
    return matrix, pos_to_idx, aa_to_idx

# Function to convert a batch of motifs to a matrix representation for vectorized scoring
def motifs_to_matrix(motifs, pos_to_idx, aa_to_idx):
    """
    Convert a list of motifs to a binary matrix indicating which amino acids are present
    at each position for each motif.
    
    Args:
        motifs: List of motif strings
        pos_to_idx: Mapping of positions to indices
        aa_to_idx: Mapping of amino acids to indices
        
    Returns:
        Binary matrix of shape (num_motifs, num_positions, num_amino_acids)
        Array indicating number of valid positions per motif
    """
    # Positions in the motif string (excluding position 0)
    motif_positions = list(range(-5, 0)) + list(range(1, 5))
    motif_indices = list(range(0, 5)) + list(range(6, 10))  # Indices in 10-char string, skipping pos 0
    
    # Initialize matrix with zeros
    # Shape: (num_motifs, num_positions, num_amino_acids)
    matrix = np.zeros((len(motifs), len(pos_to_idx), len(aa_to_idx)), dtype=bool)
    
    # Count of valid positions per motif (for filtering out those with too many Xs)
    valid_positions = np.zeros(len(motifs), dtype=int)
    
    # Track statistics for debugging
    invalid_length_count = 0
    invalid_aa_count = 0
    x_position_count = 0
    
    # Fill the matrix
    for m_idx, motif in enumerate(motifs):
        if len(motif) != 10:
            invalid_length_count += 1
            continue  # Skip invalid motifs
            
        for str_idx, pos in zip(motif_indices, motif_positions):
            aa = motif[str_idx].upper()
            
            # Skip X positions
            if aa == 'X':
                x_position_count += 1
                continue
                
            if aa in aa_to_idx and pos in pos_to_idx:
                pos_idx = pos_to_idx[pos]
                aa_idx = aa_to_idx[aa]
                matrix[m_idx, pos_idx, aa_idx] = True
                valid_positions[m_idx] += 1
            else:
                invalid_aa_count += 1
    
    # Print statistics
    if invalid_length_count > 0:
        print(f"Warning: Skipped {invalid_length_count} motifs with invalid length")
    
    if invalid_aa_count > 0:
        print(f"Warning: Found {invalid_aa_count} amino acids not in our reference set")
    
    if x_position_count > 0:
        print(f"Info: Skipped {x_position_count} 'X' positions during matrix creation")
        
    return matrix, valid_positions

# Vectorized scoring function with file-based processing to save memory
def score_motifs_vectorized(motifs, kinase_data, batch_size=5000, temp_dir="./temp_results"):
    """
    Score a list of motifs against all kinases using vectorized operations.
    Writes batch results to temporary files to save memory.
    
    Args:
        motifs: List of motif strings
        kinase_data: DataFrame with kinase scoring data
        batch_size: Number of motifs to process in each batch
        temp_dir: Directory to store temporary batch results
        
    Returns:
        Path to the final results file
    """
    # Create temp directory if it doesn't exist
    temp_path = Path(temp_dir)
    if not temp_path.exists():
        temp_path.mkdir(parents=True)
    
    # Parse column names and create kinase matrix
    position_aa_map = parse_column_names(kinase_data)
    kinase_matrix, pos_to_idx, aa_to_idx = create_kinase_matrix(kinase_data, position_aa_map)
    
    # Check for duplicate motifs
    motif_counts = {}
    for m in motifs:
        motif_counts[m] = motif_counts.get(m, 0) + 1
    
    duplicate_motifs = [m for m, count in motif_counts.items() if count > 1]
    if duplicate_motifs:
        print(f"Found {len(duplicate_motifs)} duplicate motifs out of {len(motifs)} total motifs.")
        # Process unique motifs to save computation
        unique_motifs = list(set(motifs))
        print(f"Will process {len(unique_motifs)} unique motifs and then map results back.")
    else:
        unique_motifs = motifs
    
    # Write header to results file
    final_results_path = os.path.join(temp_dir, "final_results.csv")
    with open(final_results_path, 'w') as f:
        # Write header (motif name and kinase names)
        f.write("motif," + ",".join(kinase_data.index) + "\n")
    
    # Process in batches to manage memory
    total_motifs = len(unique_motifs)
    num_batches = (total_motifs + batch_size - 1) // batch_size
    
    print(f"Processing {total_motifs} unique motifs in {num_batches} batches of size {batch_size}")
    
    # Dictionary to store the results for duplicate motifs
    results_map = {}
    
    # Process each batch
    for batch_idx in tqdm(range(num_batches), desc="Processing motifs in batches"):
        start_idx = batch_idx * batch_size
        end_idx = min(start_idx + batch_size, total_motifs)
        
        batch_motifs = unique_motifs[start_idx:end_idx]
        
        # Convert motifs to matrix representation
        motif_matrix, valid_positions = motifs_to_matrix(batch_motifs, pos_to_idx, aa_to_idx)
        
        # Calculate scores using vectorized operations
        # We'll process one position at a time to reduce memory usage
        batch_scores = np.ones((len(batch_motifs), len(kinase_data)))
        
        # For each position
        for pos_idx in range(len(pos_to_idx)):
            # Create temporary scores for this position
            pos_scores = np.ones((len(batch_motifs), len(kinase_data)))
            
            # For each amino acid
            for aa_idx in range(len(aa_to_idx)):
                # Get motifs that have this amino acid at this position
                mask = motif_matrix[:, pos_idx, aa_idx]
                if np.any(mask):
                    # Get scores for this amino acid at this position for all kinases
                    aa_scores = kinase_matrix[:, pos_idx, aa_idx]
                    
                    # Apply to selected motifs (broadcasting)
                    pos_scores[mask] *= aa_scores
            
            # Multiply into the overall scores
            batch_scores *= pos_scores
        
        # Set scores to None for motifs with too few valid positions
        invalid_mask = valid_positions < 4
        batch_scores[invalid_mask] = np.nan
        
        # Create a batch results DataFrame
        batch_df = pd.DataFrame(
            batch_scores,
            index=batch_motifs,
            columns=kinase_data.index
        )
        
        # Write batch results to CSV
        batch_file = os.path.join(temp_dir, f"batch_{batch_idx}.csv")
        batch_df.to_csv(batch_file)
        
        # Store the results for unique motifs that might be duplicated later
        if duplicate_motifs:
            for motif in batch_motifs:
                results_map[motif] = batch_df.loc[motif].values
        
        # Clear memory
        del batch_scores, batch_df, motif_matrix, valid_positions
        
    # If we have duplicates, create the final results by mapping
    if duplicate_motifs:
        # Write all motifs (including duplicates) to final results file
        with open(final_results_path, 'a') as f:
            for motif in tqdm(motifs, desc="Writing final results (with duplicates)"):
                # Get the scores for this motif from our map
                scores = results_map.get(motif, [np.nan] * len(kinase_data))
                # Write as CSV line
                f.write(motif + "," + ",".join(map(str, scores)) + "\n")
    else:
        # Concatenate all batch files
        print("Combining batch results...")
        with open(final_results_path, 'a') as outfile:
            for batch_idx in range(num_batches):
                batch_file = os.path.join(temp_dir, f"batch_{batch_idx}.csv")
                with open(batch_file, 'r') as infile:
                    # Skip header for all but the first file
                    next(infile)  # Skip header
                    for line in infile:
                        outfile.write(line)
                
                # Optionally, remove the temporary batch file
                os.remove(batch_file)
    
    print(f"Results saved to {final_results_path}")
    return final_results_path

def load_motifs_from_feather(file_path):
    """
    Load motifs from the Motif column of a Feather file.
    
    Args:
        file_path: Path to the Feather file
        
    Returns:
        List of motifs
    """
    try:
        import pyarrow.feather as feather
        
        # Load the Feather file
        df = feather.read_feather(file_path)
        
        # Check if Motif column exists
        if 'Fourmer' not in df.columns:
            print(f"Error: 'Fourmer' column not found in {file_path}")
            print(f"Available columns: {df.columns.tolist()}")
            return None
        
        # Extract motifs from the Motif column
        motifs = df['Fourmer'].tolist()
        
        # Print the first ten motifs
        print("First ten motifs:")
        for i, motif in enumerate(motifs[:10]):
            print(f"{i+1}. {motif}")
        
        return motifs
    except Exception as e:
        print(f"Error loading motifs from Feather file: {e}")
        return None

def validate_motifs_list(motifs):
    """
    Validate a list of motifs.
    
    Args:
        motifs: List of motifs to validate
        
    Returns:
        True if all motifs are valid, False otherwise
    """
    valid = True
    for i, motif in enumerate(motifs):
        if len(motif) != 10:
            print(f"Error: Motif at index {i} ('{motif}') has length {len(motif)}, expected 10")
            valid = False
            if i >= 10:  # Only show up to 10 errors
                print("Too many errors, stopping validation...")
                break
        else:
            # Check if position 0 (index 5) is S, T, or Y
            if motif[5].upper() not in ['S', 'T', 'Y']:
                print(f"Warning: Motif at index {i} ('{motif}') has {motif[5]} at position 0, expected S, T, or Y")
            
            for aa in motif:
                if not aa.isalpha():
                    print(f"Error: Motif at index {i} ('{motif}') contains non-alphabetic character '{aa}'")
                    valid = False
                    if i >= 10:  # Only show up to 10 errors
                        print("Too many errors, stopping validation...")
                        break
    
    return valid

def main():
    # Define file paths
    kinase_data_file = "F:/Kinome/Tyr_Kinase_Densitomitries.xlsx"
    motifs_file = "F:/Kinome/Precomputed_Data/All_Y_Motifs.feather"
    output_file = "F:/Kinome/Tyr_motif_scores.csv"
    temp_dir = "F:/Kinome/temp_results"
    
    # Load kinase data
    print("Loading kinase data...")
    kinase_data = load_kinase_data(kinase_data_file)
    
    if kinase_data is None:
        print("Failed to load kinase data. Exiting.")
        return
    
    print(f"Loaded data for {len(kinase_data)} kinases.")
    
    # Load motifs from Feather file
    print(f"Loading motifs from {motifs_file}...")
    motifs = load_motifs_from_feather(motifs_file)
    
    if motifs is None:
        print("Failed to load motifs. Exiting.")
        return
    
    print(f"Loaded {len(motifs)} motifs.")
    
    # Check for duplicate motifs
    unique_motifs = set(motifs)
    if len(unique_motifs) < len(motifs):
        print(f"Note: Found {len(motifs) - len(unique_motifs)} duplicate motifs.")
    
    # Validate motifs
    print("Validating motifs...")
    if not validate_motifs_list(motifs):
        print("Motif validation failed. Please check your motifs.")
        return
    
    # Prompt user for batch size
    try:
        batch_size = int(input("Enter batch size (default: 5000): ") or 5000)
    except ValueError:
        print("Invalid input. Using default batch size of 5000.")
        batch_size = 5000
    
    # Score motifs using file-based vectorized approach
    print("Scoring motifs using vectorized approach with file-based processing...")
    results_file = score_motifs_vectorized(
        motifs, 
        kinase_data, 
        batch_size=batch_size,
        temp_dir=temp_dir
    )
    
    # If the output file is different from the results file, copy it
    if results_file != output_file:
        print(f"Copying results to final destination: {output_file}")
        shutil.copy(results_file, output_file)
    
    print(f"Done! Results saved to {output_file}")

if __name__ == "__main__":
    main()