import pandas as pd
import numpy as np
import glob
import os
import time
import gc
import multiprocessing as mp
from functools import partial
from Bio.Align import substitution_matrices

# Load BLOSUM62 matrix once
BLOSUM62 = substitution_matrices.load("BLOSUM62")

def calculate_blosum62_score(motif1, motif2):
    """
    Calculate BLOSUM62 similarity score between two motifs.
    Ignores 'X' characters in the match calculation.
    
    Args:
        motif1: First amino acid motif
        motif2: Second amino acid motif
    
    Returns:
        Normalized similarity score
    """
    # Initialize score
    total_score = 0
    valid_positions = 0
    
    # Ensure motifs are of the same length
    min_len = min(len(motif1), len(motif2))
    
    # Calculate score only for positions that don't have X
    for i in range(min_len):
        if motif1[i] != 'X' and motif2[i] != 'X':
            try:
                total_score += BLOSUM62[motif1[i], motif2[i]]
                valid_positions += 1
            except KeyError:
                # Handle case where unusual amino acid codes are used
                pass
    
    # Return normalized score (avoid division by zero)
    if valid_positions > 0:
        return total_score / valid_positions
    else:
        return 0.0

def read_and_prep_file(file_info, motif1_col='Motif1', motif2_col='Motif2', similarity_threshold=0.3):
    """
    Process a single file in parallel with BLOSUM62 score calculation
    
    Args:
        file_info: Tuple of (file_index, file_path)
        motif1_col: Column name for first motif
        motif2_col: Column name for second motif
        similarity_threshold: Minimum similarity score to keep a row
        
    Returns:
        Dictionary with processed dataframe and stats
    """
    file_index, file_path = file_info
    try:
        file_start_time = time.time()
        file_size_mb = os.path.getsize(file_path) / (1024 * 1024)
        
        # Read the file
        df = pd.read_feather(file_path)
        total_rows = len(df)
        
        # Calculate BLOSUM62 scores with optimized approach
        print(f"  Calculating BLOSUM62 scores for file {file_index+1}...", flush=True)
        score_start_time = time.time()
        
        # Process in chunks to provide better progress updates
        chunk_size = 100000
        df['blosum62_score'] = 0.0  # Initialize the score column
        
        for chunk_start in range(0, total_rows, chunk_size):
            chunk_end = min(chunk_start + chunk_size, total_rows)
            if chunk_start % chunk_size == 0:
                print(f"    File {file_index+1}: Processing rows {chunk_start:,} to {chunk_end:,} of {total_rows:,}...", flush=True)
            
            # Use optimized vectorization with pandas apply on the chunk
            chunk_df = df.iloc[chunk_start:chunk_end]
            scores = chunk_df.apply(
                lambda row: calculate_blosum62_score(row[motif1_col], row[motif2_col]),
                axis=1
            )
            df.iloc[chunk_start:chunk_end, df.columns.get_loc('blosum62_score')] = scores
            
            # Periodically collect garbage to reduce memory pressure
            if chunk_end % (chunk_size * 5) == 0:
                gc.collect()
        
        print(f"    File {file_index+1}: Completed scoring all {total_rows:,} rows", flush=True)
        
        # Calculate statistics for different thresholds
        initial_rows = len(df)
        matches_above_03 = (df['blosum62_score'] >= 0.3).sum()
        matches_above_05 = (df['blosum62_score'] >= 0.5).sum()
        matches_above_075 = (df['blosum62_score'] >= 0.75).sum()
        matches_above_09 = (df['blosum62_score'] >= 0.9).sum()
        
        # Calculate percentage of matches above each threshold
        pct_above_03 = (matches_above_03 / initial_rows * 100) if initial_rows > 0 else 0
        pct_above_05 = (matches_above_05 / initial_rows * 100) if initial_rows > 0 else 0
        pct_above_075 = (matches_above_075 / initial_rows * 100) if initial_rows > 0 else 0
        pct_above_09 = (matches_above_09 / initial_rows * 100) if initial_rows > 0 else 0
        
        # Filter rows based on similarity threshold
        df = df[df['blosum62_score'] >= similarity_threshold]
        filtered_rows = initial_rows - len(df)
        
        score_time = time.time() - score_start_time
        
        # Return the dataframe along with info for progress reporting
        return {
            'df': df, 
            'file_index': file_index,
            'file_path': file_path,
            'size_mb': file_size_mb,
            'initial_rows': initial_rows,
            'filtered_rows': filtered_rows,
            'remaining_rows': len(df),
            'matches_above_03': matches_above_03,
            'matches_above_05': matches_above_05,
            'matches_above_075': matches_above_075,
            'matches_above_09': matches_above_09,
            'pct_above_03': pct_above_03,
            'pct_above_05': pct_above_05,
            'pct_above_075': pct_above_075,
            'pct_above_09': pct_above_09,
            'score_time': score_time,
            'total_time': time.time() - file_start_time
        }
    except Exception as e:
        import traceback
        traceback_str = traceback.format_exc()
        return {
            'df': None,
            'file_index': file_index, 
            'file_path': file_path,
            'error': f"{str(e)}\n{traceback_str}"
        }

def combine_feather_files_with_blosum62(
    input_dir, 
    output_file, 
    id1_col='ID1', 
    id2_col='ID2',
    motif1_col='Motif1', 
    motif2_col='Motif2',
    similarity_threshold=0.3,
    batch_size=40, 
    memory_threshold_gb=250, 
    workers=64
):
    """
    Combines feather files with BLOSUM62 scoring and filtering.
    
    Args:
        input_dir: Directory containing feather files
        output_file: Path to output combined feather file
        id1_col: Column name for first ID
        id2_col: Column name for second ID
        motif1_col: Column name for first motif
        motif2_col: Column name for second motif
        similarity_threshold: Minimum similarity score to keep a row
        batch_size: Number of files to process at once
        memory_threshold_gb: When memory usage exceeds this value, write intermediate results
        workers: Number of parallel workers
    """
    start_time = time.time()
    
    # Get all feather files
    feather_files = glob.glob(os.path.join(input_dir, "*.feather"))
    print(f"Found {len(feather_files)} feather files")
    
    # Sample first file to get column names
    print("Reading schema from first file...")
    sample_df = pd.read_feather(feather_files[0])
    print(f"Columns in file: {sample_df.columns.tolist()}")
    print(f"First column (ID1): {id1_col}")
    print(f"Second column (ID2): {id2_col}")
    print(f"Motif columns: {motif1_col}, {motif2_col}")
    
    # Verify columns exist
    required_cols = [id1_col, id2_col, motif1_col, motif2_col]
    for col in required_cols:
        if col not in sample_df.columns:
            raise ValueError(f"Required column '{col}' not found in the feather files")
    
    # Free memory
    del sample_df
    gc.collect()
    
    # Create a temporary directory for intermediate files
    temp_dir = os.path.dirname(output_file)
    if not temp_dir:
        temp_dir = "."  # Use current directory if no directory specified
    os.makedirs(temp_dir, exist_ok=True)
    intermediate_file = os.path.join(temp_dir, "intermediate_result.feather")
    
    # Variables to track progress
    all_filtered_rows = None
    total_rows_processed = 0
    total_rows_filtered = 0
    total_matches_above_03 = 0
    total_matches_above_05 = 0
    total_matches_above_075 = 0
    total_matches_above_09 = 0
    unique_rows_kept = 0
    intermediate_flush_count = 0
    
    # Create a summary file to store per-file statistics
    summary_file = os.path.join(temp_dir, "similarity_summary.csv")
    with open(summary_file, 'w') as f:
        f.write("file_name,total_matches,matches_above_03,pct_above_03,matches_above_05,pct_above_05,"
                "matches_above_075,pct_above_075,matches_above_09,pct_above_09\n")
    
    # Process files in batches
    total_batches = (len(feather_files) + batch_size - 1) // batch_size
    
    # Create a pool of workers for parallel file reading and processing
    pool = mp.Pool(processes=workers)
    
    for i in range(0, len(feather_files), batch_size):
        batch_start_time = time.time()
        batch_end = min(i + batch_size, len(feather_files))
        
        # Create file info list for this batch
        file_infos = [(i+idx, file_path) for idx, file_path in enumerate(feather_files[i:batch_end])]
        
        print(f"Processing batch {i//batch_size + 1}/{total_batches}: Files {i+1} to {batch_end}")
        print(f"  Reading and processing {len(file_infos)} files in parallel...")
        
        # Create a partial function with our similarity threshold
        process_func = partial(
            read_and_prep_file,
            motif1_col=motif1_col,
            motif2_col=motif2_col,
            similarity_threshold=similarity_threshold
        )
        
        # Read and process files in parallel
        read_start_time = time.time()
        batch_results = pool.map(process_func, file_infos)
        read_time = time.time() - read_start_time
        print(f"  Parallel reading and scoring completed in {read_time:.2f}s")
        
        # Process results and collect valid dataframes
        batch_dfs = []
        for result in batch_results:
            if result.get('df') is not None:
                # Print progress info for successful processing with detailed thresholds
                print(f"  File {result['file_index']+1}/{len(feather_files)}: {os.path.basename(result['file_path'])} "
                      f"({result['size_mb']:.2f} MB) - {result['initial_rows']:,} total matches checked")
                print(f"    Similarity distribution:")
                print(f"      ≥ 0.3: {result['matches_above_03']:,} matches ({result['pct_above_03']:.2f}%)")
                print(f"      ≥ 0.5: {result['matches_above_05']:,} matches ({result['pct_above_05']:.2f}%)")
                print(f"      ≥ 0.75: {result['matches_above_075']:,} matches ({result['pct_above_075']:.2f}%)")
                print(f"      ≥ 0.9: {result['matches_above_09']:,} matches ({result['pct_above_09']:.2f}%)")
                print(f"    Processing time: {result['score_time']:.2f}s (total: {result['total_time']:.2f}s)")
                # Add file statistics to summary file
                with open(summary_file, 'a') as f:
                    f.write(f"{os.path.basename(result['file_path'])},{result['initial_rows']},"
                            f"{result['matches_above_03']},{result['pct_above_03']:.2f},"
                            f"{result['matches_above_05']},{result['pct_above_05']:.2f},"
                            f"{result['matches_above_075']},{result['pct_above_075']:.2f},"
                            f"{result['matches_above_09']},{result['pct_above_09']:.2f}\n")
                
                # Add to running totals
                batch_dfs.append(result['df'])
                total_rows_processed += result['initial_rows']
                total_rows_filtered += result['filtered_rows']
                total_matches_above_03 += result['matches_above_03']
                total_matches_above_05 += result['matches_above_05']
                total_matches_above_075 += result['matches_above_075']
                total_matches_above_09 += result['matches_above_09']
            else:
                print(f"  ERROR processing file {result['file_index']+1}/{len(feather_files)}: "
                      f"{os.path.basename(result['file_path'])} - {result.get('error', 'Unknown error')}")
        
        if not batch_dfs:
            print("  No files were successfully processed in this batch. Skipping.")
            continue
            
        print(f"  Concatenating {len(batch_dfs)} files...")
        concat_start_time = time.time()
        batch_combined = pd.concat(batch_dfs, ignore_index=True)
        print(f"  Concatenation complete. Batch size: {len(batch_combined):,} rows (took {time.time() - concat_start_time:.2f}s)")
        
        # Free memory
        del batch_dfs
        gc.collect()
        
        # Convert ID columns to categorical for faster deduplication if they're string-like
        if batch_combined[id1_col].dtype == 'object':
            batch_combined[id1_col] = batch_combined[id1_col].astype('category')
        if batch_combined[id2_col].dtype == 'object':
            batch_combined[id2_col] = batch_combined[id2_col].astype('category')
        
        # We'll skip deduplication within the batch and just keep all rows
        # to perform a single deduplication at the end
        batch_unique = batch_combined
        print(f"  Rows in batch: {len(batch_unique):,}")
        
        # Free memory
        del batch_combined
        gc.collect()
        
        # Combine with previous results
        if all_filtered_rows is None:
            all_filtered_rows = batch_unique
            print(f"  First batch processed: {len(all_filtered_rows):,} unique rows")
        else:
            # Load from intermediate file if we wrote one
            if os.path.exists(intermediate_file) and all_filtered_rows is None:
                print(f"  Loading previous intermediate results from {intermediate_file}...")
                load_start = time.time()
                all_filtered_rows = pd.read_feather(intermediate_file)
                print(f"  Loaded {len(all_filtered_rows):,} rows from intermediate file in {time.time() - load_start:.2f}s")
                
                # Convert ID columns to categorical for faster deduplication
                if all_filtered_rows[id1_col].dtype == 'object':
                    all_filtered_rows[id1_col] = all_filtered_rows[id1_col].astype('category')
                if all_filtered_rows[id2_col].dtype == 'object':
                    all_filtered_rows[id2_col] = all_filtered_rows[id2_col].astype('category')
            
            # Just concatenate without deduplication (will be done at the end)
            print("  Merging with master dataset...")
            merge_start_time = time.time()
            prev_len = len(all_filtered_rows)
            
            all_filtered_rows = pd.concat([all_filtered_rows, batch_unique], ignore_index=True)
            concat_time = time.time() - merge_start_time
            print(f"  Concatenation complete in {concat_time:.2f}s.")
        
        unique_rows_kept = len(all_filtered_rows)
        memory_usage_gb = all_filtered_rows.memory_usage(deep=True).sum() / (1024**3)
        print(f"  Total unique rows so far: {unique_rows_kept:,}")
        print(f"  Memory usage of combined dataset: {memory_usage_gb:.2f} GB")
        
        # Free memory
        del batch_unique
        gc.collect()
        
        # Check if we need to flush to disk to save memory
        if memory_usage_gb > memory_threshold_gb:
            print(f"  Memory usage ({memory_usage_gb:.2f} GB) exceeds threshold ({memory_threshold_gb} GB)")
            print(f"  Writing intermediate results to disk...")
            flush_start = time.time()
            all_filtered_rows.to_feather(intermediate_file)
            print(f"  Intermediate results written in {time.time() - flush_start:.2f}s")
            # Clear memory
            del all_filtered_rows
            all_filtered_rows = None
            gc.collect()
            intermediate_flush_count += 1
        
        batch_time = time.time() - batch_start_time
        files_remaining = len(feather_files) - (i + len(file_infos))
        if i > 0 and files_remaining > 0:
            avg_time_per_batch = batch_time
            est_time_remaining = avg_time_per_batch * (files_remaining / batch_size)
            print(f"  Batch completed in {batch_time:.2f}s. Estimated time remaining: {est_time_remaining/60:.2f} minutes")
        else:
            print(f"  Batch completed in {batch_time:.2f}s")
        
        print(f"  Current progress: {(i+len(file_infos))/len(feather_files)*100:.2f}% complete")
        print("-" * 80)
    
    # Close the multiprocessing pool
    pool.close()
    pool.join()
    
    # Final deduplication if we wrote intermediate files
    if intermediate_flush_count > 0 and all_filtered_rows is None:
        print(f"Loading final intermediate results for completion...")
        all_filtered_rows = pd.read_feather(intermediate_file)
    
    # Perform the final deduplication on the entire dataset
    print("Performing final deduplication across all data...")
    final_dedup_start = time.time()
    initial_row_count = len(all_filtered_rows)
    all_filtered_rows = all_filtered_rows.drop_duplicates(subset=[id1_col, id2_col], keep='first')
    total_dupes_removed = initial_row_count - len(all_filtered_rows)
    print(f"Removed {total_dupes_removed:,} duplicates in {time.time() - final_dedup_start:.2f}s")
    
    # Write result to feather file
    print(f"Writing {len(all_filtered_rows):,} unique rows to {output_file}...")
    write_start_time = time.time()
    # Use PyArrow engine for faster writing
    all_filtered_rows.reset_index(drop=True).to_feather(output_file, compression=None)
    write_time = time.time() - write_start_time
    print(f"Write completed in {write_time:.2f}s")
    
    # Clean up intermediate file
    if os.path.exists(intermediate_file):
        os.remove(intermediate_file)
        print(f"Removed intermediate file")
    
    elapsed_time = time.time() - start_time
    # Copy the summary file to the output directory
    final_summary_file = os.path.splitext(output_file)[0] + "_summary.csv"
    try:
        import shutil
        shutil.copy2(summary_file, final_summary_file)
        print(f"Similarity statistics saved to {final_summary_file}")
    except Exception as e:
        print(f"Could not copy summary file: {e}")
    
    # Calculate percentages for the grand total
    pct_above_03 = (total_matches_above_03 / total_rows_processed * 100) if total_rows_processed > 0 else 0
    pct_above_05 = (total_matches_above_05 / total_rows_processed * 100) if total_rows_processed > 0 else 0
    pct_above_075 = (total_matches_above_075 / total_rows_processed * 100) if total_rows_processed > 0 else 0
    pct_above_09 = (total_matches_above_09 / total_rows_processed * 100) if total_rows_processed > 0 else 0
    
    print("=" * 80)
    print(f"Finished! Combined data written to {output_file}")
    print(f"Total rows processed: {total_rows_processed:,}")
    print(f"Similarity distribution across all files:")
    print(f"  ≥ 0.3:  {total_matches_above_03:,} matches ({pct_above_03:.2f}%)")
    print(f"  ≥ 0.5:  {total_matches_above_05:,} matches ({pct_above_05:.2f}%)")
    print(f"  ≥ 0.75: {total_matches_above_075:,} matches ({pct_above_075:.2f}%)")
    print(f"  ≥ 0.9:  {total_matches_above_09:,} matches ({pct_above_09:.2f}%)")
    print(f"Total unique rows after filtering and deduplication: {len(all_filtered_rows):,}")
    print(f"Intermediate file flushes: {intermediate_flush_count}")
    print(f"Total time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
    print(f"Detailed per-file statistics saved to {final_summary_file}")

if __name__ == "__main__":
    # Example usage
    combine_feather_files_with_blosum62(
        input_dir="/k3s-storage/extra-ssd/dvanderwall/batch_output/", 
        output_file="combined_STY_motif_matches_with_blosum.feather",
        id1_col='ID1',
        id2_col='ID2',
        motif1_col='Motif1',
        motif2_col='Motif2',
        similarity_threshold=0.3,
        batch_size=96,
        memory_threshold_gb=350,
        workers=96
    )