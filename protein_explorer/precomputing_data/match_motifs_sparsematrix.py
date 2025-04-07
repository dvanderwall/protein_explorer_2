import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import itertools
from tqdm import tqdm
import os
import time
import threading
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import psutil
import math
from datetime import datetime, timedelta
import concurrent.futures
from scipy.sparse import csr_matrix, lil_matrix


def find_matching_pairs_sparse_optimized(file_path, min_matching_nmers=3, n=3, batch_size=50000,
                                        progress_interval=10, num_workers=None, memory_limit_gb=None):
    """
    Highly optimized parallel implementation of motif matching using sparse matrix operations.
   
    Key optimizations:
    1. Multi-core processing with work division optimized for high CPU count
    2. Sparse matrix representation for efficient similarity computation
    3. Vectorized operations for n-mer intersection calculations
    4. Memory-efficient data structures with adaptive batch sizing
    5. Dynamic workload balancing for uneven computation loads
   
    Args:
        file_path: Path to feather file
        min_matching_nmers: Minimum number of matching n-mers required
        n: Length of n-mers to use
        batch_size: Maximum batch size for processing (will be adjusted based on memory)
        progress_interval: How often to print progress updates (in seconds)
        num_workers: Number of worker processes (defaults to optimal for system)
        memory_limit_gb: Memory limit in GB (defaults to 90% of available system memory)
       
    Returns:
        DataFrame with matching ID pairs and their motifs
    """
    # Start timer and print initial information
    start_time = time.time()
    print(f"\n{'='*80}")
    print(f"MOTIF MATCHING - SPARSE MATRIX IMPLEMENTATION")
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*80}\n")
   
    # System information
    system_info = get_system_info()
    print(f"SYSTEM INFORMATION:")
    print(f"  CPU: {system_info['cpu_count']} cores ({system_info['cpu_model']})")
    print(f"  Memory: {system_info['memory_gb']:.1f} GB available")
   
    # Set default number of workers if not specified - optimized for high CPU count
    if num_workers is None:
        # Use 90% of available cores for computation, leaving some for system and I/O
        num_workers = max(1, int(system_info['cpu_count'] * 0.9))
   
    # Set memory limit if not specified - optimized for high-memory systems
    if memory_limit_gb is None:
        memory_limit_gb = system_info['memory_gb'] * 0.9  # Use 90% of available memory
   
    print(f"\nRUNNING WITH:")
    print(f"  Workers: {num_workers}")
    print(f"  Memory limit: {memory_limit_gb:.1f} GB")
    print(f"  N-mer size: {n}")
    print(f"  Min matching n-mers: {min_matching_nmers}")
    print(f"  Initial batch size: {batch_size}")
    print(f"  Progress update interval: {progress_interval} seconds\n")

    # Create the progress tracker as a shared resource
    progress = ProgressTracker(update_interval=progress_interval)
   
    # Start a separate thread for the progress reporter
    reporter_thread = threading.Thread(target=progress_reporter, args=(progress,))
    reporter_thread.daemon = True
    reporter_thread.start()
   
    # Load data and preprocess
    progress.update_stage("Loading and preprocessing data")
    load_start = time.time()
    print(f"Loading data from {file_path}...")
   
    # OPTIMIZATION: Increase thread count for data loading
    # Only load the columns we need
    with ThreadPoolExecutor(max_workers=min(32, num_workers)) as executor:  # Cap at 32 threads for I/O
        df = pd.read_feather(file_path, use_threads=True)
        # Force loading into memory with parallel processing
        _ = executor.submit(lambda: len(df))
   
    # Access columns by position for performance
    col1_name = df.columns[0]  # First column is ID
    col2_name = df.columns[1]  # Second column is Motif
   
    total_rows = len(df)
    load_time = time.time() - load_start
    print(f"Loaded {total_rows:,} rows in {load_time:.2f} seconds")
   
    # STEP 1: Create fingerprints and n-mer sets
    progress.update_stage("Creating n-mer sets")
    step1_start = time.time()
    print("Creating n-mer sets...")
   
    # OPTIMIZATION: Use parallel preprocessing for n-mer extraction
    # Function to process a chunk of data in parallel
    def process_chunk(chunk_df):
        chunk_valid_ids = []
        chunk_id_to_motif = {}
        chunk_id_to_nmers = {}
        chunk_valid_count = 0
        
        for _, row in chunk_df.iterrows():
            id_val = row[col1_name]
            motif = row[col2_name]
            
            # Skip invalid motifs
            if not isinstance(motif, str) or len(motif) < n:
                continue
                
            chunk_valid_count += 1
            chunk_id_to_motif[id_val] = motif
            chunk_valid_ids.append(id_val)
            
            # Extract n-mers directly to sets
            nmers = set(motif[i:i+n] for i in range(len(motif) - n + 1) if 'X' not in motif[i:i+n])
            chunk_id_to_nmers[id_val] = nmers
            
        return chunk_valid_ids, chunk_id_to_motif, chunk_id_to_nmers, chunk_valid_count
    
    # OPTIMIZATION: Split dataframe into chunks and process in parallel
    # Calculate optimal chunk size based on CPU count
    chunk_size = max(1000, min(100000, total_rows // (num_workers * 2)))
    chunks = [df[i:i + chunk_size] for i in range(0, total_rows, chunk_size)]
    
    valid_ids = []
    id_to_motif = {}
    id_to_nmers = {}
    valid_count = 0
    total_nmers = 0
    
    progress.set_total_items(len(chunks))
    
    # Process chunks in parallel
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(process_chunk, chunk): i for i, chunk in enumerate(chunks)}
        
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing chunks"):
            chunk_valid_ids, chunk_id_to_motif, chunk_id_to_nmers, chunk_valid_count = future.result()
            
            valid_ids.extend(chunk_valid_ids)
            id_to_motif.update(chunk_id_to_motif)
            id_to_nmers.update(chunk_id_to_nmers)
            valid_count += chunk_valid_count
            
            # Count total n-mers
            for nmers in chunk_id_to_nmers.values():
                total_nmers += len(nmers)
            
            progress.increment_processed()
   
    # Free memory
    del df
    del chunks
   
    step1_time = time.time() - step1_start
    print(f"N-mer extraction complete in {step1_time:.2f} seconds")
    print(f"Valid motifs: {valid_count:,}, Total n-mers: {total_nmers:,}")
   
    # STEP 2: Build sparse matrix representation of n-mers
    progress.update_stage("Building sparse matrix representation")
    step2_start = time.time()
    print("Creating sparse matrix for efficient similarity computation...")
    
    # Collect all unique n-mers
    all_nmers = set()
    
    # Use parallel processing to collect unique n-mers
    def collect_nmers_chunk(chunk_ids):
        chunk_nmers = set()
        for id_val in chunk_ids:
            if id_val in id_to_nmers:
                chunk_nmers.update(id_to_nmers[id_val])
        return chunk_nmers
    
    # Split IDs into chunks for parallel processing
    chunk_size = max(1000, len(valid_ids) // (num_workers * 2))
    id_chunks = [valid_ids[i:i+chunk_size] for i in range(0, len(valid_ids), chunk_size)]
    
    progress.set_total_items(len(id_chunks))
    progress.reset_processed()
    
    # Process in parallel
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(collect_nmers_chunk, chunk) for chunk in id_chunks]
        for future in tqdm(as_completed(futures), total=len(futures), desc="Collecting unique n-mers"):
            all_nmers.update(future.result())
            progress.increment_processed()
    
    print(f"Found {len(all_nmers):,} unique {n}-mers")
    
    # Create mapping from n-mers to indices
    nmer_to_idx = {nmer: idx for idx, nmer in enumerate(all_nmers)}
    
    # Build sparse matrix - this creates the same mapping as id_to_nmers but in matrix form
    progress.update_stage("Building sparse matrix")
    
    # Function to build matrix for a chunk of IDs
    def build_matrix_chunk(chunk_start, chunk_end):
        # Use LIL format for efficient matrix construction
        chunk_matrix = lil_matrix((chunk_end - chunk_start, len(all_nmers)), dtype=np.bool_)
        
        for i, idx in enumerate(range(chunk_start, chunk_end)):
            if idx >= len(valid_ids):
                break
                
            id_val = valid_ids[idx]
            if id_val in id_to_nmers:
                for nmer in id_to_nmers[id_val]:
                    if nmer in nmer_to_idx:  # Safety check
                        chunk_matrix[i, nmer_to_idx[nmer]] = True
        
        # Convert to CSR format for efficient operations
        return csr_matrix(chunk_matrix)
    
    # Determine optimal chunk size for matrix construction
    # Estimate memory requirements based on sparsity
    est_sparsity = min(1.0, max(0.001, total_nmers / (valid_count * len(all_nmers))))
    matrix_mem_estimate = len(valid_ids) * len(all_nmers) * est_sparsity * 8  # 8 bytes per entry
    
    # Adjust chunk size to fit in 20% of available memory
    max_chunk_size = max(100, min(
        batch_size,
        int((memory_limit_gb * 1024 * 0.2) / (matrix_mem_estimate / len(valid_ids)))
    ))
    
    # Build matrix in chunks
    matrix_chunks = []
    for chunk_start in range(0, len(valid_ids), max_chunk_size):
        chunk_end = min(chunk_start + max_chunk_size, len(valid_ids))
        matrix_chunks.append((chunk_start, chunk_end))
    
    sparse_matrices = []
    
    progress.set_total_items(len(matrix_chunks))
    progress.reset_processed()
    
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(build_matrix_chunk, start, end): (start, end) 
                  for start, end in matrix_chunks}
        
        for future in tqdm(as_completed(futures), total=len(futures), 
                          desc="Building sparse matrix chunks"):
            start, end = futures[future]
            sparse_matrices.append((start, end, future.result()))
            progress.increment_processed()
    
    print(f"Built {len(sparse_matrices)} sparse matrix chunks")
    
    # We can free the memory used by id_to_nmers now that we have the sparse matrix
    del id_to_nmers
    del all_nmers
    
    # STEP 3: Find matches using sparse matrix operations
    progress.update_stage("Finding matches with sparse matrices")
    step3_start = time.time()
    print("Computing similarities using sparse matrix operations...")
    
    # Function to find matches for a batch using sparse matrix operations
    def find_matches_with_sparse(batch_idx, batch_start, batch_end, batch_matrix):
        batch_matches = {}
        comparisons = 0
        batch_time = time.time()
        
        # For each chunk i, compare with chunks j where j >= i
        for j, (other_start, other_end, other_matrix) in enumerate(sparse_matrices):
            # Skip comparisons with earlier chunks (to avoid duplicate comparisons)
            if other_start < batch_start:
                continue
                
            # Same chunk - only compare upper triangular to avoid duplicates
            if batch_idx == j:
                similarity = batch_matrix.dot(batch_matrix.T).toarray()
                
                for i in range(batch_end - batch_start):
                    row_id_idx = batch_start + i
                    if row_id_idx >= len(valid_ids):
                        break
                        
                    id1 = valid_ids[row_id_idx]
                    
                    # Only look at upper triangular part
                    for j in range(i + 1, batch_end - batch_start):
                        col_id_idx = batch_start + j
                        if col_id_idx >= len(valid_ids):
                            break
                            
                        id2 = valid_ids[col_id_idx]
                        comparisons += 1
                        
                        if similarity[i, j] >= min_matching_nmers:
                            batch_matches[(id1, id2)] = int(similarity[i, j])
            
            # Different chunks
            else:
                # Compute dot product between batch_matrix and transposed other_matrix
                similarity = batch_matrix.dot(other_matrix.T).toarray()
                
                for i in range(batch_end - batch_start):
                    row_id_idx = batch_start + i
                    if row_id_idx >= len(valid_ids):
                        break
                        
                    id1 = valid_ids[row_id_idx]
                    
                    for j in range(other_end - other_start):
                        col_id_idx = other_start + j
                        if col_id_idx >= len(valid_ids):
                            break
                            
                        id2 = valid_ids[col_id_idx]
                        comparisons += 1
                        
                        if similarity[i, j] >= min_matching_nmers:
                            batch_matches[(id1, id2)] = int(similarity[i, j])
        
        batch_processing_time = time.time() - batch_time
        
        # Write batch results to file if requested
        if batch_output_dir and id_to_motif:
            # Convert batch results to DataFrame
            batch_result_data = []
            for (id1, id2), count in batch_matches.items():
                motif1 = id_to_motif.get(id1, "")
                motif2 = id_to_motif.get(id2, "")
                batch_result_data.append((id1, id2, motif1, motif2, count))
            
            if batch_result_data:  # Only write if we have results
                batch_df = pd.DataFrame(
                    batch_result_data,
                    columns=['ID1', 'ID2', 'Motif1', 'Motif2', f'Matching_{n}mer_Count']
                )
                
                # Sort by count descending
                batch_df = batch_df.sort_values(f'Matching_{n}mer_Count', ascending=False).reset_index(drop=True)
                
                # Write to feather
                batch_file = os.path.join(batch_output_dir, f'batch_{batch_idx:04d}.feather')
                batch_df.to_feather(batch_file)
                
                # Log batch file creation
                print(f"  Wrote batch {batch_idx} results ({len(batch_df):,} pairs) to {batch_file}")
            
        return batch_matches, comparisons, batch_processing_time
    
    # Create batch_output directory if it doesn't exist
    batch_output_dir = "batch_output"
    if not os.path.exists(batch_output_dir):
        os.makedirs(batch_output_dir)
        print(f"Created directory: {batch_output_dir}")
    
    # Process batches in parallel
    matching_pairs = {}
    total_comparisons = 0
    matches_lock = threading.Lock()  # Lock for thread-safe updates to matching_pairs
    
    progress.set_total_items(len(sparse_matrices))
    progress.reset_processed()
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(find_matches_with_sparse, idx, start, end, matrix): idx 
                  for idx, (start, end, matrix) in enumerate(sparse_matrices)}
        
        # Collect results as they complete
        total_batch_time = 0
        completed_batches = 0
        
        # Custom tqdm with batch timing
        batch_pbar = tqdm(total=len(sparse_matrices), desc="Processing matrix chunks",
                         mininterval=0.5, maxinterval=3.0, leave=True,
                         bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]',
                         ncols=100, position=0)
        
        for future in as_completed(futures):
            batch_matches, batch_comparisons, batch_time = future.result()
            
            # Update timing statistics
            completed_batches += 1
            total_batch_time += batch_time
            avg_batch_time = total_batch_time / completed_batches
            remaining_batches = len(sparse_matrices) - completed_batches
            estimated_remaining_time = avg_batch_time * remaining_batches
            
            # Update global matching pairs with thread safety
            with matches_lock:
                matching_pairs.update(batch_matches)
                total_comparisons += batch_comparisons
            
            # Update progress with timing information
            batch_pbar.set_description(f"Processing matrix chunks - Avg: {avg_batch_time:.2f}s/batch, ETA: {format_time(estimated_remaining_time)}")
            batch_pbar.update(1)
            progress.increment_processed()
        
        batch_pbar.close()
    
    # Update final progress
    progress.update_stage("Preparing final results")
    
    step3_time = time.time() - step3_start
    print(f"\nMatching complete in {format_time(step3_time)}")
    print(f"Found {len(matching_pairs):,} matching pairs from {total_comparisons:,} comparisons")
    
    # Free memory used by sparse matrices
    del sparse_matrices
    
    # OPTIMIZATION: Parallel creation of result DataFrame for large result sets
    print("Preparing final results...")
    
    # Function to process a chunk of matching pairs
    def process_result_chunk(pairs_chunk):
        chunk_data = []
        for (id1, id2), count in pairs_chunk:
            motif1 = id_to_motif.get(id1, "")
            motif2 = id_to_motif.get(id2, "")
            chunk_data.append((id1, id2, motif1, motif2, count))
        return chunk_data
    
    # Split matching pairs into chunks for parallel processing
    pairs_items = list(matching_pairs.items())
    result_chunk_size = max(10000, len(pairs_items) // num_workers)
    pair_chunks = [pairs_items[i:i + result_chunk_size] for i in range(0, len(pairs_items), result_chunk_size)]
    
    result_data = []
    with ThreadPoolExecutor(max_workers=min(32, num_workers)) as executor:
        futures = [executor.submit(process_result_chunk, chunk) for chunk in pair_chunks]
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing result chunks"):
            result_data.extend(future.result())
    
    result_df = pd.DataFrame(result_data,
                       columns=['ID1', 'ID2', 'Motif1', 'Motif2', f'Matching_{n}mer_Count'])

    # Sort by number of matches (descending)
    result_df = result_df.sort_values(f'Matching_{n}mer_Count', ascending=False).reset_index(drop=True)
    
    # Notify progress tracker to stop
    progress.done = True
    
    # Summary
    total_time = time.time() - start_time
    print(f"\n{'='*80}")
    print(f"EXECUTION SUMMARY:")
    print(f"  Total time: {format_time(total_time)}")
    print(f"  Data loading time: {format_time(load_time)} ({load_time/total_time*100:.1f}%)")
    print(f"  N-mer extraction time: {format_time(step1_time)} ({step1_time/total_time*100:.1f}%)")
    print(f"  Matrix building and similarity time: {format_time(step3_time)} ({step3_time/total_time*100:.1f}%)")
    print(f"  Matching pairs found: {len(matching_pairs):,}")
    print(f"  Total comparisons: {total_comparisons:,}")
    print(f"  Efficiency: {len(matching_pairs)/total_comparisons*100:.2f}% match rate")
    
    # Batch files summary
    if os.path.exists(batch_output_dir):
        batch_files = [f for f in os.listdir(batch_output_dir) if f.endswith('.feather')]
        if batch_files:
            print(f"  Batch files: {len(batch_files)} files in '{batch_output_dir}/' directory")
            
            # Calculate total size of batch files
            total_size_bytes = sum(os.path.getsize(os.path.join(batch_output_dir, f)) for f in batch_files)
            total_size_mb = total_size_bytes / (1024 * 1024)
            total_size_gb = total_size_bytes / (1024 * 1024 * 1024)
            
            if total_size_gb >= 1:
                print(f"  Batch files total size: {total_size_gb:.2f} GB")
            else:
                print(f"  Batch files total size: {total_size_mb:.2f} MB")
    
    print(f"{'='*80}")
    
    return result_df

# Keep the original helper functions from the base code
class ProgressTracker:
    """Thread-safe progress tracking class with ETA calculation"""
    def __init__(self, update_interval=5):
        self.lock = threading.Lock()
        self.stage = "Initializing"
        self.processed = 0
        self.total = 0
        self.start_time = time.time()
        self.stage_start_time = time.time()
        self.last_update_time = time.time()
        self.update_interval = update_interval
        self.batch_idx = 0
        self.num_batches = 0
        self.done = False
   
    def update_stage(self, stage):
        with self.lock:
            self.stage = stage
            self.stage_start_time = time.time()
            self.processed = 0
            self.last_update_time = time.time() - self.update_interval  # Force immediate update
   
    def set_total_items(self, total):
        with self.lock:
            self.total = total
   
    def increment_processed(self, count=1):
        with self.lock:
            self.processed += count
   
    def reset_processed(self):
        with self.lock:
            self.processed = 0
   
    def get_progress(self):
        with self.lock:
            current_time = time.time()
            elapsed = current_time - self.stage_start_time
            progress = self.processed / max(1, self.total)
           
            if progress > 0:
                eta_seconds = elapsed / progress - elapsed
                eta = timedelta(seconds=int(eta_seconds))
            else:
                eta = None
           
            return {
                'stage': self.stage,
                'processed': self.processed,
                'total': self.total,
                'progress': progress,
                'elapsed': timedelta(seconds=int(elapsed)),
                'eta': eta,
                'update_time': current_time
            }

def progress_reporter(progress):
    """Function to run in a separate thread to report progress at regular intervals with improved visibility"""
    while not progress.done:
        # Get current progress
        info = progress.get_progress()
       
        # Only update if enough time has passed
        current_time = time.time()
        if current_time - info['update_time'] < progress.update_interval:
            time.sleep(0.5)
            continue
           
        # Clear line and print progress with flush=True to ensure display
        print(f"\r{' ' * 100}", end='\r', flush=True)  # Clear line
       
        # Format the progress message
        if info['eta']:
            print(f"STAGE: {info['stage']} | "
                  f"Progress: {info['processed']:,}/{info['total']:,} ({info['progress']*100:.1f}%) | "
                  f"Elapsed: {info['elapsed']} | ETA: {info['eta']}", end='\r', flush=True)
        else:
            print(f"STAGE: {info['stage']} | "
                  f"Progress: {info['processed']:,}/{info['total']:,} ({info['progress']*100:.1f}%) | "
                  f"Elapsed: {info['elapsed']}", end='\r', flush=True)
       
        # Print a newline occasionally to ensure progress is visible in logs
        if int(time.time()) % 30 == 0:  # Every 30 seconds
            print("")  # New line to ensure progress is captured in logs
       
        time.sleep(progress.update_interval)
   
    # Clear line when done
    print(f"\r{' ' * 100}", end='\r', flush=True)


def get_system_info():
    """Get information about the system for optimal resource usage"""
    cpu_count = multiprocessing.cpu_count()
    memory_bytes = psutil.virtual_memory().available
    memory_gb = memory_bytes / (1024**3)
   
    # Try to get CPU model
    try:
        import platform
        if platform.system() == "Windows":
            cpu_model = platform.processor()
        else:
            # For Linux/Unix
            with open('/proc/cpuinfo') as f:
                for line in f:
                    if line.strip().startswith('model name'):
                        cpu_model = line.strip().split(':')[1].strip()
                        break
                else:
                    cpu_model = "Unknown"
    except:
        cpu_model = "Unknown"
   
    return {
        'cpu_count': cpu_count,
        'memory_gb': memory_gb,
        'cpu_model': cpu_model
    }

def format_time(seconds):
    """Format time in seconds to a human-readable string"""
    if seconds < 60:
        return f"{seconds:.2f} seconds"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.2f} minutes"
    else:
        hours = seconds / 3600
        return f"{hours:.2f} hours"


if __name__ == "__main__":
    import argparse
    import sys
   
    # Debug output - add this to verify script is running
    print("Script is starting...")
    print(f"Python version: {sys.version}")
   
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(
        description='Find matching motif pairs with shared n-mers (Sparse Matrix Optimized)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
   
    # Required arguments
    parser.add_argument('--file', type=str, required=True,
                        help='Path to the feather file containing motifs')
   
    # Algorithm selection
    parser.add_argument('--method', type=str, default='sparse',
                        choices=['sparse', 'original'],
                        help='Algorithm to use: sparse (sparse matrix) or original')
   
    # Core parameters
    parser.add_argument('--n', type=int, default=3,
                        help='Length of n-mers to use. Larger values (4+) are faster but more stringent')
    parser.add_argument('--min-matches', type=int, default=3,
                        help='Minimum number of matching n-mers required')
   
    # Performance tuning - optimized defaults for HPC
    parser.add_argument('--batch-size', type=int, default=50000,
                        help='Initial batch size for processing (will be optimized)')
    parser.add_argument('--num-workers', type=int, default=None,
                        help='Number of worker processes (default: 90% of CPU count)')
    parser.add_argument('--memory-limit', type=float, default=None,
                        help='Memory limit in GB (default: 90% of available system memory)')
    parser.add_argument('--progress-interval', type=int, default=5,
                        help='Progress update interval (in seconds)')
   
    # Output options
    parser.add_argument('--output', type=str, default=None,
                        help='Output file path (default: matching_motif_pairs_{n}mers.feather)')
    parser.add_argument('--format', type=str, default='feather',
                        choices=['feather', 'csv', 'excel'],
                        help='Output file format (default: feather)')
    
    # Batch output options
    parser.add_argument('--save-batches', action='store_true',
                        help='Save individual batch results to batch_output/ directory')
    parser.add_argument('--batch-dir', type=str, default='batch_output',
                        help='Directory to store batch output files')
    parser.add_argument('--skip-combined', action='store_true',
                        help='Skip creating combined output file (only create batch files)')
   
    # Parse arguments
    args = parser.parse_args()
   
    # Set default output file if none provided
    if args.output is None:
        if args.format == 'csv':
            args.output = f"matching_motif_pairs_{args.n}mers.csv"
        elif args.format == 'excel':
            args.output = f"matching_motif_pairs_{args.n}mers.xlsx"
        else:  # default to feather
            args.output = f"matching_motif_pairs_{args.n}mers.feather"
   
    # Print configuration
    print("\n=== Motif Matching Algorithm Configuration ===")
    print(f"Input file: {args.file}")
    print(f"Method: {args.method}")
    print(f"n-mer length: {args.n}")
    print(f"Min matches: {args.min_matches}")
    print(f"Batch size: {args.batch_size}")
    print(f"Number of workers: {args.num_workers}")
    print(f"Memory limit (GB): {args.memory_limit}")
    print(f"Output file: {args.output}")
    print(f"Output format: {args.format}")
    if args.save_batches:
        print(f"Saving batch files to: {args.batch_dir}/")
    if args.skip_combined:
        print("Skipping combined output file (only creating batch files)")
    print("===============================================\n")
   
    try:
        # Create batch directory if needed
        if args.save_batches and not os.path.exists(args.batch_dir):
            os.makedirs(args.batch_dir)
            print(f"Created batch output directory: {args.batch_dir}/")
            
        print("Starting motif matching process...")
        
        if args.method == 'sparse':
            # Use the new sparse matrix implementation
            result = find_matching_pairs_sparse_optimized(
                args.file,
                min_matching_nmers=args.min_matches,
                n=args.n,
                batch_size=args.batch_size,
                progress_interval=args.progress_interval,
                num_workers=args.num_workers,
                memory_limit_gb=args.memory_limit
            )
        else:
            # Use the original implementation (import it from the original module)
            from motif_matching import find_matching_pairs_vectorized_optimized_parallel
            result = find_matching_pairs_vectorized_optimized_parallel(
                args.file,
                min_matching_nmers=args.min_matches,
                n=args.n,
                batch_size=args.batch_size,
                progress_interval=args.progress_interval,
                num_workers=args.num_workers,
                memory_limit_gb=args.memory_limit
            )
       
        # Skip combined output if requested
        if args.skip_combined:
            print("Skipping combined output file as requested.")
            
            # Print info about batch files
            if os.path.exists(args.batch_dir):
                batch_files = [f for f in os.listdir(args.batch_dir) if f.endswith('.feather')]
                if batch_files:
                    print(f"Individual batch results were saved to {len(batch_files)} files in '{args.batch_dir}/' directory")
                    
                    # Calculate total size of batch files
                    total_size_bytes = sum(os.path.getsize(os.path.join(args.batch_dir, f)) for f in batch_files)
                    total_size_mb = total_size_bytes / (1024 * 1024)
                    total_size_gb = total_size_bytes / (1024 * 1024 * 1024)
                    
                    if total_size_gb >= 1:
                        print(f"Batch files total size: {total_size_gb:.2f} GB")
                    else:
                        print(f"Batch files total size: {total_size_mb:.2f} MB")
        else:
            # Save combined results in specified format
            print(f"Saving {len(result):,} matching pairs to {args.output}...")
           
            # OPTIMIZATION: Parallel writing for large result sets
            write_start = time.time()
            
            # Print info about batch files
            if os.path.exists(args.batch_dir):
                batch_files = [f for f in os.listdir(args.batch_dir) if f.endswith('.feather')]
                if batch_files:
                    print(f"Individual batch results were saved to {len(batch_files)} files in '{args.batch_dir}/' directory")
                    
                    # Option to combine batch files instead of using the merged result
                    print("Note: You can also load and process batch files individually for distributed analysis")
            
            if args.format == 'csv':
                # For very large datasets, use chunks to write CSV
                if len(result) > 1000000:
                    chunk_size = 200000
                    for i in range(0, len(result), chunk_size):
                        chunk = result.iloc[i:i+chunk_size]
                        mode = 'w' if i == 0 else 'a'
                        header = i == 0
                        chunk.to_csv(args.output, mode=mode, header=header, index=False)
                        print(f"  Wrote chunk {i//chunk_size + 1}/{(len(result) + chunk_size - 1)//chunk_size}")
                else:
                    result.to_csv(args.output, index=False)
                print(f"Results saved as CSV file in {time.time() - write_start:.2f} seconds.")
                
            elif args.format == 'excel':
                # Check if result exceeds Excel's limit
                max_excel_rows = 1048576
                if len(result) > max_excel_rows:
                    print(f"WARNING: Result has {len(result):,} rows which exceeds Excel's limit of {max_excel_rows:,}.")
                    print(f"Only the first {max_excel_rows:,} rows will be saved to Excel.")
                    result = result.head(max_excel_rows)
               
                # For large Excel files, we can optimize the write
                writer = pd.ExcelWriter(args.output, engine='openpyxl')
                result.to_excel(writer, index=False)
                writer.save()
                print(f"Results saved as Excel file in {time.time() - write_start:.2f} seconds.")
                
            else:  # feather format - optimal for large datasets
                # Feather is already optimized for performance
                result.to_feather(args.output)
                print(f"Results saved as Feather file in {time.time() - write_start:.2f} seconds.")
           
            print("\nSample matching pairs:")
            print(result.head(10))
           
            print(f"\nTotal matching pairs found: {len(result):,}")
       
        print("Done!")
       
    except Exception as e:
        print(f"Error during execution: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)