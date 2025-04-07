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


def find_matching_pairs_vectorized_optimized_parallel(file_path, min_matching_nmers=3, n=3, batch_size=50000,
                                                    progress_interval=10, num_workers=None, memory_limit_gb=None, 
                                                    subbatch_size=20):
    """
    Highly optimized parallel implementation of motif matching for high-performance computing.
   
    Key optimizations:
    1. Multi-core processing with work division optimized for 118 CPUs
    2. Shared memory for large data structures
    3. Adaptive batch sizing based on 512GB memory
    4. Dynamic workload balancing for uneven computation loads
    5. Enhanced vectorization with numpy
    6. Memory-efficient data structures
   
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
    print(f"MOTIF MATCHING - HIGH PERFORMANCE COMPUTING IMPLEMENTATION")
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
    print(f"  Subbatch size for vectorization: {subbatch_size}")
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
   
    # STEP 1: Create optimized fingerprints and n-mer sets
    progress.update_stage("Creating fingerprints and n-mer sets")
    step1_start = time.time()
    print("Creating fingerprints and n-mer sets...")
   
    # OPTIMIZATION: Use parallel preprocessing for fingerprinting
    # Pre-allocate memory for better performance
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # 20 standard amino acids
    aa_to_idx = {aa: i for i, aa in enumerate(amino_acids)}
   
    # Function to process a chunk of data in parallel
    def process_chunk(chunk_df):
        chunk_valid_ids = []
        chunk_fingerprints = []
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
            
            # Create a boolean fingerprint
            fingerprint = np.zeros(len(amino_acids), dtype=np.bool_)
            
            for aa in motif:
                if aa in aa_to_idx:  # Skip 'X'
                    fingerprint[aa_to_idx[aa]] = True
            
            chunk_fingerprints.append(fingerprint)
            
            # Extract n-mers directly to sets
            nmers = set(motif[i:i+n] for i in range(len(motif) - n + 1) if 'X' not in motif[i:i+n])
            chunk_id_to_nmers[id_val] = nmers
            
        return chunk_valid_ids, chunk_fingerprints, chunk_id_to_motif, chunk_id_to_nmers, chunk_valid_count
    
    # OPTIMIZATION: Split dataframe into chunks and process in parallel
    # Calculate optimal chunk size based on CPU count
    chunk_size = max(1000, min(100000, total_rows // (num_workers * 2)))
    chunks = [df[i:i + chunk_size] for i in range(0, total_rows, chunk_size)]
    
    valid_ids = []
    fingerprints = []
    id_to_motif = {}
    id_to_nmers = {}
    valid_count = 0
    total_nmers = 0
    
    progress.set_total_items(len(chunks))
    
    # Process chunks in parallel
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(process_chunk, chunk): i for i, chunk in enumerate(chunks)}
        
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing chunks"):
            chunk_valid_ids, chunk_fingerprints, chunk_id_to_motif, chunk_id_to_nmers, chunk_valid_count = future.result()
            
            valid_ids.extend(chunk_valid_ids)
            fingerprints.extend(chunk_fingerprints)
            id_to_motif.update(chunk_id_to_motif)
            id_to_nmers.update(chunk_id_to_nmers)
            valid_count += chunk_valid_count
            
            # Count total n-mers
            for nmers in chunk_id_to_nmers.values():
                total_nmers += len(nmers)
            
            progress.increment_processed()
   
    # Convert fingerprints to numpy array for vectorized operations
    fingerprints = np.array(fingerprints, dtype=np.bool_)
   
    # Free memory
    del df
    del chunks
   
    step1_time = time.time() - step1_start
    print(f"Fingerprinting complete in {step1_time:.2f} seconds")
    print(f"Valid motifs: {valid_count:,}, Total n-mers: {total_nmers:,}")
   
    # STEP 2: Find matches using optimized parallel processing
    progress.update_stage("Finding matches (multi-core processing)")
    step2_start = time.time()
    print("\nStarting parallel matching process...")
   
    # Number of IDs
    num_ids = len(valid_ids)
   
    # OPTIMIZATION: Better memory estimation and utilization for high-memory systems
    # Calculate memory required per batch entry
    estimated_memory_per_id = estimate_memory_usage(fingerprints, id_to_nmers)
    print(f"Estimated memory per ID: {estimated_memory_per_id:.2f} MB")
   
    # Calculate optimal batch size based on available memory
    # Increased safety factor for high-memory systems
    optimal_batch_size = max(100, min(
        batch_size,
        int((memory_limit_gb * 1024) / (estimated_memory_per_id * num_workers) * 0.8)
    ))
   
    # OPTIMIZATION: Adjust batch size for fewer, larger batches on high-CPU systems
    comp_batch_size = min(optimal_batch_size, max(5000, num_ids // (num_workers * 4)))
    num_batches = (num_ids + comp_batch_size - 1) // comp_batch_size
   
    print(f"Processing {num_ids:,} IDs in {num_batches:,} batches of size {comp_batch_size:,}")
    print(f"Using {num_workers} worker processes")
   
    # OPTIMIZATION: Pre-compute min n-mer threshold to avoid doing this in the loop
    min_nmer_thresholds = {}
    for id_val, nmers in id_to_nmers.items():
        min_nmer_thresholds[id_val] = len(nmers) >= min_matching_nmers
   
    # Dictionary to store matching pairs (for final results)
    matching_pairs = {}
    matches_lock = threading.Lock()  # Lock for thread-safe updates to matching_pairs
   
    # Set up progress tracking for batches
    progress.set_total_items(num_batches)
    progress.reset_processed()

    # OPTIMIZATION: Implement dynamic batch sizing for work balancing
    # We'll create variable-sized batches based on the number of n-mers
    # to better balance the workload across processes
    
    # First, calculate the total workload
    id_to_workload = {}
    for i, id_val in enumerate(valid_ids):
        # Workload approximation: number of n-mers × position in the list
        # (later positions have more comparisons)
        workload = len(id_to_nmers[id_val]) * (num_ids - i)
        id_to_workload[i] = workload
    
    # Sort indices by workload (descending)
    sorted_indices = sorted(range(num_ids), key=lambda i: id_to_workload[i], reverse=True)
    
    # Create balanced batches
    balanced_batches = []
    batch_workloads = [0] * num_workers
    batch_assignments = [[] for _ in range(num_workers)]
    
    # Assign IDs to workers using a greedy approach
    for idx in sorted_indices:
        # Find the worker with the lowest current workload
        min_workload_worker = batch_workloads.index(min(batch_workloads))
        batch_assignments[min_workload_worker].append(idx)
        batch_workloads[min_workload_worker] += id_to_workload[idx]
    
    # Convert worker assignments to batches
    for worker_batch in batch_assignments:
        if worker_batch:  # Skip empty batches
            # Further split large worker batches
            if len(worker_batch) > comp_batch_size:
                for i in range(0, len(worker_batch), comp_batch_size):
                    sub_batch = worker_batch[i:i+comp_batch_size]
                    balanced_batches.append(sub_batch)
            else:
                balanced_batches.append(worker_batch)
    
    num_batches = len(balanced_batches)
    print(f"Created {num_batches} dynamically balanced work batches")
   
    # Create batch_output directory if it doesn't exist
    batch_output_dir = "batch_output"
    if not os.path.exists(batch_output_dir):
        os.makedirs(batch_output_dir)
        print(f"Created directory: {batch_output_dir}")
    
    # Process batches in parallel with improved batch distribution
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit jobs for each batch
        batch_jobs = []
        for batch_idx, id_indices in enumerate(balanced_batches):
            # Submit the batch job with indices
            job = executor.submit(
                process_batch_by_indices_vectorized,  # Using the new vectorized function
                batch_idx=batch_idx,
                num_batches=num_batches,
                valid_ids=valid_ids,
                id_indices=id_indices,  # Pass indices instead of start/end
                fingerprints=fingerprints,
                id_to_nmers=id_to_nmers,
                min_nmer_thresholds=min_nmer_thresholds,
                num_ids=num_ids,
                min_matching_nmers=min_matching_nmers,
                id_to_motif=id_to_motif,  # Pass motif data for batch file output
                batch_output_dir=batch_output_dir,  # Pass output directory
                n=n,  # Pass n-mer size for column naming
                subbatch_size=subbatch_size  # New parameter for controlling vectorization
            )
            batch_jobs.append(job)
    
        # Collect results as they complete
        total_comparisons = 0
        completed_batches = 0
        total_batch_time = 0

        # Custom tqdm with batch timing
        batch_pbar = tqdm(total=num_batches, desc="Processing batches",
                        mininterval=0.5, maxinterval=3.0, leave=True,
                        bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]',
                        ncols=100, position=0)

        for job in as_completed(batch_jobs):
            batch_results, batch_comparisons, batch_time = job.result()
        
            # Update timing statistics
            completed_batches += 1
            total_batch_time += batch_time
            avg_batch_time = total_batch_time / completed_batches
            remaining_batches = num_batches - completed_batches
            estimated_remaining_time = avg_batch_time * remaining_batches
        
            # Update global matching pairs with thread safety
            with matches_lock:
                matching_pairs.update(batch_results)
                total_comparisons += batch_comparisons
        
            # Update progress with timing information
            batch_pbar.set_description(f"Processing batches - Avg: {avg_batch_time:.2f}s/batch, ETA: {format_time(estimated_remaining_time)}")
            batch_pbar.update(1)
            progress.increment_processed()

        batch_pbar.close()

    # Update final progress
    progress.update_stage("Preparing final results")
   
    step2_time = time.time() - step2_start
    print(f"\nMatching complete in {format_time(step2_time)}")
    print(f"Found {len(matching_pairs):,} matching pairs from {total_comparisons:,} comparisons")
   
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
    print(f"  Fingerprinting time: {format_time(step1_time)} ({step1_time/total_time*100:.1f}%)")
    print(f"  Matching time: {format_time(step2_time)} ({step2_time/total_time*100:.1f}%)")
    print(f"  Matching pairs found: {len(matching_pairs):,}")
    print(f"  Total comparisons: {total_comparisons:,}")
    print(f"  Efficiency: {len(matching_pairs)/total_comparisons*100:.2f}% match rate")
    
    # Batch files summary
    batch_output_dir = "batch_output"
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


def process_batch_by_indices_vectorized(batch_idx, num_batches, valid_ids, id_indices, fingerprints, id_to_nmers,
                     min_nmer_thresholds, num_ids, min_matching_nmers, id_to_motif=None, 
                     batch_output_dir=None, n=3, subbatch_size=10):
    """
    Process a batch of comparisons using specific ID indices with subbatch vectorization.
    This enhanced version processes multiple IDs at once using vectorized operations.
    
    Args:
        batch_idx: Index of this batch
        num_batches: Total number of batches
        valid_ids: List of all valid IDs
        id_indices: Indices of IDs to process in this batch
        fingerprints: Numpy array of all fingerprints
        id_to_nmers: Dictionary mapping IDs to their n-mer sets
        min_nmer_thresholds: Dictionary with boolean flags for IDs with enough n-mers
        num_ids: Total number of IDs
        min_matching_nmers: Minimum number of n-mers required for a match
        id_to_motif: Dictionary mapping IDs to their motifs (for output)
        batch_output_dir: Directory to save batch results
        n: Size of n-mers
        subbatch_size: Number of IDs to process at once in a subbatch
        
    Returns:
        Tuple of (matching_pairs_dict, comparisons_count, processing_time)
    """
    batch_start_time = time.time()
    print(f"[BATCH {batch_idx}] Starting batch processing with {len(id_indices)} IDs using subbatch vectorization")
    
    # Dictionary to store matches found in this batch
    batch_matching_pairs = {}
    batch_comparisons = 0
    
    # Split the id_indices into subbatches for vectorized processing
    subbatches = [id_indices[i:i+subbatch_size] for i in range(0, len(id_indices), subbatch_size)]
    print(f"[BATCH {batch_idx}] Created {len(subbatches)} subbatches of size up to {subbatch_size}")
    
    # Process each subbatch 
    for subbatch_idx, subbatch in enumerate(subbatches):
        subbatch_start_time = time.time()
        
        # Print progress for subbatches
        if subbatch_idx % 5 == 0 or subbatch_idx == len(subbatches) - 1:
            elapsed = time.time() - batch_start_time
            print(f"[BATCH {batch_idx}] Processing subbatch {subbatch_idx+1}/{len(subbatches)} ({(subbatch_idx+1)/len(subbatches)*100:.1f}%) in {elapsed:.2f}s")
        
        # Skip IDs that don't have enough n-mers based on the threshold
        valid_subbatch = []
        for idx in subbatch:
            id_val = valid_ids[idx]
            if min_nmer_thresholds[id_val]:
                valid_subbatch.append(idx)
        
        if not valid_subbatch:
            continue
            
        # Extract fingerprints for all valid IDs in this subbatch
        subbatch_ids = [valid_ids[idx] for idx in valid_subbatch]
        
        # For each ID in the subbatch, process it against later IDs
        for position, idx in enumerate(valid_subbatch):
            id1 = valid_ids[idx]
            fingerprint1 = fingerprints[idx]
            nmers1 = id_to_nmers[id1]
            
            # Only compare with IDs that come after this one in the full list
            compare_indices = [i for i in range(idx + 1, num_ids)]
            
            if not compare_indices:
                continue
                
            # Get comparison fingerprints all at once
            compare_fingerprints = fingerprints[compare_indices]
            
            # VECTORIZED OPERATION: Find all shared amino acids at once
            # This applies fingerprint1 against all comparison fingerprints simultaneously
            matches = np.logical_and(compare_fingerprints, fingerprint1)
            shared_aas = np.sum(matches, axis=1)
            
            # Early filtering - only consider potential matches with enough shared amino acids
            potential_match_positions = np.where(shared_aas >= min_matching_nmers)[0]
            
            # For each potential match, check n-mer intersection
            match_count = 0
            for pos_idx, pos in enumerate(potential_match_positions):
                try:
                    compare_idx = compare_indices[pos]
                    batch_comparisons += 1
                    
                    id2 = valid_ids[compare_idx]
                    if not min_nmer_thresholds[id2]:
                        continue
                    
                    nmers2 = id_to_nmers[id2]
                    
                    # Use set intersection to find common n-mers
                    common_count = len(nmers1 & nmers2)
                    
                    if common_count >= min_matching_nmers:
                        batch_matching_pairs[(id1, id2)] = common_count
                        match_count += 1
                except Exception as e:
                    print(f"[BATCH {batch_idx}] ERROR processing match {id1}->{compare_idx}: {e}")
                    raise
            
            # Occasionally print progress for this ID in the subbatch
            if position % (max(1, subbatch_size // 5)) == 0:
                print(f"[BATCH {batch_idx}][SUBBATCH {subbatch_idx+1}] ID {id1}: Found {match_count} matches from {len(potential_match_positions)} candidates")
        
        subbatch_time = time.time() - subbatch_start_time
        if subbatch_idx % 5 == 0 or subbatch_idx == len(subbatches) - 1:
            print(f"[BATCH {batch_idx}] Subbatch {subbatch_idx+1} completed in {subbatch_time:.2f}s")
    
    # Calculate processing time for this batch
    batch_processing_time = time.time() - batch_start_time
    print(f"[BATCH {batch_idx}] Completed in {batch_processing_time:.2f}s with {batch_comparisons} comparisons and {len(batch_matching_pairs)} matches")
    
    # Write batch results to feather file if requested
    if batch_output_dir and id_to_motif:
        print(f"[BATCH {batch_idx}] Starting to write results to file")
        try:
            # Convert batch results to DataFrame
            batch_result_data = []
            for (id1, id2), count in batch_matching_pairs.items():
                motif1 = id_to_motif.get(id1, "")
                motif2 = id_to_motif.get(id2, "")
                batch_result_data.append((id1, id2, motif1, motif2, count))
            
            print(f"[BATCH {batch_idx}] Converted {len(batch_matching_pairs)} pairs to DataFrame format")
            
            if batch_result_data:  # Only write if we have results
                batch_df = pd.DataFrame(
                    batch_result_data,
                    columns=['ID1', 'ID2', 'Motif1', 'Motif2', f'Matching_{n}mer_Count']
                )
                
                # Sort by count descending
                batch_df = batch_df.sort_values(f'Matching_{n}mer_Count', ascending=False).reset_index(drop=True)
                print(f"[BATCH {batch_idx}] Sorted DataFrame with {len(batch_df)} rows")
                
                # Write to feather
                batch_file = os.path.join(batch_output_dir, f'batch_{batch_idx:04d}.feather')
                print(f"[BATCH {batch_idx}] Writing to {batch_file}")
                batch_df.to_feather(batch_file)
                
                # Log batch file creation
                print(f"[BATCH {batch_idx}] Successfully wrote {len(batch_df):,} pairs to {batch_file}")
        except Exception as e:
            print(f"[BATCH {batch_idx}] ERROR writing batch results: {e}")
            raise
    
    # Return results for this batch
    print(f"[BATCH {batch_idx}] Returning results ({len(batch_matching_pairs)} pairs)")
    return batch_matching_pairs, batch_comparisons, batch_processing_time




def process_batch_by_indices(batch_idx, num_batches, valid_ids, id_indices, fingerprints, id_to_nmers,
                     min_nmer_thresholds, num_ids, min_matching_nmers, id_to_motif=None, 
                     batch_output_dir=None, n=3):
    """
    Process a batch of comparisons using specific ID indices.
    This improved version takes a list of indices rather than a range.
    Also writes batch results to individual feather files.
    
    Returns:
        Tuple of (matching_pairs_dict, comparisons_count, processing_time)
    """
    batch_start_time = time.time()
    print(f"[BATCH {batch_idx}] Starting batch processing with {len(id_indices)} IDs")
    
    # Dictionary to store matches found in this batch
    batch_matching_pairs = {}
    batch_comparisons = 0
    
    # Process each ID in the batch
    for idx_pos, idx in enumerate(id_indices):
        # Print progress every 100 indices or so
        if idx_pos % 100 == 0:
            elapsed = time.time() - batch_start_time
            print(f"[BATCH {batch_idx}] Processed {idx_pos}/{len(id_indices)} indices ({idx_pos/len(id_indices)*100:.1f}%) in {elapsed:.2f}s")
        
        id1 = valid_ids[idx]
        if not min_nmer_thresholds[id1]:
            continue
        
        fingerprint1 = fingerprints[idx]
        nmers1 = id_to_nmers[id1]
        
        # Only compare with IDs that come after this one in the full list
        # This ensures we don't do duplicate comparisons
        compare_indices = [i for i in range(idx + 1, num_ids)]
        
        print(f"[BATCH {batch_idx}] ID {id1}: Comparing with {len(compare_indices)} potential matches")
        
        if not compare_indices:
            continue
        
        # Get compare fingerprints
        try:
            compare_fingerprints = fingerprints[compare_indices]
            print(f"[BATCH {batch_idx}] ID {id1}: Got comparison fingerprints shape {compare_fingerprints.shape}")
        except Exception as e:
            print(f"[BATCH {batch_idx}] ERROR getting fingerprints: {e}")
            raise
        
        # OPTIMIZATION: Vectorized fingerprint matching
        # Use logical AND to find shared amino acids
        try:
            matches = np.logical_and(compare_fingerprints, fingerprint1)
            print(f"[BATCH {batch_idx}] ID {id1}: Performed logical AND operation")
            
            # Count shared amino acids
            shared_aas = np.sum(matches, axis=1)
            print(f"[BATCH {batch_idx}] ID {id1}: Counted shared amino acids")
            
            # OPTIMIZATION: Early filtering - only consider potential matches
            # that have enough shared amino acids to potentially meet the criteria
            potential_match_positions = np.where(shared_aas >= min_matching_nmers)[0]
            print(f"[BATCH {batch_idx}] ID {id1}: Found {len(potential_match_positions)} potential matches after filtering")
        except Exception as e:
            print(f"[BATCH {batch_idx}] ERROR in vectorized matching: {e}")
            raise
        
        # For each potential match, check n-mer intersection
        match_count = 0
        for pos_idx, pos in enumerate(potential_match_positions):
            # Print progress on potential matches occasionally
            if len(potential_match_positions) > 100 and pos_idx % 100 == 0:
                print(f"[BATCH {batch_idx}] ID {id1}: Processed {pos_idx}/{len(potential_match_positions)} potential matches")
                
            try:
                compare_idx = compare_indices[pos]
                batch_comparisons += 1
                
                id2 = valid_ids[compare_idx]
                if not min_nmer_thresholds[id2]:
                    continue
                
                nmers2 = id_to_nmers[id2]
                
                # OPTIMIZATION: Use & operator for set intersection (faster than method call)
                common_count = len(nmers1 & nmers2)
                
                if common_count >= min_matching_nmers:
                    batch_matching_pairs[(id1, id2)] = common_count
                    match_count += 1
            except Exception as e:
                print(f"[BATCH {batch_idx}] ERROR processing match {id1}->{compare_idx}: {e}")
                raise
        
        print(f"[BATCH {batch_idx}] ID {id1}: Found {match_count} matches")
    
    # Calculate processing time for this batch
    batch_processing_time = time.time() - batch_start_time
    print(f"[BATCH {batch_idx}] Completed in {batch_processing_time:.2f}s with {batch_comparisons} comparisons and {len(batch_matching_pairs)} matches")
    
    # Write batch results to feather file if requested
    if batch_output_dir and id_to_motif:
        print(f"[BATCH {batch_idx}] Starting to write results to file")
        try:
            # Convert batch results to DataFrame
            batch_result_data = []
            for (id1, id2), count in batch_matching_pairs.items():
                motif1 = id_to_motif.get(id1, "")
                motif2 = id_to_motif.get(id2, "")
                batch_result_data.append((id1, id2, motif1, motif2, count))
            
            print(f"[BATCH {batch_idx}] Converted {len(batch_matching_pairs)} pairs to DataFrame format")
            
            if batch_result_data:  # Only write if we have results
                batch_df = pd.DataFrame(
                    batch_result_data,
                    columns=['ID1', 'ID2', 'Motif1', 'Motif2', f'Matching_{n}mer_Count']
                )
                
                # Sort by count descending
                batch_df = batch_df.sort_values(f'Matching_{n}mer_Count', ascending=False).reset_index(drop=True)
                print(f"[BATCH {batch_idx}] Sorted DataFrame with {len(batch_df)} rows")
                
                # Write to feather
                batch_file = os.path.join(batch_output_dir, f'batch_{batch_idx:04d}.feather')
                print(f"[BATCH {batch_idx}] Writing to {batch_file}")
                batch_df.to_feather(batch_file)
                
                # Log batch file creation
                print(f"[BATCH {batch_idx}] Successfully wrote {len(batch_df):,} pairs to {batch_file}")
        except Exception as e:
            print(f"[BATCH {batch_idx}] ERROR writing batch results: {e}")
            raise
    
    # Return results for this batch
    print(f"[BATCH {batch_idx}] Returning results ({len(batch_matching_pairs)} pairs)")
    return batch_matching_pairs, batch_comparisons, batch_processing_time

   
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

def estimate_memory_usage(fingerprints, id_to_nmers, sample_size=1000):
    """Estimate memory usage per ID"""
    # Estimate memory for fingerprints (numpy arrays)
    fingerprint_size = 0
    if len(fingerprints) > 0:
        fingerprint_size = fingerprints[0].nbytes
   
    # Estimate memory for n-mer sets (sample a few)
    nmer_set_sizes = []
    sampled_ids = list(id_to_nmers.keys())[:min(sample_size, len(id_to_nmers))]
    for id_val in sampled_ids:
        # Approximate memory usage of a set - each string is ~40 bytes overhead + actual string size
        nmers = id_to_nmers[id_val]
        if nmers:
            avg_nmer_len = sum(len(nmer) for nmer in nmers) / len(nmers)
            nmer_size = 40 + avg_nmer_len  # Bytes per n-mer
            nmer_set_sizes.append(nmer_size * len(nmers))
   
    avg_nmer_set_size = sum(nmer_set_sizes) / max(1, len(nmer_set_sizes))
   
    # Total memory per ID in MB
    total_memory_per_id = (fingerprint_size + avg_nmer_set_size) / (1024**2)
   
    return total_memory_per_id

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
        description='Find matching motif pairs with shared n-mers (High Performance Computing Optimized)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
   
    # Required arguments
    parser.add_argument('--file', type=str, required=True,
                        help='Path to the feather file containing motifs')
   
    # Algorithm selection - set default to 'parallel' since you want to use the parallel version
    parser.add_argument('--method', type=str, default='parallel',
                        choices=['parallel'],
                        help='Algorithm to use: parallel (high-performance parallel)')
   
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
    
    parser.add_argument('--subbatch-size', type=int, default=20,
                    help='Size of subbatches for vectorized processing (default: 20)')
   
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
        result = find_matching_pairs_vectorized_optimized_parallel(
            args.file,
            min_matching_nmers=args.min_matches,
            n=args.n,
            batch_size=args.batch_size,
            progress_interval=args.progress_interval,
            num_workers=args.num_workers,
            memory_limit_gb=args.memory_limit,
            subbatch_size=args.subbatch_size  # Just pass subbatch_size
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