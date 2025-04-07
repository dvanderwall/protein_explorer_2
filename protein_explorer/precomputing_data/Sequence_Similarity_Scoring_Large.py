import numpy as np
from collections import defaultdict
import itertools
from Bio.SubsMat import MatrixInfo
import random
from tqdm import tqdm

class AminoAcidLSH:
    """
    Locality Sensitive Hashing for amino acid sequences using MinHash and BLOSUM-based similarities.
    """
    
    def __init__(self, n_hash_functions=100, n_bands=20, seed=42):
        """
        Initialize the LSH object.
        
        Args:
            n_hash_functions: Number of hash functions to use
            n_bands: Number of bands for LSH
            seed: Random seed for reproducibility
        """
        self.n_hash_functions = n_hash_functions
        self.n_bands = n_bands
        self.rows_per_band = n_hash_functions // n_bands
        self.seed = seed
        self.band_buckets = [defaultdict(set) for _ in range(n_bands)]
        
        # Load BLOSUM62 matrix
        self.blosum62 = MatrixInfo.blosum62
        
        # Define amino acids (including X for unknown)
        self.amino_acids = list("ACDEFGHIKLMNPQRSTVWYX")
        
        # Generate hash functions
        random.seed(seed)
        self.permutations = [list(range(len(self.amino_acids) * 9)) for _ in range(n_hash_functions)]
        for p in self.permutations:
            random.shuffle(p)
    
    def _seq_to_kmer_set(self, seq):
        """Convert sequence to a set of k-mers (here using 3-mers)"""
        # Replace any rare characters with 'X'
        seq = ''.join(aa if aa in self.amino_acids else 'X' for aa in seq)
        
        # Generate 3-mers
        kmers = set()
        for i in range(len(seq) - 2):
            kmer = seq[i:i+3]
            if 'X' not in kmer:  # Skip kmers with unknown residues
                kmers.add(kmer)
        return kmers
    
    def _minhash_signature(self, seq):
        """Generate MinHash signature for a sequence"""
        kmer_set = self._seq_to_kmer_set(seq)
        if not kmer_set:
            return [float('inf')] * self.n_hash_functions
        
        # Convert k-mers to integers
        kmer_indices = []
        for kmer in kmer_set:
            for i, aa in enumerate(kmer):
                idx = self.amino_acids.index(aa) * 9 + i
                kmer_indices.append(idx)
        
        # Apply MinHash
        signature = []
        for permutation in self.permutations:
            min_hash = float('inf')
            for idx in kmer_indices:
                hash_val = permutation[idx % len(permutation)]
                min_hash = min(min_hash, hash_val)
            signature.append(min_hash)
        
        return signature
    
    def index_sequences(self, sequences):
        """
        Index sequences for fast candidate pair generation.
        
        Args:
            sequences: List of amino acid sequences
        """
        print(f"Indexing {len(sequences)} sequences...")
        self.sequences = sequences
        self.signatures = []
        
        # Generate signatures for all sequences
        for seq in tqdm(sequences):
            sig = self._minhash_signature(seq)
            self.signatures.append(sig)
            
        # Populate band buckets
        for i, sig in enumerate(tqdm(self.signatures)):
            for band_idx in range(self.n_bands):
                # Get band slice from signature
                start = band_idx * self.rows_per_band
                end = start + self.rows_per_band
                band = tuple(sig[start:end])
                
                # Add to bucket
                self.band_buckets[band_idx][band].add(i)
    
    def get_candidate_pairs(self, min_band_matches=1):
        """
        Generate candidate pairs that are potentially similar.
        
        Args:
            min_band_matches: Minimum number of bands that must match
            
        Returns:
            List of (seq_idx1, seq_idx2) tuples for candidate pairs
        """
        print("Generating candidate pairs...")
        candidate_matches = defaultdict(int)
        
        # Count band matches for each potential pair
        for band_idx in range(self.n_bands):
            for bucket in self.band_buckets[band_idx].values():
                if len(bucket) > 1:
                    for a, b in itertools.combinations(bucket, 2):
                        candidate_matches[(min(a, b), max(a, b))] += 1
        
        # Filter candidates by minimum band matches
        candidate_pairs = []
        for pair, matches in candidate_matches.items():
            if matches >= min_band_matches:
                candidate_pairs.append(pair)
                
        print(f"Found {len(candidate_pairs)} candidate pairs")
        return candidate_pairs
    
    def calculate_blosum_similarity(self, seq1, seq2):
        """
        Calculate BLOSUM62 similarity score between two sequences.
        
        Args:
            seq1, seq2: Amino acid sequences of equal length
            
        Returns:
            BLOSUM62 similarity score
        """
        score = 0
        for i in range(min(len(seq1), len(seq2))):
            a, b = seq1[i], seq2[i]
            if a == 'X' or b == 'X':
                continue
            
            pair = (a, b) if a <= b else (b, a)
            if pair in self.blosum62:
                score += self.blosum62[pair]
        return score
    
    def filter_candidates_by_blosum(self, candidate_pairs, threshold=30):
        """
        Filter candidate pairs by BLOSUM62 similarity.
        
        Args:
            candidate_pairs: List of (seq_idx1, seq_idx2) tuples
            threshold: Minimum BLOSUM62 score to keep a pair
            
        Returns:
            List of (seq_idx1, seq_idx2) tuples that pass the threshold
        """
        print(f"Filtering {len(candidate_pairs)} candidate pairs by BLOSUM score...")
        filtered_pairs = []
        
        for idx1, idx2 in tqdm(candidate_pairs):
            seq1 = self.sequences[idx1]
            seq2 = self.sequences[idx2]
            score = self.calculate_blosum_similarity(seq1, seq2)
            
            if score >= threshold:
                filtered_pairs.append((idx1, idx2))
                
        print(f"Retained {len(filtered_pairs)} pairs after filtering")
        return filtered_pairs


def batch_process_large_dataset(sequences, site_ids=None, batch_size=100000, n_hash_functions=100, n_bands=20, blosum_threshold=30):
    """
    Process a large dataset in batches to find similar sequence pairs.
    
    Args:
        sequences: List of amino acid sequences
        site_ids: List of site IDs corresponding to the sequences (optional)
        batch_size: Size of each batch
        n_hash_functions: Number of hash functions for LSH
        n_bands: Number of bands for LSH
        blosum_threshold: Minimum BLOSUM score to keep a pair
        
    Returns:
        If site_ids provided: List of (site_id1, site_id2, seq1, seq2) tuples for similar pairs
        Otherwise: List of (seq_idx1, seq_idx2) tuples for similar sequence pairs
    """
    import time
    from datetime import datetime, timedelta
    
    all_pairs = []
    total_pairs_found = 0
    n_sequences = len(sequences)
    n_batches = (n_sequences + batch_size - 1) // batch_size
    
    print(f"========== STARTING BATCH PROCESSING ==========")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Processing {n_sequences:,} sequences in {n_batches} batches (batch size: {batch_size:,})")
    print(f"LSH parameters: {n_hash_functions} hash functions, {n_bands} bands")
    print(f"BLOSUM similarity threshold: {blosum_threshold}")
    print(f"===============================================")
    
    start_time = time.time()
    batch_comparison_count = 0
    total_batch_comparisons = (n_batches * (n_batches + 1)) // 2  # Total number of batch comparisons
    
    for i in range(n_batches):
        batch_start_time = time.time()
        print(f"\n[Batch {i+1}/{n_batches}] Processing sequences {i*batch_size+1:,}-{min((i+1)*batch_size, n_sequences):,}")
        
        start_i = i * batch_size
        end_i = min(start_i + batch_size, n_sequences)
        batch_i = sequences[start_i:end_i]
        batch_i_size = len(batch_i)
        
        print(f"[Batch {i+1}/{n_batches}] Loading {batch_i_size:,} sequences into LSH index")
        
        # When comparing batch i with itself
        lsh = AminoAcidLSH(n_hash_functions=n_hash_functions, n_bands=n_bands)
        lsh.index_sequences(batch_i)
        
        print(f"[Batch {i+1}/{n_batches}] Finding candidate pairs within batch")
        candidates = lsh.get_candidate_pairs()
        print(f"[Batch {i+1}/{n_batches}] Found {len(candidates):,} candidate pairs within batch")
        
        print(f"[Batch {i+1}/{n_batches}] Filtering candidates by BLOSUM similarity")
        filtered = lsh.filter_candidates_by_blosum(candidates, threshold=blosum_threshold)
        print(f"[Batch {i+1}/{n_batches}] Retained {len(filtered):,} pairs after BLOSUM filtering")
        
        # Adjust indices to global sequence indices
        adjusted_pairs = [(start_i + idx1, start_i + idx2) for idx1, idx2 in filtered]
        all_pairs.extend(adjusted_pairs)
        total_pairs_found += len(adjusted_pairs)
        batch_comparison_count += 1
        
        # Progress update
        elapsed_time = time.time() - start_time
        progress = batch_comparison_count / total_batch_comparisons
        estimated_total_time = elapsed_time / progress if progress > 0 else 0
        estimated_time_remaining = estimated_total_time - elapsed_time
        
        print(f"[Progress] {progress:.1%} complete, {batch_comparison_count}/{total_batch_comparisons} batch comparisons")
        print(f"[Progress] Elapsed time: {timedelta(seconds=int(elapsed_time))}")
        print(f"[Progress] Estimated time remaining: {timedelta(seconds=int(estimated_time_remaining))}")
        print(f"[Progress] Total pairs found so far: {total_pairs_found:,}")
        
        # Compare with other batches
        for j in range(i+1, n_batches):
            cross_batch_start_time = time.time()
            print(f"\n[Batch {i+1}×{j+1}/{n_batches}] Comparing batch {i+1} with batch {j+1}")
            
            start_j = j * batch_size
            end_j = min(start_j + batch_size, n_sequences)
            batch_j = sequences[start_j:end_j]
            batch_j_size = len(batch_j)
            
            print(f"[Batch {i+1}×{j+1}/{n_batches}] Processing cross-batch comparison: {batch_i_size:,} × {batch_j_size:,} sequences")
            
            # Create cross-batch LSH
            cross_batch_lsh = AminoAcidLSH(n_hash_functions=n_hash_functions, n_bands=n_bands)
            cross_batch_sequences = batch_i + batch_j
            print(f"[Batch {i+1}×{j+1}/{n_batches}] Loading {len(cross_batch_sequences):,} sequences into LSH index")
            cross_batch_lsh.index_sequences(cross_batch_sequences)
            
            # Get and filter candidates
            print(f"[Batch {i+1}×{j+1}/{n_batches}] Finding candidate pairs between batches")
            cross_candidates = cross_batch_lsh.get_candidate_pairs()
            print(f"[Batch {i+1}×{j+1}/{n_batches}] Found {len(cross_candidates):,} total candidate pairs")
            
            # Keep only pairs that cross the batch boundary
            cross_batch_candidates = []
            for idx1, idx2 in cross_candidates:
                if (idx1 < len(batch_i) and idx2 >= len(batch_i)) or (idx1 >= len(batch_i) and idx2 < len(batch_i)):
                    if idx1 >= len(batch_i):
                        idx1, idx2 = idx2, idx1  # Ensure idx1 is from batch_i
                    cross_batch_candidates.append((idx1, idx2 - len(batch_i)))
            
            print(f"[Batch {i+1}×{j+1}/{n_batches}] Extracted {len(cross_batch_candidates):,} cross-batch candidate pairs")
            
            # Filter by BLOSUM score
            print(f"[Batch {i+1}×{j+1}/{n_batches}] Filtering by BLOSUM similarity (threshold={blosum_threshold})")
            filtered_cross = []
            blosum_start_time = time.time()
            
            # Process in chunks and show progress
            chunk_size = max(1, len(cross_batch_candidates) // 10)  # Show progress in 10% increments
            for chunk_idx in range(0, len(cross_batch_candidates), chunk_size):
                chunk = cross_batch_candidates[chunk_idx:chunk_idx + chunk_size]
                chunk_pairs = []
                
                for idx1, idx2 in tqdm(chunk):
                    seq1 = batch_i[idx1]
                    seq2 = batch_j[idx2]
                    score = cross_batch_lsh.calculate_blosum_similarity(seq1, seq2)
                    
                    if score >= blosum_threshold:
                        chunk_pairs.append((start_i + idx1, start_j + idx2))
                
                filtered_cross.extend(chunk_pairs)
                
                # Progress update for BLOSUM filtering
                percent_done = min(100, int((chunk_idx + len(chunk)) / len(cross_batch_candidates) * 100))
                print(f"[Batch {i+1}×{j+1}/{n_batches}] BLOSUM filtering: {percent_done}% complete, found {len(filtered_cross):,} similar pairs so far")
            
            blosum_time = time.time() - blosum_start_time
            print(f"[Batch {i+1}×{j+1}/{n_batches}] BLOSUM filtering completed in {timedelta(seconds=int(blosum_time))}")
            print(f"[Batch {i+1}×{j+1}/{n_batches}] Retained {len(filtered_cross):,} similar pairs after BLOSUM filtering")
            
            all_pairs.extend(filtered_cross)
            total_pairs_found += len(filtered_cross)
            batch_comparison_count += 1
            
            cross_batch_time = time.time() - cross_batch_start_time
            print(f"[Batch {i+1}×{j+1}/{n_batches}] Cross-batch comparison completed in {timedelta(seconds=int(cross_batch_time))}")
            
            # Overall progress update
            elapsed_time = time.time() - start_time
            progress = batch_comparison_count / total_batch_comparisons
            estimated_total_time = elapsed_time / progress if progress > 0 else 0
            estimated_time_remaining = estimated_total_time - elapsed_time
            
            print(f"[Progress] {progress:.1%} complete, {batch_comparison_count}/{total_batch_comparisons} batch comparisons")
            print(f"[Progress] Elapsed time: {timedelta(seconds=int(elapsed_time))}")
            print(f"[Progress] Estimated time remaining: {timedelta(seconds=int(estimated_time_remaining))}")
            print(f"[Progress] Total pairs found so far: {total_pairs_found:,}")
        
        batch_time = time.time() - batch_start_time
        print(f"[Batch {i+1}/{n_batches}] Batch processing completed in {timedelta(seconds=int(batch_time))}")
    
    total_time = time.time() - start_time
    
    print(f"\n========== BATCH PROCESSING COMPLETED ==========")
    print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total processing time: {timedelta(seconds=int(total_time))}")
    print(f"Total similar pairs found: {total_pairs_found:,}")
    print(f"================================================")
    
    # If site_ids are provided, convert index pairs to site_id pairs with sequences
    if site_ids:
        print(f"Converting {len(all_pairs):,} index pairs to site_id pairs with sequences...")
        result_pairs = []
        for idx1, idx2 in all_pairs:
            result_pairs.append((
                site_ids[idx1],
                site_ids[idx2],
                sequences[idx1],
                sequences[idx2]
            ))
        return result_pairs
    else:
        return all_pairs


# Load and process the Excel file containing STY motifs
def main():
    import pandas as pd
    import os
    
    # Load the Excel file
    file_path = 'C:/Users/mz30/protein_explorer/protein_explorer/precomputing_data/All_STY_Motifs.xlsx'
    print(f"Loading data from {file_path}")
    
    # Check if file exists
    if not os.path.exists(file_path):
        print(f"Error: File {file_path} not found!")
        return
    
    try:
        # Read the Excel file - assuming first column is site ID and second column is motif
        df = pd.read_excel(file_path)
        
        # Output file structure information
        print(f"File loaded successfully. Columns found: {df.columns.tolist()}")
        print(f"Total number of sequences: {len(df)}")
        
        # Assuming the second column contains the motifs
        motif_column = df.columns[1]
        print(f"Using column '{motif_column}' for motif sequences")
        
        # Extract motif sequences
        sequences = df[motif_column].tolist()
        site_ids = df[df.columns[0]].tolist()
        
        # Basic data validation
        print(f"Sample sequences (first 5):")
        for i in range(min(5, len(sequences))):
            print(f"  {site_ids[i]}: {sequences[i]}")
        
        # Check sequence lengths
        seq_lengths = [len(seq) for seq in sequences]
        print(f"Sequence length statistics: min={min(seq_lengths)}, max={max(seq_lengths)}, most common={max(set(seq_lengths), key=seq_lengths.count)}")
        
        # Process data and find similar pairs
        similar_pairs = batch_process_large_dataset(
            sequences=sequences,
            site_ids=site_ids,  # Pass the site IDs to track them
            batch_size=50000,  # Adjust based on available memory
            n_hash_functions=100,
            n_bands=20,
            blosum_threshold=30
        )
        
        print(f"Found {len(similar_pairs)} similar sequence pairs")
        
        # Output some examples
        print("Sample of similar pairs found:")
        for i, (site_id1, site_id2, seq1, seq2) in enumerate(similar_pairs[:10]):
            print(f"Pair {i+1}: {site_id1} <-> {site_id2}")
            print(f"  {seq1} <-> {seq2}")
        
        # Save results to CSV
        output_path = os.path.join(os.path.dirname(file_path), "similar_sty_motif_pairs.csv")
        with open(output_path, 'w') as f:
            f.write("site_id_1,site_id_2,motif_1,motif_2\n")
            for site_id1, site_id2, seq1, seq2 in similar_pairs:
                f.write(f"{site_id1},{site_id2},{seq1},{seq2}\n")
        
        print(f"Results saved to {output_path}")
        
    except Exception as e:
        print(f"Error processing file: {str(e)}")


if __name__ == "__main__":
    main()