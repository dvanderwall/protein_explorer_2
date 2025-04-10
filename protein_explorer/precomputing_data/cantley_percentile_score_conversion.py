import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.feather as feather
import time
import os
import uuid
import gc
import psutil

# Helper function to display memory usage
def get_memory_usage():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024  # MB

# Function to convert a full column to percentiles using ECDF
def column_to_percentile(series):
    if pd.api.types.is_numeric_dtype(series):
        # Create ECDF function using numpy for better performance
        sorted_data = np.sort(series.values)
        n = len(sorted_data)
        ecdf_y = np.arange(1, n + 1) / n
        
        # Use numpy's searchsorted for faster lookups
        idx = np.searchsorted(sorted_data, series.values)
        # Handle edge cases where value equals max
        idx[idx == n] = n - 1
        return pd.Series(ecdf_y[idx], index=series.index)
    else:
        # Return non-numeric columns unchanged
        return series

# Create a temp folder
temp_folder = f"temp_percentile_{uuid.uuid4().hex[:8]}"
os.makedirs(temp_folder, exist_ok=True)
print(f"Created temporary folder: {temp_folder}")

# File paths
input_file = "F:/Kinome/Tyr_motif_scores.csv"
output_file = "F:/Kinome/Tyrosine_motif_scores_percentile.feather"

# First, get column names and number of rows
print("Getting column information...")
with open(input_file, 'r') as f:
    header = f.readline().strip().split(',')

# Count rows (more memory efficient than pd.read_csv)
def count_rows(file_path):
    count = 0
    with open(file_path, 'r') as f:
        # Skip header
        next(f)
        # Count lines
        for _ in f:
            count += 1
    return count

print("Counting rows (this might take a while)...")
row_count = count_rows(input_file)
print(f"File has {row_count:,} rows and {len(header)} columns")

# Create an empty result dataframe with just the first column (ID column)
print("Reading ID column...")
id_col = pd.read_csv(input_file, usecols=[0])
result_df = id_col.copy()
print(f"Loaded ID column with {len(result_df):,} rows")

# Process columns in batches of 20
start_time = time.time()
cols_to_process = header[1:]  # Skip ID column
total_cols = len(cols_to_process)
batch_size = 40

print(f"Processing {total_cols} columns in batches of {batch_size}...")

# Calculate number of batches
num_batches = (total_cols + batch_size - 1) // batch_size

for batch_idx in range(num_batches):
    batch_start = time.time()
    start_col = batch_idx * batch_size
    end_col = min(start_col + batch_size, total_cols)
    
    batch_cols = cols_to_process[start_col:end_col]
    print(f"\nBatch {batch_idx+1}/{num_batches}, processing columns {start_col+1}-{end_col} of {total_cols}")
    
    # Read this batch of columns
    cols_to_read = [header[0]] + batch_cols  # Include ID column
    print(f"Reading {len(batch_cols)} columns...")
    try:
        batch_data = pd.read_csv(input_file, usecols=cols_to_read)
        print(f"Loaded batch with {len(batch_data):,} rows and {len(batch_cols)} columns")
        
        # Process each column in the batch
        for col_name in batch_cols:
            col_start = time.time()
            print(f"  Processing column: '{col_name}'", end="", flush=True)
            
            # Check if numeric
            if pd.api.types.is_numeric_dtype(batch_data[col_name]):
                print(f" (numeric)", end="", flush=True)
                # Calculate percentiles
                result_df[col_name] = column_to_percentile(batch_data[col_name])
            else:
                print(f" (non-numeric)", end="", flush=True)
                # Copy as-is
                result_df[col_name] = batch_data[col_name]
            
            col_end = time.time()
            print(f" - {col_end - col_start:.2f}s")
        
        # Clean up batch data
        del batch_data
        gc.collect()
        print(f"  Memory usage after batch: {get_memory_usage():.2f} MB")
        
    except Exception as e:
        print(f"ERROR processing batch: {str(e)}")
        # Save what we have so far
        error_file = f"{temp_folder}/error_at_batch_{batch_idx+1}.feather"
        feather.write_feather(result_df, error_file)
        print(f"Saved progress to {error_file}")
        continue
    
    batch_end = time.time()
    batch_time = batch_end - batch_start
    remaining_batches = num_batches - (batch_idx + 1)
    est_remaining_time = batch_time * remaining_batches
    
    print(f"Batch {batch_idx+1} completed in {batch_time:.2f}s")
    print(f"Estimated time remaining: {est_remaining_time/60:.1f} minutes")
    
    # Save checkpoint after each batch
    checkpoint_file = f"{temp_folder}/batch_checkpoint_{batch_idx+1}.feather"
    feather.write_feather(result_df, checkpoint_file)
    print(f"Saved checkpoint after batch {batch_idx+1}")

# Save final result
print(f"\nSaving final result to {output_file}...")
feather.write_feather(result_df, output_file)

# Clean up temporary files
print("Cleaning up temporary files...")
for file in os.listdir(temp_folder):
    os.remove(os.path.join(temp_folder, file))
os.rmdir(temp_folder)

end_time = time.time()
print(f"\nAll processing completed in {(end_time - start_time)/60:.2f} minutes")
print(f"Result saved to {output_file}")