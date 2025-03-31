import h5py
import pandas as pd
import os

def inspect_h5_and_convert(h5_filepath, feather_filepath=None):
    """
    Print first 10 rows of HDF5 file, convert to Feather, and print first 10 rows of Feather
    
    Args:
        h5_filepath: Path to the HDF5 file
        feather_filepath: Path for the output Feather file (optional)
    """
    # Generate feather filepath if not provided
    if feather_filepath is None:
        feather_filepath = h5_filepath.replace('.h5', '.feather')
    
    print(f"\nInspecting H5 file: {h5_filepath}\n")
    
    # Open the H5 file and print first 10 rows
    with h5py.File(h5_filepath, 'r') as h5_file:
        # Check for the 'similarities' dataset
        if 'similarities' not in h5_file:
            print("Error: Could not find 'similarities' dataset in the H5 file")
            print("Available datasets:", list(h5_file.keys()))
            return None
        
        # Get the dataset
        dataset = h5_file['similarities']
        num_rows = dataset.shape[0]
        
        print(f"Dataset has {num_rows:,} rows")
        print("\nFirst 10 rows from HDF5 file:")
        
        # Read first 10 rows
        first_10_rows = dataset[:10]
        
        # Print as table - handle bytes by decoding to strings
        print("-" * 80)
        print(f"{'ID1':<40} | {'ID2':<40} | {'Similarity'}")
        print("-" * 80)
        
        for row in first_10_rows:
            # Convert bytes to strings if needed
            id1 = row['id1'].decode('utf-8') if isinstance(row['id1'], bytes) else row['id1']
            id2 = row['id2'].decode('utf-8') if isinstance(row['id2'], bytes) else row['id2']
            print(f"{id1:<40} | {id2:<40} | {row['similarity']:.6f}")
        
        print("-" * 80)
        
        # Convert to DataFrame and write to Feather
        print(f"\nConverting to Feather file: {feather_filepath}")
        
        # Read the entire dataset (or use chunked approach if it's very large)
        print("Reading data and converting to DataFrame...")
        
        # For very large files, use chunked approach
        chunk_size = 100000
        
        if num_rows > chunk_size:
            print(f"Large file detected. Processing in chunks of {chunk_size:,} rows...")
            
            # Initialize list to store dataframes
            df_list = []
            
            # Calculate number of chunks
            num_chunks = (num_rows + chunk_size - 1) // chunk_size
            
            # Process in chunks
            for i in range(num_chunks):
                start_idx = i * chunk_size
                end_idx = min((i + 1) * chunk_size, num_rows)
                
                print(f"Processing chunk {i+1}/{num_chunks} (rows {start_idx:,} to {end_idx:,})...")
                
                # Read chunk
                chunk_data = dataset[start_idx:end_idx]
                
                # Convert to DataFrame and handle bytes
                chunk_df = pd.DataFrame({
                    'id1': [x.decode('utf-8') if isinstance(x, bytes) else x for x in chunk_data['id1']],
                    'id2': [x.decode('utf-8') if isinstance(x, bytes) else x for x in chunk_data['id2']],
                    'similarity': chunk_data['similarity']
                })
                
                df_list.append(chunk_df)
            
            # Combine all chunks
            print("Combining chunks...")
            df = pd.concat(df_list, ignore_index=True)
        else:
            # For smaller files, read all at once
            data = dataset[:]
            df = pd.DataFrame({
                'id1': [x.decode('utf-8') if isinstance(x, bytes) else x for x in data['id1']],
                'id2': [x.decode('utf-8') if isinstance(x, bytes) else x for x in data['id2']],
                'similarity': data['similarity']
            })
        
        # Write to Feather
        print(f"Writing {len(df):,} rows to Feather file...")
        df.to_feather(feather_filepath)
        
    # Read and print first 10 rows from the Feather file
    print("\nReading first 10 rows from the converted Feather file:")
    feather_df = pd.read_feather(feather_filepath)
    
    print("\nFirst 10 rows from Feather file:")
    print("-" * 80)
    print(feather_df.head(10).to_string(index=False))
    print("-" * 80)
    
    print("\nConversion complete!")
    print(f"Original H5 file: {h5_filepath}")
    print(f"New Feather file: {feather_filepath}")
    print(f"H5 file size: {os.path.getsize(h5_filepath) / (1024 * 1024):.2f} MB")
    print(f"Feather file size: {os.path.getsize(feather_filepath) / (1024 * 1024):.2f} MB")
    
    return feather_filepath

# Set your file path
h5_filepath = "C:/Users/mz30/protein_explorer/phosphosite_similarity_results.h5"
feather_filepath = h5_filepath.replace('.h5', '.feather')

# Run the conversion and inspection
inspect_h5_and_convert(h5_filepath, feather_filepath)