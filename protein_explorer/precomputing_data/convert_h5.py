import h5py
import pandas as pd
import os

def convert_h5_to_feather(h5_filepath, feather_filepath=None):
    """
    Convert HDF5 file to Feather format using Python
    
    Args:
        h5_filepath: Path to the HDF5 file
        feather_filepath: Path for the output Feather file (optional)
    """
    # Generate feather filepath if not provided
    if feather_filepath is None:
        feather_filepath = h5_filepath.replace('.h5', '.feather')
    
    print(f"Opening H5 file: {h5_filepath}")
    
    # Open the H5 file
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
        
        # Read entire dataset at once if it's not too large
        # Or process in chunks if needed
        data = dataset[:]
        
        # Convert to pandas DataFrame
        print("Converting to DataFrame...")
        df = pd.DataFrame({
            'id1': data['id1'],
            'id2': data['id2'],
            'similarity': data['similarity']
        })
        
        # Write to Feather
        print(f"Writing to Feather file: {feather_filepath}")
        df.to_feather(feather_filepath)
    
    print("Conversion complete!")
    print(f"Original H5 file: {h5_filepath}")
    print(f"New Feather file: {feather_filepath}")
    
    return feather_filepath

# Example usage - save this as convert_h5.py and run it
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python convert_h5.py <h5_filepath> [feather_filepath]")
        sys.exit(1)
    
    h5_filepath = sys.argv[1]
    feather_filepath = sys.argv[2] if len(sys.argv) > 2 else None
    
    convert_h5_to_feather(h5_filepath, feather_filepath)