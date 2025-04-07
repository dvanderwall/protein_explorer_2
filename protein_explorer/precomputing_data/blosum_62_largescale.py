import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sqlite3
import tempfile
import time

def analyze_blosum_with_sqlite(feather_file, db_file=None, sample_rate=0.01):
    """
    Analyze blosum scores using SQLite to handle large files.
    
    Args:
        feather_file: Path to the feather file
        db_file: Path to SQLite database (temporary if None)
        sample_rate: Fraction of rows to sample (0.01 = 1%)
    """

    # Clear the existing empty database
    db_path = os.path.join(tempfile.gettempdir(), "blosum_analysis.db")
    conn = sqlite3.connect(db_path)
    conn.execute("DROP TABLE IF EXISTS blosum_scores")
    conn.commit()
    conn.close()
    print(f"Cleared empty table from {db_path}")

    print(f"Processing file: {feather_file}")
    
    # Create temporary database if needed
    if db_file is None:
        db_file = os.path.join(tempfile.gettempdir(), "blosum_analysis.db")
    
    print(f"Using SQLite database: {db_file}")
    
    # Connect to the database
    conn = sqlite3.connect(db_file)
    
    # Import data in chunks
    chunk_size = 1000
    table_exists = False
    
    try:
        # Check if table already exists
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='blosum_scores'")
        table_exists = cursor.fetchone() is not None
        
        if not table_exists:
            print(f"Importing data with {sample_rate*100:.2f}% sampling rate...")
            
            # Create table
            conn.execute("CREATE TABLE blosum_scores (score REAL)")
            
            # Process the feather file in chunks
            chunk_start = 0
            total_imported = 0
            
            while True:
                print(f"Reading chunk starting at row {chunk_start}...")
                try:
                    # Read a chunk with just the blosum score column
                    chunk = pd.read_feather(
                        feather_file, 
                        columns=['blosum62_score']
                    ).iloc[chunk_start:chunk_start + chunk_size]
                    
                    # If empty, we've reached the end
                    if len(chunk) == 0:
                        break
                    
                    # Sample the data
                    if sample_rate < 1.0:
                        chunk = chunk.sample(frac=sample_rate)
                    
                    # Insert into database
                    chunk.to_sql('blosum_scores', conn, if_exists='append', index=False)
                    
                    # Update counters
                    total_imported += len(chunk)
                    chunk_start += chunk_size
                    
                    print(f"Imported {total_imported} rows so far")
                    
                except Exception as e:
                    print(f"Error reading chunk: {str(e)}")
                    # Try to continue with next chunk
                    chunk_start += chunk_size
            
            print(f"Finished importing data. Total rows in database: {total_imported}")
        else:
            print("Using existing database table")
            
        # Create index for better performance
        print("Creating index...")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_score ON blosum_scores(score)")
        
        # Calculate basic statistics
        print("Calculating statistics...")
        cursor = conn.cursor()
        
        # Count
        cursor.execute("SELECT COUNT(*) FROM blosum_scores")
        count = cursor.fetchone()[0]
        
        # Mean
        cursor.execute("SELECT AVG(score) FROM blosum_scores")
        mean = cursor.fetchone()[0]
        
        # Min/Max
        cursor.execute("SELECT MIN(score), MAX(score) FROM blosum_scores")
        min_score, max_score = cursor.fetchone()
        
        # Percentiles
        percentiles = []
        for p in [10, 25, 50, 75, 90, 95, 99]:
            cursor.execute(f"""
                SELECT score FROM blosum_scores
                ORDER BY score
                LIMIT 1
                OFFSET (SELECT COUNT(*) * {p/100.0} FROM blosum_scores)
            """)
            percentiles.append((p, cursor.fetchone()[0]))
        
        # Calculate histogram
        bins = 100
        bin_width = (max_score - min_score) / bins
        
        histogram = []
        for i in range(bins):
            bin_start = min_score + i * bin_width
            bin_end = bin_start + bin_width
            
            cursor.execute(f"""
                SELECT COUNT(*) FROM blosum_scores
                WHERE score >= {bin_start} AND score < {bin_end}
            """)
            
            bin_count = cursor.fetchone()[0]
            histogram.append((bin_start, bin_end, bin_count))
        
        # Print results
        print("\n----- BLOSUM62 Score Analysis -----")
        print(f"Total rows analyzed: {count:,}")
        print(f"Score range: {min_score:.4f} to {max_score:.4f}")
        print(f"Mean score: {mean:.4f}")
        
        print("\nPercentiles:")
        for p, v in percentiles:
            print(f"{p}th percentile: {v:.4f}")
        
        print("\nHistogram:")
        max_count = max(h[2] for h in histogram)
        for bin_start, bin_end, bin_count in histogram:
            if bin_count > 0:
                bar_len = int(50 * bin_count / max_count)
                print(f"{bin_start:.2f}-{bin_end:.2f}: {'#' * bar_len} ({bin_count:,})")
        
        # Generate a simple plot if matplotlib is available
        try:
            bin_centers = [h[0] + (h[1]-h[0])/2 for h in histogram]
            bin_counts = [h[2] for h in histogram]
            
            plt.figure(figsize=(12, 6))
            plt.bar(bin_centers, bin_counts, width=bin_width*0.9)
            plt.xlabel('BLOSUM62 Score')
            plt.ylabel('Frequency')
            plt.title('Distribution of BLOSUM62 Scores')
            plt.grid(axis='y', alpha=0.3)
            
            # Add mean line
            plt.axvline(mean, color='r', linestyle='dashed', linewidth=1, label=f'Mean: {mean:.3f}')
            
            # Add median line (50th percentile)
            median = next(v for p, v in percentiles if p == 50)
            plt.axvline(median, color='g', linestyle='dashed', linewidth=1, label=f'Median: {median:.3f}')
            
            plt.legend()
            
            # Save the plot
            output_path = os.path.join(os.path.dirname(feather_file), "blosum_score_histogram.png")
            plt.savefig(output_path, dpi=300)
            print(f"\nSaved histogram to: {output_path}")
            
        except Exception as plot_error:
            print(f"Could not generate plot: {str(plot_error)}")
            print("Text-based histogram shown instead.")
        
    finally:
        # Close the connection
        conn.close()
        print(f"Analysis complete. Database file: {db_file}")

# Example usage
if __name__ == "__main__":
    file_path = "F:/Kinome/combined_STY_motif_matches_with_blosum.feather"
    # Use a sample rate of 1% to process the data faster
    analyze_blosum_with_sqlite(file_path, sample_rate=0.01)