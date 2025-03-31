import pymysql
import os
import logging
import time
import pandas as pd
import pyarrow.feather as feather
from tqdm import tqdm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Database connection parameters
DB_HOST = os.environ.get('DB_HOST', '35.245.113.195')
DB_USER = os.environ.get('DB_USER', 'root')
DB_PASS = os.environ.get('DB_PASS', '@Bismark6')
DB_NAME = os.environ.get('DB_NAME', 'kinoplex-db')
DB_PORT = int(os.environ.get('DB_PORT', '3306'))

def create_similarity_table_from_feather(feather_file_path):
    """
    Create a new Sequence_Similarity_Edges_4 table from feather file.
    
    Args:
        feather_file_path: Path to feather file containing similarity data
    """
    connection = None
    try:
        # Load the feather file
        logger.info(f"Loading feather file from {feather_file_path}...")
        df = feather.read_feather(feather_file_path)
        logger.info(f"Loaded {len(df)} rows of similarity data")
        
        # Print the actual column names from the feather file
        logger.info(f"Feather file contains columns: {df.columns.tolist()}")
        
        # Rename columns to match the database table
        # Assuming the feather file has columns that need to be mapped to ID1, ID2, Similarity
        # Replace these with the actual column names from your feather file
        column_mapping = {
            df.columns[0]: 'ID1',   # First column maps to ID1
            df.columns[1]: 'ID2',   # Second column maps to ID2
            df.columns[2]: 'Similarity'  # Third column maps to Similarity
        }
        
        logger.info(f"Mapping columns: {column_mapping}")
        df = df.rename(columns=column_mapping)
        
        # Connect to the database
        logger.info(f"Connecting to database {DB_NAME} at {DB_HOST}...")
        connection = pymysql.connect(
            host=DB_HOST,
            user=DB_USER,
            password=DB_PASS,
            database=DB_NAME,
            port=DB_PORT,
            charset='utf8mb4'
        )
        
        with connection.cursor() as cursor:
            # Check if the table already exists
            cursor.execute("SHOW TABLES LIKE 'Sequence_Similarity_Edges_4'")
            if cursor.fetchone():
                logger.warning("Table Sequence_Similarity_Edges_4 already exists")
                response = input("Do you want to drop and recreate it? (y/n): ")
                if response.lower() != 'y':
                    logger.info("Operation cancelled")
                    return
                cursor.execute("DROP TABLE Sequence_Similarity_Edges_4")
                logger.info("Dropped existing table")
            
            # Create the new table
            cursor.execute("""
                CREATE TABLE Sequence_Similarity_Edges_4 (
                    ID1 VARCHAR(100) NOT NULL,
                    ID2 VARCHAR(100) NOT NULL,
                    Similarity FLOAT NOT NULL
                )
            """)
            logger.info("Created new table Sequence_Similarity_Edges_4")
            
            # Prepare data for insertion
            logger.info("Preparing data for insertion...")
            # Convert DataFrame to list of tuples for batch insertion
            data_tuples = list(zip(df['ID1'], df['ID2'], df['Similarity']))
            
            # Insert data in batches to avoid memory issues
            batch_size = 10000
            total_batches = len(data_tuples) // batch_size + (1 if len(data_tuples) % batch_size > 0 else 0)
            
            logger.info(f"Inserting data in {total_batches} batches...")
            for i in tqdm(range(0, len(data_tuples), batch_size)):
                batch = data_tuples[i:i + batch_size]
                # Use executemany for efficient batch insertion
                cursor.executemany(
                    "INSERT INTO Sequence_Similarity_Edges_4 (ID1, ID2, Similarity) VALUES (%s, %s, %s)",
                    [(str(t[0]), str(t[1]), float(t[2])) for t in batch]
                )
                connection.commit()  # Commit each batch
            
            logger.info("All data inserted successfully")
            
            # Create the same indexes as the original table
            logger.info("Creating indexes...")
            start_time = time.time()
            cursor.execute("""
                CREATE INDEX idx_id1_sim ON Sequence_Similarity_Edges_4(ID1(50), Similarity)
            """)
            logger.info(f"Created idx_id1_sim index in {time.time() - start_time:.2f} seconds")
            
            start_time = time.time()
            cursor.execute("""
                CREATE INDEX idx_id2_sim ON Sequence_Similarity_Edges_4(ID2(50), Similarity)
            """)
            logger.info(f"Created idx_id2_sim index in {time.time() - start_time:.2f} seconds")
            
            # Final commit
            connection.commit()
            logger.info("Successfully created Sequence_Similarity_Edges_4 with all data and indexes")
            
            # Print table summary
            cursor.execute("SELECT COUNT(*) FROM Sequence_Similarity_Edges_4")
            count = cursor.fetchone()[0]
            logger.info(f"Table now contains {count} rows")
            
            # Verify indexes
            cursor.execute("SHOW INDEX FROM Sequence_Similarity_Edges_4")
            indexes = cursor.fetchall()
            logger.info(f"Table has {len(indexes)} indexes:")
            for idx in indexes:
                logger.info(f"  - {idx[2]} on column {idx[4]}")

    except Exception as e:
        logger.error(f"Error creating similarity table: {e}")
        import traceback
        logger.error(traceback.format_exc())
        if connection:
            connection.rollback()
    finally:
        if connection:
            connection.close()
            logger.info("Database connection closed")

if __name__ == "__main__":
    feather_file_path = "C:/Users/mz30/protein_explorer/phosphosite_similarity_results.feather"
    create_similarity_table_from_feather(feather_file_path)