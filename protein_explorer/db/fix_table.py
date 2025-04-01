import os
import logging
import pymysql
import time
import sys
import pandas as pd

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Database connection parameters
DB_HOST = os.environ.get('DB_HOST', '35.245.113.195')
DB_USER = os.environ.get('DB_USER', 'root')
DB_PASS = os.environ.get('DB_PASS', '@Bismark6')
DB_NAME = os.environ.get('DB_NAME', 'kinoplex-db')
DB_PORT = int(os.environ.get('DB_PORT', '3306'))

# Path to the original data file
DATA_FILE = r'C:/Users/mz30/protein_explorer/data/Combined_Kinome_10A_Master_Filtered_2.feather'

def fix_structural_similarity_edges():
    """Fix column alignment in the Structural_Similarity_Edges table."""
    connection = None
    
    try:
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
            # Create backup if it doesn't exist
            logger.info("Creating backup if it doesn't exist...")
            cursor.execute("CREATE TABLE IF NOT EXISTS Structural_Similarity_Edges_Backup AS SELECT * FROM Structural_Similarity_Edges")
            
            # Check if this is a fresh restoration or if we need to use the backup
            cursor.execute("SELECT COUNT(*) FROM Structural_Similarity_Edges_Backup")
            backup_count = cursor.fetchone()[0]
            logger.info(f"Backup table has {backup_count} rows")
            
            # Load original data from feather file
            logger.info(f"Reading original data from {DATA_FILE}...")
            df = pd.read_feather(DATA_FILE)
            logger.info(f"Original data has {len(df)} rows and {len(df.columns)} columns")
            logger.info(f"Original columns: {list(df.columns)}")
            
            # Drop existing table
            logger.info("Dropping existing table...")
            cursor.execute("DROP TABLE IF EXISTS Structural_Similarity_Edges")
            
            # Create new table with correct schema
            columns = []
            for col in df.columns:
                # Map pandas dtypes to MySQL data types
                dtype = df[col].dtype
                if pd.api.types.is_float_dtype(dtype):
                    sql_type = "DOUBLE"
                elif pd.api.types.is_integer_dtype(dtype):
                    sql_type = "INT"
                elif pd.api.types.is_bool_dtype(dtype):
                    sql_type = "BOOL"
                else:
                    # For string/object types, use VARCHAR(255)
                    sql_type = "VARCHAR(255)"
                
                # Escape column name to handle special characters
                escaped_col = f"`{col}`"
                columns.append(f"{escaped_col} {sql_type}")
            
            # Create table with explicit columns from the feather file
            create_table_query = f"""
            CREATE TABLE Structural_Similarity_Edges (
                id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
                {', '.join(columns)}
            )
            """
            
            logger.info("Creating new table with correct schema...")
            cursor.execute(create_table_query)
            
            # Create a temporary CSV file
            import tempfile
            temp_csv = os.path.join(tempfile.gettempdir(), 'struct_edges.csv')
            logger.info(f"Creating temporary CSV file: {temp_csv}")
            df.to_csv(temp_csv, index=False)
            
            # Use LOAD DATA INFILE for fast import with explicit column names
            logger.info("Loading data with explicit column mapping...")
            load_query = f"""
            LOAD DATA LOCAL INFILE '{temp_csv.replace('\\', '/')}' 
            INTO TABLE Structural_Similarity_Edges
            FIELDS TERMINATED BY ',' 
            OPTIONALLY ENCLOSED BY '"' 
            LINES TERMINATED BY '\\n'
            IGNORE 1 LINES
            ({', '.join([f'`{col}`' for col in df.columns])})
            """
            
            try:
                cursor.execute(load_query)
                connection.commit()
                logger.info("Data loaded successfully")
            except Exception as e:
                logger.error(f"Error during LOAD DATA INFILE: {e}")
                
                # Alternative approach using pandas to_sql
                logger.info("Trying alternative import method via SQLAlchemy...")
                
                # Reset the table
                cursor.execute("TRUNCATE TABLE Structural_Similarity_Edges")
                connection.commit()
                
                # Use SQLAlchemy for bulk import
                from sqlalchemy import create_engine
                
                password = DB_PASS.replace('@', '%40')  # URL encode special characters
                engine = create_engine(f"mysql+pymysql://{DB_USER}:{password}@{DB_HOST}:{DB_PORT}/{DB_NAME}")
                
                # Import in chunks
                chunk_size = 1000000
                for i in range(0, len(df), chunk_size):
                    chunk = df.iloc[i:i + chunk_size]
                    chunk.to_sql('Structural_Similarity_Edges', engine, if_exists='append', index=False)
                    logger.info(f"Imported chunk {i//chunk_size + 1} ({i} to {min(i+chunk_size, len(df))}) rows")
                
                logger.info("Data imported via SQLAlchemy")
            
            # Remove temporary CSV
            if os.path.exists(temp_csv):
                os.remove(temp_csv)
                logger.info("Temporary CSV file removed")
            
            # Create indexes
            logger.info("Creating indexes...")
            indexes = [
                "CREATE INDEX idx_query_rmsd ON Structural_Similarity_Edges(`Query`(50), `RMSD`)",
                "CREATE INDEX idx_target ON Structural_Similarity_Edges(`Target`(50))"
            ]
            
            for idx_sql in indexes:
                try:
                    cursor.execute(idx_sql)
                    connection.commit()
                    logger.info(f"Created index: {idx_sql}")
                except Exception as e:
                    logger.error(f"Error creating index: {e}")
            
            # Verify the data
            cursor.execute("SELECT COUNT(*) FROM Structural_Similarity_Edges")
            count = cursor.fetchone()[0]
            logger.info(f"New table has {count} rows")
            
            # Show sample data
            cursor.execute("SELECT * FROM Structural_Similarity_Edges LIMIT 5")
            sample = cursor.fetchall()
            logger.info("Sample data from new table:")
            for row in sample:
                logger.info(f"  {row[:5]}...")  # Show first few fields
            
            logger.info("Table column alignment fixed successfully")
            
    except Exception as e:
        logger.error(f"Error fixing table: {e}")
        if connection:
            connection.rollback()
        raise
    finally:
        if connection:
            connection.close()
            logger.info("Database connection closed")

if __name__ == "__main__":
    fix_structural_similarity_edges()