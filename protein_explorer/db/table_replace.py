import os
import logging
import pymysql
import time
import sys
import pandas as pd
import tempfile

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Database connection parameters - using the same as in your db.py
DB_HOST = os.environ.get('DB_HOST', '35.245.113.195')
DB_USER = os.environ.get('DB_USER', 'root')
DB_PASS = os.environ.get('DB_PASS', '@Bismark6')
DB_NAME = os.environ.get('DB_NAME', 'kinoplex-db')
DB_PORT = int(os.environ.get('DB_PORT', '3306'))

# Path to the new data file
NEW_DATA_FILE = r'C:/Users/mz30/protein_explorer/data/Combined_Kinome_10A_Master_Filtered_2.feather'

def update_structural_similarity_table():
    """Replace the Structural_Similarity_Edges table with new data using faster CSV import."""
    connection = None
    temp_csv = None
    
    try:
        # Step 1: Read Feather and convert to CSV for faster import
        logger.info(f"Reading feather file and converting to CSV...")
        start_time = time.time()
        df = pd.read_feather(NEW_DATA_FILE)
        logger.info(f"Loaded feather file with {len(df)} rows and {len(df.columns)} columns")
        logger.info(f"Column names: {list(df.columns)}")
        
        # Create a temp CSV file
        temp_csv = os.path.join(tempfile.gettempdir(), 'temp_structure_data.csv')
        df.to_csv(temp_csv, index=False)
        logger.info(f"Converted to CSV in {time.time() - start_time:.2f} seconds: {temp_csv}")
        
        # Connect to the database
        logger.info(f"Connecting to database {DB_NAME} at {DB_HOST}...")
        connection = pymysql.connect(
            host=DB_HOST,
            user=DB_USER,
            password=DB_PASS,
            database=DB_NAME,
            port=DB_PORT,
            charset='utf8mb4',
            local_infile=True  # Important for LOAD DATA LOCAL INFILE
        )
        
        with connection.cursor() as cursor:
            # Create backup of original table
            logger.info("Creating backup of current Structural_Similarity_Edges table...")
            cursor.execute("CREATE TABLE IF NOT EXISTS Structural_Similarity_Edges_Backup AS SELECT * FROM Structural_Similarity_Edges")
            
            # Check the backup
            cursor.execute("SELECT COUNT(*) FROM Structural_Similarity_Edges_Backup")
            backup_count = cursor.fetchone()[0]
            logger.info(f"Backup created with {backup_count} rows")
            
            # Drop the old table and create new one with same name but new columns
            logger.info("Dropping temporary table if it exists...")
            cursor.execute("DROP TABLE IF EXISTS Structural_Similarity_Edges_New")
            
            # Create new table schema based on DataFrame columns
            columns = []
            for col in df.columns:
                # Get the data type from the DataFrame
                dtype = df[col].dtype
                
                # Map pandas dtypes to MySQL data types
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
            
            # Create table with new schema
            create_table_query = f"""
            CREATE TABLE Structural_Similarity_Edges_New (
                id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
                {', '.join(columns)}
            ) ENGINE=InnoDB
            """
            
            logger.info(f"Creating new table with schema based on feather file...")
            cursor.execute(create_table_query)
            
            # Use LOAD DATA INFILE for fast import
            logger.info("Importing data using LOAD DATA INFILE...")
            start_time = time.time()
            
            # Prepare the LOAD DATA INFILE statement
            load_query = f"""
            LOAD DATA LOCAL INFILE '{temp_csv.replace('\\', '/')}' 
            INTO TABLE Structural_Similarity_Edges_New
            FIELDS TERMINATED BY ',' 
            OPTIONALLY ENCLOSED BY '"' 
            LINES TERMINATED BY '\\n'
            IGNORE 1 LINES
            """
            
            try:
                cursor.execute(load_query)
                connection.commit()
                logger.info(f"Data import completed in {time.time() - start_time:.2f} seconds")
            except Exception as e:
                logger.error(f"Error during LOAD DATA INFILE: {e}")
                logger.info("Attempting to use alternative import method...")
                
                # If LOAD DATA INFILE fails, try using sqlalchemy as a fallback
                try:
                    from sqlalchemy import create_engine
                    
                    # Create engine
                    connection_string = f"mysql+pymysql://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
                    engine = create_engine(connection_string)
                    
                    # Import data in chunks
                    logger.info("Importing data using SQLAlchemy in chunks...")
                    start_time = time.time()
                    
                    # Truncate the table first
                    cursor.execute("TRUNCATE TABLE Structural_Similarity_Edges_New")
                    connection.commit()
                    
                    # Import in chunks of 10000 rows
                    df.to_sql('Structural_Similarity_Edges_New', engine, 
                             if_exists='append', index=False, 
                             chunksize=10000, method='multi')
                    
                    logger.info(f"SQLAlchemy import completed in {time.time() - start_time:.2f} seconds")
                except Exception as e2:
                    logger.error(f"Alternative import also failed: {e2}")
                    raise
            
            # Verify the imported data
            cursor.execute("SELECT COUNT(*) FROM Structural_Similarity_Edges_New")
            new_count = cursor.fetchone()[0]
            logger.info(f"New table has {new_count} rows")
            
            if new_count == 0:
                raise ValueError("No data was imported into the new table. Aborting replacement.")
            
            # Replace the original table with the new one
            logger.info("Replacing original table with new table...")
            cursor.execute("RENAME TABLE Structural_Similarity_Edges TO Structural_Similarity_Edges_Old, Structural_Similarity_Edges_New TO Structural_Similarity_Edges")
            connection.commit()
            
            # Create indexes (if needed)
            logger.info("Creating indices on the new table...")
            index_statements = []
            
            # Check if expected columns exist for the standard indexes
            if 'Query' in df.columns and 'RMSD' in df.columns:
                index_statements.append("CREATE INDEX idx_query_rmsd ON Structural_Similarity_Edges(`Query`(50), `RMSD`)")
            
            if 'Target' in df.columns:
                index_statements.append("CREATE INDEX idx_target ON Structural_Similarity_Edges(`Target`(50))")
            
            # Execute index creation
            start_time = time.time()
            for statement in index_statements:
                try:
                    logger.info(f"Executing: {statement}")
                    cursor.execute(statement)
                    connection.commit()
                except pymysql.err.OperationalError as e:
                    if "Duplicate key name" in str(e):
                        logger.info(f"Index already exists: {e}")
                    else:
                        logger.error(f"Error creating index: {e}")
            
            end_time = time.time()
            logger.info(f"Index creation completed in {end_time - start_time:.2f} seconds")
            
            # Check if user wants to drop the old table
            logger.info("Do you want to drop the old table Structural_Similarity_Edges_Old? (yes/no)")
            choice = input().strip().lower()
            if choice == 'yes':
                cursor.execute("DROP TABLE Structural_Similarity_Edges_Old")
                connection.commit()
                logger.info("Old table dropped")
            else:
                logger.info("Old table retained for safety")
            
            # Analyze table for query optimization
            logger.info("Analyzing table for optimal query performance...")
            cursor.execute("ANALYZE TABLE Structural_Similarity_Edges")
            analyze_result = cursor.fetchall()
            for result in analyze_result:
                logger.info(f"Table analysis: {result}")
            
            logger.info("Table update completed successfully")
            
    except Exception as e:
        logger.error(f"Error updating structural similarity table: {e}")
        if connection:
            connection.rollback()
        logger.info("Rolling back changes. The original table remains unchanged.")
        raise
    finally:
        # Clean up
        if connection:
            connection.close()
            logger.info("Database connection closed")
        
        # Remove temporary CSV file
        if temp_csv and os.path.exists(temp_csv):
            try:
                os.remove(temp_csv)
                logger.info(f"Temporary CSV file removed")
            except Exception as e:
                logger.warning(f"Failed to remove temporary CSV file: {e}")

if __name__ == "__main__":
    try:
        update_structural_similarity_table()
    except Exception as e:
        logger.error(f"Script execution failed: {e}")
        sys.exit(1)