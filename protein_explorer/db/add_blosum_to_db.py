"""
Script to add a central motif column to the Phosphosite_Supplementary_Data table
and create necessary indexes for optimized BLOSUM matching.
"""

import os
import sys
import time
import logging
from sqlalchemy import text


# Debug the import path
print("Current working directory:", os.getcwd())
print("Python path:", sys.path)

# Add the current directory and parent directory to the path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
grandparent_dir = os.path.dirname(parent_dir)

print("Script location:", __file__)
print("Current directory:", current_dir)
print("Parent directory:", parent_dir)
print("Grandparent directory:", grandparent_dir)

# Add parent directory to path (should be the protein_explorer folder)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

# Direct import with explicit path
try:
    from db import execute_query, init_db, engine
    print("Successfully imported from db module in the same directory")
except ImportError:
    try:
        # Try relative import
        from .db import execute_query, init_db, engine
        print("Successfully imported from .db module")
    except ImportError:
        try:
            # Try with full path
            sys.path.insert(0, current_dir)
            from db import execute_query, init_db, engine
            print("Successfully imported from db module with current_dir in path")
        except ImportError:
            print("Still can't import - trying one more approach")
            try:
                # Find the actual db.py file and import directly
                import importlib.util
                import os
                
                # Search for db.py
                for root, dirs, files in os.walk(grandparent_dir):
                    if 'db.py' in files:
                        print(f"Found db.py at: {os.path.join(root, 'db.py')}")
                        spec = importlib.util.spec_from_file_location("db", os.path.join(root, "db.py"))
                        db = importlib.util.module_from_spec(spec)
                        spec.loader.exec_module(db)
                        execute_query = db.execute_query
                        init_db = db.init_db
                        print("Successfully imported using direct file location")
                        break
                else:
                    print("Could not find db.py anywhere in the directory structure")
            except Exception as e:
                print(f"Final import attempt failed: {e}")
                sys.exit(1)

# Configure logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def execute_non_query(query, params=None):
    """Execute a non-query SQL command (INSERT, UPDATE, CREATE, etc.) that doesn't return rows."""
    try:
        with engine.connect() as conn:
            if params:
                conn.execute(text(query), params)
            else:
                conn.execute(text(query))
            conn.commit()
        return True
    except Exception as e:
        logger.error(f"Database error executing command: {e}")
        logger.error(f"Query was: {query}")
        logger.error(f"Parameters were: {params}")
        return False

def add_central_motif_column():
    """Add a new column with precomputed central motifs (-4:+4) to the database."""
    try:
        logger.info("Checking if central_motif column already exists")
        check_query = """
            SELECT COUNT(*) as count FROM INFORMATION_SCHEMA.COLUMNS 
            WHERE TABLE_NAME = 'Phosphosite_Supplementary_Data' 
            AND COLUMN_NAME = 'central_motif'
        """
        df = execute_query(check_query)
        
        if df.empty:
            logger.error("Could not check if column exists - table may not exist")
            return False
            
        if df.iloc[0]['count'] > 0:
            logger.info("central_motif column already exists. Skipping column creation.")
        else:
            logger.info("Adding central_motif column to Phosphosite_Supplementary_Data")
            alter_query = """
                ALTER TABLE Phosphosite_Supplementary_Data 
                ADD COLUMN central_motif VARCHAR(9) AFTER `SITE_+/-7_AA`;
            """
            if execute_non_query(alter_query):
                logger.info("central_motif column added successfully")
            else:
                logger.error("Failed to add central_motif column")
                return False
        
        # Update all rows with the central motif
        logger.info("Updating central_motif values for all rows...")
        start_time = time.time()
        
        # First update pass - use the standard central window extraction
        update_query = """
            UPDATE Phosphosite_Supplementary_Data
            SET central_motif = 
                CASE 
                    WHEN `SITE_+/-7_AA` IS NOT NULL AND LENGTH(`SITE_+/-7_AA`) >= 9 
                    THEN SUBSTRING(`SITE_+/-7_AA`, LENGTH(`SITE_+/-7_AA`)/2 - 4, 9)
                    ELSE NULL
                END
            WHERE `SITE_+/-7_AA` IS NOT NULL AND (central_motif IS NULL OR central_motif = '');
        """
        if execute_non_query(update_query):
            logger.info("Central motif values updated successfully")
        else:
            logger.error("Failed to update central motif values")
            return False
        
        end_time = time.time()
        logger.info(f"Column values updated in {end_time - start_time:.2f} seconds")
        
        # Check if index exists
        logger.info("Checking if index on central_motif already exists")
        index_check_query = """
            SELECT COUNT(*) as count FROM INFORMATION_SCHEMA.STATISTICS
            WHERE TABLE_NAME = 'Phosphosite_Supplementary_Data' 
            AND INDEX_NAME = 'idx_central_motif'
        """
        df = execute_query(index_check_query)
        
        if df.iloc[0]['count'] > 0:
            logger.info("Index idx_central_motif already exists")
        else:
            # Create an index on this column for faster lookups
            logger.info("Creating index on central_motif column...")
            index_query = """
                CREATE INDEX idx_central_motif ON Phosphosite_Supplementary_Data(central_motif);
            """
            if execute_non_query(index_query):
                logger.info("Index created successfully")
            else:
                logger.error("Failed to create index")
                return False
        
        # Count how many rows have valid central_motif values
        count_query = """
            SELECT COUNT(*) as count FROM Phosphosite_Supplementary_Data
            WHERE central_motif IS NOT NULL;
        """
        df = execute_query(count_query)
        logger.info(f"Total rows with valid central_motif: {df.iloc[0]['count']}")
        
        # Create additional indices for faster filtering
        logger.info("Creating additional indices for optimized BLOSUM matching...")
        
        # Index for motif start (first 3 characters)
        start_index_check_query = """
            SELECT COUNT(*) as count FROM INFORMATION_SCHEMA.STATISTICS
            WHERE TABLE_NAME = 'Phosphosite_Supplementary_Data' 
            AND INDEX_NAME = 'idx_motif_start'
        """
        df = execute_query(start_index_check_query)
        
        if df.iloc[0]['count'] == 0:
            start_index_query = """
                CREATE INDEX idx_motif_start ON Phosphosite_Supplementary_Data(LEFT(central_motif, 3));
            """
            if execute_non_query(start_index_query):
                logger.info("Index on motif start created successfully")
            else:
                logger.error("Failed to create index on motif start")
        else:
            logger.info("Index on motif start already exists")
        
        # Index for motif middle (middle 3 characters)
        middle_index_check_query = """
            SELECT COUNT(*) as count FROM INFORMATION_SCHEMA.STATISTICS
            WHERE TABLE_NAME = 'Phosphosite_Supplementary_Data' 
            AND INDEX_NAME = 'idx_motif_middle'
        """
        df = execute_query(middle_index_check_query)
        
        if df.iloc[0]['count'] == 0:
            middle_index_query = """
                CREATE INDEX idx_motif_middle ON Phosphosite_Supplementary_Data(SUBSTRING(central_motif, 4, 3));
            """
            if execute_non_query(middle_index_query):
                logger.info("Index on motif middle created successfully")
            else:
                logger.error("Failed to create index on motif middle")
        else:
            logger.info("Index on motif middle already exists")
        
        # Index for motif end (last 3 characters)
        end_index_check_query = """
            SELECT COUNT(*) as count FROM INFORMATION_SCHEMA.STATISTICS
            WHERE TABLE_NAME = 'Phosphosite_Supplementary_Data' 
            AND INDEX_NAME = 'idx_motif_end'
        """
        df = execute_query(end_index_check_query)
        
        if df.iloc[0]['count'] == 0:
            end_index_query = """
                CREATE INDEX idx_motif_end ON Phosphosite_Supplementary_Data(RIGHT(central_motif, 3));
            """
            if execute_non_query(end_index_query):
                logger.info("Index on motif end created successfully")
            else:
                logger.error("Failed to create index on motif end")
        else:
            logger.info("Index on motif end already exists")
            
        logger.info("Central motif column added, populated, and indexed successfully")
        return True
    except Exception as e:
        logger.error(f"Error adding central motif column: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def print_central_motif_samples(limit=20):
    """Print sample values from the central_motif column for verification."""
    try:
        sample_query = """
            SELECT PhosphositeID, `SITE_+/-7_AA` as full_motif, central_motif
            FROM Phosphosite_Supplementary_Data
            WHERE central_motif IS NOT NULL
            LIMIT :limit
        """
        df = execute_query(sample_query, {"limit": limit})
        
        if df.empty:
            print("No central motif data found")
            return
        
        print("\n===== CENTRAL MOTIF SAMPLES =====")
        print(f"{'Site ID':<30} {'Full Motif':<20} {'Central Motif':<10}")
        print("-" * 60)
        
        for _, row in df.iterrows():
            site_id = row['PhosphositeID']
            full_motif = row['full_motif']
            central_motif = row['central_motif']
            
            print(f"{site_id:<30} {full_motif:<20} {central_motif:<10}")
        
        print("=" * 60)
        
        # Also show stats on central motif lengths
        length_query = """
            SELECT LENGTH(central_motif) as length, COUNT(*) as count
            FROM Phosphosite_Supplementary_Data
            WHERE central_motif IS NOT NULL
            GROUP BY LENGTH(central_motif)
            ORDER BY LENGTH(central_motif)
        """
        df_lengths = execute_query(length_query)
        
        if not df_lengths.empty:
            print("\n===== CENTRAL MOTIF LENGTH DISTRIBUTION =====")
            print(f"{'Length':<10} {'Count':<10}")
            print("-" * 20)
            
            for _, row in df_lengths.iterrows():
                length = row['length']
                count = row['count']
                print(f"{length:<10} {count:<10}")
            
            print("=" * 20)
        
    except Exception as e:
        print(f"Error printing central motif samples: {e}")

if __name__ == "__main__":
    logger.info("Initializing database connection...")
    init_db()
    
    logger.info("Starting database modifications...")
    
    # Add central_motif column to the main table
    if add_central_motif_column():
        logger.info("Database modifications completed: central_motif column added and indexed.")
        
        # Print sample motifs to verify
        print_central_motif_samples(20)
    else:
        logger.error("Database modifications failed: could not add central_motif column.")
        sys.exit(1)
    
    logger.info("All database modifications completed successfully.")