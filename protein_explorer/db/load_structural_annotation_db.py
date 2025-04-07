#!/usr/bin/env python
"""
Script to bulk load structural annotations from CSV into MySQL database.
Uses LOAD DATA INFILE for optimal performance with multiple fallback options.

Usage:
  python add_structural_annotations.py
"""

import os
import sys
import time
import logging
import subprocess
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('structural_load.log')
    ]
)

logger = logging.getLogger(__name__)

# Database connection parameters
DB_USER = os.environ.get('DB_USER', 'root')
DB_PASS = os.environ.get('DB_PASS', '@Bismark6')
DB_HOST = os.environ.get('DB_HOST', '35.245.113.195')
DB_PORT = os.environ.get('DB_PORT', '3306')
DB_NAME = os.environ.get('DB_NAME', 'kinoplex-db')

# Path to the CSV file
CSV_FILE_PATH = "C:/Users/mz30/protein_explorer/protein_explorer/precomputing_data/STY_Structural_Annotations.csv"
# Path for temporary CSV file
TEMP_CSV_PATH = "C:/Users/mz30/protein_explorer/protein_explorer/precomputing_data/temp_structural_data.csv"

def prepare_csv_for_import():
    """Prepare CSV data for import by handling any problematic formatting."""
    try:
        logger.info(f"Reading CSV file: {CSV_FILE_PATH}")
        
        # Read the CSV file
        df = pd.read_csv(CSV_FILE_PATH)
        
        # Log column information
        logger.info(f"CSV columns: {df.columns.tolist()}")
        logger.info(f"CSV rows: {len(df)}")
        
        # Create combined PhosphositeID column for better indexing if not already exists
        if 'PhosphositeID' not in df.columns:
            if 'UniProtID' in df.columns and 'ResidueNumber' in df.columns:
                df['PhosphositeID'] = df['UniProtID'] + '_' + df['ResidueNumber'].astype(str)
                logger.info("Created PhosphositeID column")
        
        # Handle any NULL values or special characters that might cause issues
        for column in df.columns:
            if df[column].dtype == 'object':  # String columns
                df[column] = df[column].fillna('')  # Replace NaN with empty string
        
        # Save to a temporary CSV file for import
        logger.info(f"Writing prepared data to: {TEMP_CSV_PATH}")
        df.to_csv(TEMP_CSV_PATH, index=False, quoting=1)  # quoting=1 for quotes around text fields
        
        logger.info(f"Successfully prepared {len(df)} rows for import")
        return True
    except Exception as e:
        logger.error(f"Error preparing CSV data: {e}")
        return False

def try_connect_mysql():
    """Try to connect to MySQL using available drivers."""
    # Try mysql.connector
    try:
        import mysql.connector
        conn = mysql.connector.connect(
            host=DB_HOST,
            user=DB_USER,
            password=DB_PASS,
            database=DB_NAME,
            port=int(DB_PORT),
            allow_local_infile=True  # Important for LOAD DATA LOCAL INFILE
        )
        logger.info("Connected to MySQL using mysql.connector")
        return conn, "connector"
    except ImportError:
        logger.warning("mysql.connector not available")
    except Exception as e:
        logger.warning(f"Could not connect with mysql.connector: {e}")
    
    # Try MySQLdb
    try:
        import MySQLdb
        conn = MySQLdb.connect(
            host=DB_HOST,
            user=DB_USER,
            passwd=DB_PASS,
            db=DB_NAME,
            port=int(DB_PORT),
            local_infile=1  # Enable local infile
        )
        logger.info("Connected to MySQL using MySQLdb")
        return conn, "mysqldb"
    except ImportError:
        logger.warning("MySQLdb not available")
    except Exception as e:
        logger.warning(f"Could not connect with MySQLdb: {e}")
    
    # Try pymysql
    try:
        import pymysql
        conn = pymysql.connect(
            host=DB_HOST,
            user=DB_USER,
            password=DB_PASS,
            database=DB_NAME,
            port=int(DB_PORT),
            local_infile=True  # Enable local infile
        )
        logger.info("Connected to MySQL using pymysql")
        return conn, "pymysql"
    except ImportError:
        logger.warning("pymysql not available")
    except Exception as e:
        logger.warning(f"Could not connect with pymysql: {e}")
    
    logger.error("Could not connect to MySQL with any available driver")
    return None, None

def create_table(cursor):
    """Create the table with proper structure and indexes."""
    logger.info("Creating table...")
    
    create_table_query = """
    CREATE TABLE IF NOT EXISTS `STY_Structural_Annotations` (
        `UniProtID` VARCHAR(50) NOT NULL,
        `ResidueNumber` INT NOT NULL,
        `Site` VARCHAR(50) NOT NULL,
        `ResidueType` CHAR(1) NOT NULL,
        `Motif` VARCHAR(100),
        `pLDDT` FLOAT,
        `NeighborCount` INT,
        `ChainID` VARCHAR(10),
        `SecondaryStructure` VARCHAR(50),
        `SecondaryStructureDesc` VARCHAR(100),
        `HydroxylExposure` FLOAT,
        `BackboneContacts` INT,
        `HydrogenBonds` INT,
        `MotifPLDDT` FLOAT,
        `SeqDistToHelix` FLOAT,
        `SpatialDistToHelix` FLOAT,
        `SeqDistToStrand` FLOAT,
        `SpatialDistToStrand` FLOAT,
        `SeqDistToTurn` FLOAT,
        `SpatialDistToTurn` FLOAT,
        `hydrophobic_count` INT,
        `polar_count` INT,
        `positive_count` INT,
        `negative_count` INT,
        `net_charge` FLOAT,
        `small_count` INT,
        `medium_count` INT,
        `large_count` INT,
        `hydrophobic_ratio` FLOAT,
        `charged_ratio` FLOAT,
        `HSE_CA_U` INT,
        `HSE_CA_D` INT,
        `HSE_CA_RATIO` FLOAT,
        `HSE_CB_U` INT,
        `HSE_CB_D` INT,
        `HSE_CB_RATIO` FLOAT,
        `PhosphositeID` VARCHAR(50) NOT NULL,
        PRIMARY KEY (`Site`),
        INDEX `idx_uniprot_residue` (`UniProtID`, `ResidueNumber`),
        INDEX `idx_phosphosite_id` (`PhosphositeID`),
        INDEX `idx_residue_type` (`ResidueType`),
        INDEX `idx_plddt` (`pLDDT`),
        INDEX `idx_secondary_structure` (`SecondaryStructure`)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci
    """
    
    cursor.execute(create_table_query)
    cursor.connection.commit()
    logger.info("Table created successfully")

def bulk_load_data(connection, driver_type):
    """Bulk load data from CSV file into MySQL table."""
    try:
        cursor = connection.cursor()
        
        # Check if table exists
        try:
            cursor.execute("SELECT COUNT(*) FROM `STY_Structural_Annotations`")
            count = cursor.fetchone()[0]
            logger.info(f"Table exists with {count} rows")
            
            # Ask user what to do
            choice = input("Table already exists. Choose an option:\n"
                           "1. Drop and recreate table\n"
                           "2. Truncate table (faster but keeps structure)\n"
                           "3. Append to existing table\n"
                           "4. Exit without changes\n"
                           "Enter choice (1-4): ")
            
            if choice == '1':
                logger.info("Dropping table...")
                cursor.execute("DROP TABLE IF EXISTS `STY_Structural_Annotations`")
                connection.commit()
                # Recreate table
                create_table(cursor)
            elif choice == '2':
                logger.info("Truncating table...")
                cursor.execute("TRUNCATE TABLE `STY_Structural_Annotations`")
                connection.commit()
            elif choice == '3':
                logger.info("Will append data to existing table")
            elif choice == '4':
                logger.info("Exiting without changes")
                return False
            else:
                logger.warning("Invalid choice. Exiting.")
                return False
        except Exception as e:
            logger.info(f"Table doesn't exist yet: {e}")
            # Create the table
            create_table(cursor)
        
        # Start the bulk load process
        logger.info("Starting bulk load process...")
        
        # Disable keys/indexes for faster loading
        logger.info("Disabling keys for faster loading...")
        cursor.execute("ALTER TABLE `STY_Structural_Annotations` DISABLE KEYS")
        connection.commit()
        
        # Try different approaches for setting local_infile
        try:
            # Try to set global variable - may fail if user doesn't have privileges
            logger.info("Trying to set GLOBAL local_infile...")
            cursor.execute("SET GLOBAL local_infile = 1")
        except Exception as e:
            logger.warning(f"Could not set GLOBAL local_infile: {e}")
            try:
                # Try to set session variable - less privileges required
                logger.info("Trying to set SESSION local_infile...")
                cursor.execute("SET SESSION local_infile = 1")
            except Exception as e2:
                logger.warning(f"Could not set SESSION local_infile: {e2}")
                logger.info("Continuing anyway - connection may have local_infile enabled already")
        
        # Prepare the LOAD DATA INFILE statement with proper path formatting
        csv_path = TEMP_CSV_PATH.replace('\\', '/')
        
        load_query = """
        LOAD DATA LOCAL INFILE '{}'
        INTO TABLE `STY_Structural_Annotations`
        FIELDS TERMINATED BY ','
        ENCLOSED BY '"'
        LINES TERMINATED BY '\\n'
        IGNORE 1 LINES;
        """.format(csv_path)
        
        # Execute the bulk load
        start_time = time.time()
        logger.info("Executing LOAD DATA INFILE...")
        try:
            cursor.execute(load_query)
            connection.commit()
        except Exception as e:
            logger.error(f"LOAD DATA INFILE failed: {e}")
            
            # If we failed, try alternative approach with batch inserts
            if input("Do you want to try batch insertion instead? (y/n): ").lower() == 'y':
                logger.info("Attempting batch insertion...")
                
                # Read CSV in chunks
                import pandas as pd
                chunk_size = 5000
                chunks = pd.read_csv(TEMP_CSV_PATH, chunksize=chunk_size)
                
                # Track progress
                total_rows = 0
                
                for i, chunk in enumerate(chunks):
                    # Create placeholders for SQL query
                    columns = chunk.columns.tolist()
                    placeholders = ', '.join(['%s'] * len(columns))
                    column_names = ', '.join(f'`{col}`' for col in columns)
                    
                    # Convert DataFrame chunk to list of tuples
                    values = [tuple(x) for x in chunk.values]
                    
                    # Build insert query
                    insert_query = f"""
                    INSERT INTO `STY_Structural_Annotations` ({column_names})
                    VALUES ({placeholders})
                    """
                    
                    # Execute batch insert
                    logger.info(f"Inserting chunk {i+1} with {len(values)} rows...")
                    cursor.executemany(insert_query, values)
                    connection.commit()
                    
                    total_rows += len(values)
                    logger.info(f"Inserted {total_rows} rows so far...")
                
                logger.info(f"Completed batch insertion of {total_rows} rows")
            else:
                logger.error("Bulk load failed and batch insertion was declined.")
                return False
        
        # Enable keys/indexes
        logger.info("Re-enabling keys and building indexes...")
        cursor.execute("ALTER TABLE `STY_Structural_Annotations` ENABLE KEYS")
        connection.commit()
        
        # Get final row count
        cursor.execute("SELECT COUNT(*) FROM `STY_Structural_Annotations`")
        final_count = cursor.fetchone()[0]
        
        # Log completion
        elapsed_time = time.time() - start_time
        logger.info(f"Bulk load completed in {elapsed_time:.2f} seconds")
        logger.info(f"Loaded {final_count} rows into STY_Structural_Annotations")
        
        return True
    except Exception as e:
        logger.error(f"Error during bulk load: {e}")
        return False

def try_direct_mysql_command():
    """Try using the mysql command-line client directly."""
    try:
        logger.info("Attempting to use mysql command-line client...")
        
        # Create SQL file with commands
        sql_file_path = "C:/Users/mz30/protein_explorer/protein_explorer/precomputing_data/load_structural_data.sql"
        csv_path = TEMP_CSV_PATH.replace('\\', '/')
        
        with open(sql_file_path, 'w') as f:
            f.write("""
            -- Create or replace the table
            DROP TABLE IF EXISTS `STY_Structural_Annotations`;
            
            CREATE TABLE `STY_Structural_Annotations` (
                `UniProtID` VARCHAR(50) NOT NULL,
                `ResidueNumber` INT NOT NULL,
                `Site` VARCHAR(50) NOT NULL,
                `ResidueType` CHAR(1) NOT NULL,
                `Motif` VARCHAR(100),
                `pLDDT` FLOAT,
                `NeighborCount` INT,
                `ChainID` VARCHAR(10),
                `SecondaryStructure` VARCHAR(50),
                `SecondaryStructureDesc` VARCHAR(100),
                `HydroxylExposure` FLOAT,
                `BackboneContacts` INT,
                `HydrogenBonds` INT,
                `MotifPLDDT` FLOAT,
                `SeqDistToHelix` FLOAT,
                `SpatialDistToHelix` FLOAT,
                `SeqDistToStrand` FLOAT,
                `SpatialDistToStrand` FLOAT,
                `SeqDistToTurn` FLOAT,
                `SpatialDistToTurn` FLOAT,
                `hydrophobic_count` INT,
                `polar_count` INT,
                `positive_count` INT,
                `negative_count` INT,
                `net_charge` FLOAT,
                `small_count` INT,
                `medium_count` INT,
                `large_count` INT,
                `hydrophobic_ratio` FLOAT,
                `charged_ratio` FLOAT,
                `HSE_CA_U` INT,
                `HSE_CA_D` INT,
                `HSE_CA_RATIO` FLOAT,
                `HSE_CB_U` INT,
                `HSE_CB_D` INT,
                `HSE_CB_RATIO` FLOAT,
                `PhosphositeID` VARCHAR(50) NOT NULL,
                PRIMARY KEY (`Site`),
                INDEX `idx_uniprot_residue` (`UniProtID`, `ResidueNumber`),
                INDEX `idx_phosphosite_id` (`PhosphositeID`),
                INDEX `idx_residue_type` (`ResidueType`),
                INDEX `idx_plddt` (`pLDDT`),
                INDEX `idx_secondary_structure` (`SecondaryStructure`)
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
            
            -- Disable keys for faster loading
            ALTER TABLE `STY_Structural_Annotations` DISABLE KEYS;
            
            -- Load data
            LOAD DATA LOCAL INFILE '{}' 
            INTO TABLE `STY_Structural_Annotations`
            FIELDS TERMINATED BY ','
            ENCLOSED BY '"'
            LINES TERMINATED BY '\\n'
            IGNORE 1 LINES;
            
            -- Re-enable keys
            ALTER TABLE `STY_Structural_Annotations` ENABLE KEYS;
            
            -- Show row count
            SELECT COUNT(*) FROM `STY_Structural_Annotations`;
            """.format(csv_path))
        
        # Try to find the mysql executable
        mysql_exe = "mysql"  # Default assumes it's in PATH
        
        # Try common install locations
        possible_paths = [
            r"C:\Program Files\MySQL\MySQL Server 8.0\bin\mysql.exe",
            r"C:\Program Files (x86)\MySQL\MySQL Server 8.0\bin\mysql.exe",
            r"C:\xampp\mysql\bin\mysql.exe"
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                mysql_exe = path
                break
        
        # Build command
        cmd = [
            mysql_exe,
            f"--host={DB_HOST}",
            f"--user={DB_USER}",
            f"--password={DB_PASS}",
            f"--port={DB_PORT}",
            f"--database={DB_NAME}",
            "--local-infile=1",
            f"--execute=source {sql_file_path}"
        ]
        
        # Execute command
        logger.info(f"Executing: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            logger.info("MySQL command executed successfully")
            logger.info(f"Output: {result.stdout}")
            return True
        else:
            logger.error(f"MySQL command failed: {result.stderr}")
            return False
    except Exception as e:
        logger.error(f"Error executing MySQL command: {e}")
        return False

def main():
    """Main function to execute the bulk load process."""
    logger.info("Starting bulk loading process")
    
    # Prepare CSV for import
    if not os.path.exists(TEMP_CSV_PATH):
        if not prepare_csv_for_import():
            logger.error("Failed to prepare CSV. Exiting.")
            return False
    else:
        logger.info(f"Using existing prepared CSV file: {TEMP_CSV_PATH}")
    
    # Try to use Python's MySQL libraries first
    connection, driver_type = try_connect_mysql()
    
    if connection:
        # Bulk load using Python
        success = bulk_load_data(connection, driver_type)
        connection.close()
        
        if success:
            logger.info("Bulk load completed successfully using Python")
            
            # Clean up temporary file if desired
            if input("Delete temporary CSV file? (y/n): ").lower() == 'y':
                try:
                    os.remove(TEMP_CSV_PATH)
                    logger.info(f"Deleted temporary file: {TEMP_CSV_PATH}")
                except Exception as e:
                    logger.error(f"Could not delete temporary file: {e}")
            
            return True
        else:
            logger.warning("Bulk load failed using Python, trying direct mysql command...")
    else:
        logger.warning("Could not connect to MySQL with Python drivers, trying direct mysql command...")
    
    # If Python approach failed, try direct mysql command
    success = try_direct_mysql_command()
    
    if success:
        logger.info("Bulk load completed successfully using mysql command")
        return True
    else:
        logger.error("All bulk load methods failed")
        return False

if __name__ == "__main__":
    try:
        success = main()
        if success:
            logger.info("Script completed successfully")
            sys.exit(0)
        else:
            logger.error("Script failed")
            sys.exit(1)
    except KeyboardInterrupt:
        logger.info("Script interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unhandled exception: {e}", exc_info=True)
        sys.exit(1)