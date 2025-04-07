#!/usr/bin/env python
"""
Script to bulk load sequence similarity data into the MySQL database.
Uses LOAD DATA INFILE for optimal performance with the 70 million row dataset.

Usage:
  python add_new_sequence_sim_table.py
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
        logging.FileHandler('bulk_load.log')
    ]
)

logger = logging.getLogger(__name__)

# Database connection parameters
DB_USER = os.environ.get('DB_USER', 'root')
DB_PASS = os.environ.get('DB_PASS', '@Bismark6')
DB_HOST = os.environ.get('DB_HOST', '35.245.113.195')
DB_PORT = os.environ.get('DB_PORT', '3306')
DB_NAME = os.environ.get('DB_NAME', 'kinoplex-db')

# Path to the feather file
FEATHER_FILE_PATH = "C:/Users/mz30/Downloads/Motif_Similarity_Blosum62_Filtered_0.25_cutoff.feather"
# Path for temporary CSV file
CSV_FILE_PATH = "C:/Users/mz30/Downloads/temp_similarity_data.csv"

def convert_feather_to_csv():
    """Convert feather file to CSV format for bulk loading."""
    try:
        logger.info(f"Converting feather file to CSV: {FEATHER_FILE_PATH}")
        
        # Check if we can use pyarrow
        try:
            import pyarrow.feather as feather
            df = feather.read_feather(FEATHER_FILE_PATH)
        except ImportError:
            # Fall back to pandas
            df = pd.read_feather(FEATHER_FILE_PATH)
        
        # Check required columns
        required_columns = ["ID1", "ID2", "Similarity"]
        if not all(col in df.columns for col in required_columns):
            logger.error(f"Missing required columns in feather file. Found: {df.columns.tolist()}")
            return False
        
        # Write to CSV with appropriate format for MySQL LOAD DATA INFILE
        logger.info(f"Writing data to CSV: {CSV_FILE_PATH}")
        df.to_csv(CSV_FILE_PATH, index=False, quoting=1)  # quoting=1 for quotes around text fields
        
        logger.info(f"Successfully converted {len(df)} rows to CSV")
        return True
    except Exception as e:
        logger.error(f"Error converting feather to CSV: {e}")
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

def bulk_load_data(connection, driver_type):
    """Bulk load data from CSV file into MySQL table."""
    try:
        cursor = connection.cursor()
        
        # Check if table exists
        try:
            cursor.execute("SELECT COUNT(*) FROM `Sequence_Similarity_Edges_Full`")
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
                cursor.execute("DROP TABLE IF EXISTS `Sequence_Similarity_Edges_Full`")
                connection.commit()
                # Recreate table
                create_table(cursor)
            elif choice == '2':
                logger.info("Truncating table...")
                cursor.execute("TRUNCATE TABLE `Sequence_Similarity_Edges_Full`")
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
        cursor.execute("ALTER TABLE `Sequence_Similarity_Edges_Full` DISABLE KEYS")
        connection.commit()
        
        # Set local_infile to 1
        cursor.execute("SET GLOBAL local_infile = 1")
        
        # Prepare the LOAD DATA INFILE statement with proper path formatting
        csv_path = CSV_FILE_PATH.replace('\\', '/')
        
        load_query = """
        LOAD DATA LOCAL INFILE '{}'
        INTO TABLE `Sequence_Similarity_Edges_Full`
        FIELDS TERMINATED BY ','
        ENCLOSED BY '"'
        LINES TERMINATED BY '\\n'
        IGNORE 1 LINES
        (ID1, ID2, Similarity);
        """.format(csv_path)
        
        # Execute the bulk load
        start_time = time.time()
        logger.info("Executing LOAD DATA INFILE...")
        cursor.execute(load_query)
        connection.commit()
        
        # Enable keys/indexes
        logger.info("Re-enabling keys and building indexes...")
        cursor.execute("ALTER TABLE `Sequence_Similarity_Edges_Full` ENABLE KEYS")
        connection.commit()
        
        # Get final row count
        cursor.execute("SELECT COUNT(*) FROM `Sequence_Similarity_Edges_Full`")
        final_count = cursor.fetchone()[0]
        
        # Log completion
        elapsed_time = time.time() - start_time
        logger.info(f"Bulk load completed in {elapsed_time:.2f} seconds")
        logger.info(f"Loaded {final_count} rows into Sequence_Similarity_Edges_Full")
        
        return True
    except Exception as e:
        logger.error(f"Error during bulk load: {e}")
        return False

def create_table(cursor):
    """Create the table with proper structure and indexes."""
    logger.info("Creating table...")
    
    create_table_query = """
    CREATE TABLE IF NOT EXISTS `Sequence_Similarity_Edges_Full` (
        `ID1` VARCHAR(50) NOT NULL,
        `ID2` VARCHAR(50) NOT NULL,
        `Similarity` FLOAT NOT NULL,
        INDEX `idx_id1_sim_full` (`ID1`, `Similarity`),
        INDEX `idx_id2_sim_full` (`ID2`, `Similarity`)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci
    """
    
    cursor.execute(create_table_query)
    cursor.connection.commit()
    logger.info("Table created successfully")

def try_direct_mysql_command():
    """Try using the mysql command-line client directly."""
    try:
        logger.info("Attempting to use mysql command-line client...")
        
        # Create SQL file with commands
        sql_file_path = "C:/Users/mz30/Downloads/load_data.sql"
        csv_path = CSV_FILE_PATH.replace('\\', '/')
        
        with open(sql_file_path, 'w') as f:
            f.write("""
            -- Create or replace the table
            DROP TABLE IF EXISTS `Sequence_Similarity_Edges_Full`;
            
            CREATE TABLE `Sequence_Similarity_Edges_Full` (
                `ID1` VARCHAR(50) NOT NULL,
                `ID2` VARCHAR(50) NOT NULL,
                `Similarity` FLOAT NOT NULL,
                INDEX `idx_id1_sim_full` (`ID1`, `Similarity`),
                INDEX `idx_id2_sim_full` (`ID2`, `Similarity`)
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
            
            -- Disable keys for faster loading
            ALTER TABLE `Sequence_Similarity_Edges_Full` DISABLE KEYS;
            
            -- Load data
            LOAD DATA LOCAL INFILE '{}' 
            INTO TABLE `Sequence_Similarity_Edges_Full`
            FIELDS TERMINATED BY ','
            ENCLOSED BY '"'
            LINES TERMINATED BY '\\n'
            IGNORE 1 LINES
            (ID1, ID2, Similarity);
            
            -- Re-enable keys
            ALTER TABLE `Sequence_Similarity_Edges_Full` ENABLE KEYS;
            
            -- Show row count
            SELECT COUNT(*) FROM `Sequence_Similarity_Edges_Full`;
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
    
    # Convert feather to CSV
    if not os.path.exists(CSV_FILE_PATH):
        if not convert_feather_to_csv():
            logger.error("Failed to convert feather to CSV. Exiting.")
            return False
    else:
        logger.info(f"Using existing CSV file: {CSV_FILE_PATH}")
    
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
                    os.remove(CSV_FILE_PATH)
                    logger.info(f"Deleted temporary file: {CSV_FILE_PATH}")
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