import os
import logging
import pymysql
import time

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Database connection parameters
DB_HOST = os.environ.get('DB_HOST', '35.245.113.195')
DB_USER = os.environ.get('DB_USER', 'root')
DB_PASS = os.environ.get('DB_PASS', '@Bismark6')
DB_NAME = os.environ.get('DB_NAME', 'kinoplex-db')
DB_PORT = int(os.environ.get('DB_PORT', '3306'))

# List of index creation statements with prefix lengths for TEXT/BLOB columns
# For text columns, typically using 50-255 characters works well for indexing
INDEX_STATEMENTS = [
    # Phosphosite table
    "CREATE INDEX idx_phosphositeID ON Phosphosite_Supplementary_Data(PhosphositeID(50))",
    "CREATE INDEX idx_uniprot_id ON Phosphosite_Supplementary_Data(uniprot_id(20))",
    "CREATE INDEX idx_gene ON Phosphosite_Supplementary_Data(GENE(20))",
    
    # Structural similarity table
    "CREATE INDEX idx_query_rmsd ON Structural_Similarity_Edges(Query(50), RMSD)",
    "CREATE INDEX idx_target ON Structural_Similarity_Edges(Target(50))",
    
    # Sequence similarity table
    "CREATE INDEX idx_id1_sim ON Sequence_Similarity_Edges_7(ID1(50), Similarity)",
    "CREATE INDEX idx_id2_sim ON Sequence_Similarity_Edges_7(ID2(50), Similarity)",
    
    # Kinase scores tables
    "CREATE INDEX idx_node_structure ON Structure_Kinase_Scores(node(50))",
    "CREATE INDEX idx_node_sequence ON Sequence_Kinase_Scores(node(50))"
]

# Verification queries
VERIFICATION_QUERIES = [
    "SHOW INDEX FROM Phosphosite_Supplementary_Data",
    "SHOW INDEX FROM Structural_Similarity_Edges",
    "SHOW INDEX FROM Sequence_Similarity_Edges_7",
    "SHOW INDEX FROM Structure_Kinase_Scores",
    "SHOW INDEX FROM Sequence_Kinase_Scores"
]

def create_indexes():
    """Connect to the database and create indexes with error handling for TEXT/BLOB columns"""
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
        
        # Create cursor
        with connection.cursor() as cursor:
            # Execute each index creation statement
            for statement in INDEX_STATEMENTS:
                try:
                    logger.info(f"Executing: {statement}")
                    start_time = time.time()
                    cursor.execute(statement)
                    end_time = time.time()
                    logger.info(f"Index created in {end_time - start_time:.2f} seconds")
                except pymysql.err.IntegrityError as e:
                    # Index already exists - this is fine
                    logger.info(f"Index already exists: {e}")
                except pymysql.err.OperationalError as e:
                    # Handle other errors
                    if "Duplicate key name" in str(e):
                        logger.info(f"Index already exists: {e}")
                    else:
                        logger.error(f"Error creating index with '{statement}': {e}")
            
            # Commit changes
            connection.commit()
            logger.info("All indexes created and committed successfully")
            
            # Verify indexes
            logger.info("Verifying indexes...")
            for query in VERIFICATION_QUERIES:
                table_name = query.split(" ")[-1]
                logger.info(f"Checking indexes for table: {table_name}")
                cursor.execute(query)
                results = cursor.fetchall()
                index_count = len(results)
                logger.info(f"Found {index_count} indexes for {table_name}")
                
                # Print index details
                if index_count > 0:
                    logger.info(f"Index details for {table_name}:")
                    for result in results:
                        # Extract relevant index information based on MySQL's SHOW INDEX format
                        try:
                            index_name = result[2]  # Key_name is usually in the 3rd position
                            column_name = result[4]  # Column_name is usually in the 5th position
                            logger.info(f"  - Index: {index_name}, Column: {column_name}")
                        except IndexError:
                            logger.info(f"  - {result}")  # Fallback if structure is different
        
        logger.info("Index verification complete!")
        
    except Exception as e:
        logger.error(f"Error creating indexes: {e}")
        if connection:
            connection.rollback()
    finally:
        if connection:
            connection.close()
            logger.info("Database connection closed")

if __name__ == "__main__":
    create_indexes()