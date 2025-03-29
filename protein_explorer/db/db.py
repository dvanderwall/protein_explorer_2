"""
Database access module for KinoPlex.

This module provides functions to connect to the database and execute queries against
the phosphosite tables. It includes optimized batch query functionality.
"""

import os
import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Union, Tuple
import sqlalchemy
from sqlalchemy import create_engine, text, select, or_, Table, MetaData, Column
from sqlalchemy.pool import QueuePool
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker, Session
from contextlib import contextmanager
import re
import time
import urllib.parse
import pymysql

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Database connection parameters
INSTANCE_CONNECTION_NAME = os.environ.get('CLOUD_SQL_CONNECTION_NAME', 'future-alcove-454817-e6:us-east4:kinoplex-db')
DB_USER = os.environ.get('DB_USER', 'root')  # Replace with your database user if not root
DB_PASS = os.environ.get('DB_PASS', '@Bismark6')  # Replace with your actual password
DB_NAME = os.environ.get('DB_NAME', 'kinoplex-db')  # This should be the database name, not instance name
DB_HOST = os.environ.get('DB_HOST', '35.245.113.195')  # Use the public IP for direct connections
DB_PORT = os.environ.get('DB_PORT', '3306')

# Table specifications
TABLES = {
    "Phosphosite_Supplementary_Data": {
        "id_column": "PhosphositeID",
        "columns": [
            "GENE", "PROTEIN", "ACC_ID", "MOD_RSD", "DOMAIN", "SITE_+/-7_AA", "MS_LIT",
            "Residue", "Residue_Number", "PhosphositeID", "uniprot_id", "is_known_phosphosite",
            "polar_aa_percent", "nonpolar_aa_percent", "acidic_aa_percent", "basic_aa_percent",
            "site_plddt", "motif_plddt", "nearby_count", "surface_accessibility",
            "DISEASE", "ALTERATION", "Disease_DOMAIN", "DISEASE_PMIDs", "DISEASE_NOTES",
            "Regulatory_Domain", "ON_FUNCTION", "ON_PROCESS", "ON_PROT_INTERACT", "ON_OTHER_INTERACT",
            "REGULATORY_PMIDs", "REGULATORY_NOTES", "SUB_GENE_ID", "KINASE_1", "KINASE_2",
            "KINASE_3", "KINASE_4", "KINASE_5", "IN_VIVO_RXN_1", "IN_VIVO_RXN_2", 
            "IN_VIVO_RXN_3", "IN_VIVO_RXN_4", "IN_VIVO_RXN_5", "IN_VITRO_RXN_1", "IN_VITRO_RXN_2", 
            "IN_VITRO_RXN_3", "IN_VITRO_RXN_4", "IN_VITRO_RXN_5", "StructuralSimAvailable", "HasKnownKinase"
        ]
    },
    "Sequence_Kinase_Scores": {
        "id_column": "node",
        "columns": [
            "node", "label", "AMPKA1", "AMPKA2", "AMPKB1", "AMPKG2", "ARAF", "ASK1", "ATM", "ATR", 
            "Akt1", "Akt2", "Akt3", "AlphaK1", "AurA", "AurB", "AurC", "BMPR1B", "BRAF", "BRSK1 iso2", 
            "BRSK2", "BUB1", "CAMK1A", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G", "CAMK4", "CAMKK1", "CAMKK2", 
            "CASK", "CDC7", "CDK1", "CDK12", "CDK13", "CDK14", "CDK19", "CDK2", "CDK3", "CDK4", 
            "CDK5", "CDK6", "CDK7", "CDK8", "CDK9", "CDKL5", "CK1A", "CK1D", "CK1E", "CK1G1", 
            "CK1G2", "CK2A1", "CK2A2", "CK2B", "CLK1", "CLK2", "ChaK1", "ChaK2", "Chk1", "Chk2", 
            "Cot", "DAPK1", "DAPK3", "DNAPK", "DYRK1A", "DYRK1B", "DYRK2", "EEF2K", "ERK1", "ERK2", 
            "ERK3", "ERK5", "ERK7", "FAM20C", "GCK", "GCN2", "GRK2", "GRK3", "GRK4", "GRK5", 
            "GRK6", "GSK3A", "GSK3B", "HGK", "HIPK2", "HPK1", "IKKA", "IKKB", "IKKE", "ILK", 
            "IRAK1", "IRAK4", "JNK1", "JNK2", "JNK3", "KHS2", "KIS", "LATS1", "LATS2", "LKB1", 
            "LRRK1", "LRRK2", "MAPKAPK2", "MAPKAPK5", "MARK1", "MARK2", "MARK3", "MEK1", "MEKK1", "MEKK3", 
            "MEKK6", "MELK", "MINK", "MKK3", "MKK4", "MKK6", "MKK7", "MLK2", "MLK3", "MRCKA", 
            "MSK1", "MSK2", "MST1", "MST2", "MST3", "MST4", "MYO3A", "Mer", "Mnk1", "Mnk2", 
            "NDR1", "NDR2", "NEK1", "NEK2", "NEK6", "NEK7", "NEK9", "NLK", "Nik", "NuaK1", 
            "OSR1", "P38A", "P38B", "P38D", "P38G", "P70S6KB", "PAK1", "PAK2", "PAK4", "PAK5", 
            "PAK6", "PASK", "PBK", "PDHK1", "PDHK2", "PDHK4", "PDK1", "PGK1", "PIK3CA", "PIK3CG", 
            "PINK1", "PKACA", "PKACB", "PKCA", "PKCB", "PKCB iso2", "PKCD", "PKCE", "PKCG", "PKCH", 
            "PKCI", "PKCT", "PKCZ", "PKG1", "PKG1 iso2", "PKG2", "PKM", "PKN1", "PKR", "PLK1", 
            "PLK2", "PLK3", "PLK4", "PRKD1", "PRKD2", "PRKD3", "Pim1", "Pim2", "Pim3", "Pyk2", 
            "QIK", "RAF1", "RIPK1", "RIPK2", "RIPK3", "ROCK1", "ROCK2", "RSK2", "RSK3", "SGK1", 
            "SGK3", "SLK", "SMG1", "SRPK1", "STLK3", "TAK1", "TAO1", "TAO2", "TBK1", "TGFBR1", 
            "TGFBR2", "TLK1", "TNIK", "TTBK2", "TTK", "ULK1", "ULK2", "ULK3", "VRK1", "VRK2", 
            "WNK1", "WNK4", "YSK1", "ZAK", "mTOR", "p70S6K", "p90RSK", "smMLCK"
        ]
    },
    "Structure_Kinase_Scores": {
        "id_column": "node",
        "columns": [
            "node", "label", "AMPKA1", "AMPKA2", "AMPKB1", "AMPKG2", "ARAF", "ASK1", "ATM", "ATR", 
            "Akt1", "Akt2", "Akt3", "AlphaK1", "AurA", "AurB", "AurC", "BMPR1B", "BRAF", "BRSK1 iso2", 
            "BRSK2", "BUB1", "CAMK1A", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G", "CAMK4", "CAMKK1", "CAMKK2", 
            "CASK", "CDC7", "CDK1", "CDK12", "CDK13", "CDK14", "CDK19", "CDK2", "CDK3", "CDK4", 
            "CDK5", "CDK6", "CDK7", "CDK8", "CDK9", "CDKL5", "CK1A", "CK1D", "CK1E", "CK1G1", 
            "CK1G2", "CK2A1", "CK2A2", "CK2B", "CLK1", "CLK2", "ChaK1", "ChaK2", "Chk1", "Chk2", 
            "Cot", "DAPK1", "DAPK3", "DNAPK", "DYRK1A", "DYRK1B", "DYRK2", "EEF2K", "ERK1", "ERK2", 
            "ERK3", "ERK5", "ERK7", "FAM20C", "GCK", "GCN2", "GRK2", "GRK3", "GRK4", "GRK5", 
            "GRK6", "GSK3A", "GSK3B", "HGK", "HIPK2", "HPK1", "IKKA", "IKKB", "IKKE", "ILK", 
            "IRAK1", "IRAK4", "JNK1", "JNK2", "JNK3", "KHS2", "KIS", "LATS1", "LATS2", "LKB1", 
            "LRRK1", "LRRK2", "MAPKAPK2", "MAPKAPK5", "MARK1", "MARK2", "MARK3", "MEK1", "MEKK1", "MEKK3", 
            "MEKK6", "MELK", "MINK", "MKK3", "MKK4", "MKK6", "MKK7", "MLK2", "MLK3", "MRCKA", 
            "MSK1", "MSK2", "MST1", "MST2", "MST3", "MST4", "MYO3A", "Mer", "Mnk1", "Mnk2", 
            "NDR1", "NDR2", "NEK1", "NEK2", "NEK6", "NEK7", "NEK9", "NLK", "Nik", "NuaK1", 
            "OSR1", "P38A", "P38B", "P38D", "P38G", "P70S6KB", "PAK1", "PAK2", "PAK4", "PAK5", 
            "PAK6", "PASK", "PBK", "PDHK1", "PDHK2", "PDHK4", "PDK1", "PGK1", "PIK3CA", "PIK3CG", 
            "PINK1", "PKACA", "PKACB", "PKCA", "PKCB", "PKCB iso2", "PKCD", "PKCE", "PKCG", "PKCH", 
            "PKCI", "PKCT", "PKCZ", "PKG1", "PKG1 iso2", "PKG2", "PKM", "PKN1", "PKR", "PLK1", 
            "PLK2", "PLK3", "PLK4", "PRKD1", "PRKD2", "PRKD3", "Pim1", "Pim2", "Pim3", "Pyk2", 
            "QIK", "RAF1", "RIPK1", "RIPK2", "RIPK3", "ROCK1", "ROCK2", "RSK2", "RSK3", "SGK1", 
            "SGK3", "SLK", "SMG1", "SRPK1", "STLK3", "TAK1", "TAO1", "TAO2", "TBK1", "TGFBR1", 
            "TGFBR2", "TLK1", "TNIK", "TTBK2", "TTK", "ULK1", "ULK2", "ULK3", "VRK1", "VRK2", 
            "WNK1", "WNK4", "YSK1", "ZAK", "mTOR", "p70S6K", "p90RSK", "smMLCK"
        ]
    },
    "Structural_Similarity_Edges": {
        "columns": [
            "Query", "Target", "Mapped_Site", "RMSD", "Mapped_Residues", 
            "FoldDisco_Score", "Max_Residues", "NumberMapped", "Map_Ratio"
        ]
    },
    "Sequence_Similarity_Edges_7": {
        "columns": [
            "ID1", "ID2", "Similarity"
        ]
    }
}

# Query templates for MySQL
QUERY_TEMPLATES = {
    "get_phosphosite_data": """
        SELECT * FROM `Phosphosite_Supplementary_Data`
        WHERE `PhosphositeID` = :site_id
    """,
    
    "get_phosphosites_batch": """
        SELECT * FROM `Phosphosite_Supplementary_Data`
        WHERE `PhosphositeID` IN :site_ids
    """,
    
    "find_structural_matches": """
        SELECT * FROM `Structural_Similarity_Edges`
        WHERE `Query` = :site_id AND `RMSD` <= :rmsd_threshold
        ORDER BY `RMSD` ASC
    """,
    
    "find_structural_matches_batch": """
        SELECT * FROM `Structural_Similarity_Edges`
        WHERE `Query` IN :site_ids AND `RMSD` <= :rmsd_threshold
        ORDER BY `RMSD` ASC
    """,
    
    "find_sequence_matches": """
        SELECT * FROM `Sequence_Similarity_Edges_7`
        WHERE (`ID1` = :site_id OR `ID2` = :site_id) AND `Similarity` >= :min_similarity
        ORDER BY `Similarity` DESC
    """,
    
    "find_sequence_matches_batch": """
        SELECT * FROM `Sequence_Similarity_Edges_7`
        WHERE (`ID1` IN :site_ids OR `ID2` IN :site_ids) AND `Similarity` >= :min_similarity
        ORDER BY `Similarity` DESC
    """,
    
    "get_structure_kinase_scores": """
        SELECT * FROM `Structure_Kinase_Scores`
        WHERE `node` = :site_id
    """,
    
    "get_sequence_kinase_scores": """
        SELECT * FROM `Sequence_Kinase_Scores`
        WHERE `node` = :site_id
    """,
    
    "get_structure_kinase_scores_batch": """
        SELECT * FROM `Structure_Kinase_Scores`
        WHERE `node` IN :site_ids
    """,
    
    "get_sequence_kinase_scores_batch": """
        SELECT * FROM `Sequence_Kinase_Scores`
        WHERE `node` IN :site_ids
    """,
    
    # Additional utility queries
    "check_database_connection": """
        SELECT 1
    """,
    
    "get_table_names": """
        SHOW TABLES
    """,
    
    "get_column_names": """
        SHOW COLUMNS FROM :table_name
    """,
    
    "count_phosphosites": """
        SELECT COUNT(*) as count FROM `Phosphosite_Supplementary_Data`
    """,
    
    "count_structural_matches": """
        SELECT COUNT(*) as count FROM `Structural_Similarity_Edges`
    """,
    
    "count_sequence_matches": """
        SELECT COUNT(*) as count FROM `Sequence_Similarity_Edges_7`
    """,
    
    "get_kinases_list": """
        SELECT COLUMN_NAME
        FROM INFORMATION_SCHEMA.COLUMNS
        WHERE TABLE_NAME = 'Structure_Kinase_Scores'
        AND COLUMN_NAME NOT IN ('node', 'label')
        ORDER BY COLUMN_NAME
    """,
    
    "get_known_kinase": """
        SELECT KINASE_1, KINASE_2, KINASE_3, KINASE_4, KINASE_5
        FROM `Phosphosite_Supplementary_Data`
        WHERE `PhosphositeID` = :site_id
    """,
    
    "find_sites_by_protein": """
        SELECT * FROM `Phosphosite_Supplementary_Data`
        WHERE `uniprot_id` = :uniprot_id
    """,
    
    "find_sites_by_gene": """
        SELECT * FROM `Phosphosite_Supplementary_Data`
        WHERE `GENE` = :gene_symbol
    """,
    
    "get_top_similar_sites": """
        SELECT `ID1`, `ID2`, `Similarity` 
        FROM `Sequence_Similarity_Edges_7`
        WHERE (`ID1` = :site_id OR `ID2` = :site_id)
        AND `Similarity` >= :min_similarity
        ORDER BY `Similarity` DESC
        LIMIT :limit
    """,
    
    "get_top_structural_matches": """
        SELECT * 
        FROM `Structural_Similarity_Edges`
        WHERE `Query` = :site_id
        ORDER BY `RMSD` ASC
        LIMIT :limit
    """,
    
    "get_top_kinases": """
        SELECT :kinase_columns
        FROM `:score_table`
        WHERE `node` = :site_id
    """
}

# Database connection URL
DB_PASS = urllib.parse.quote_plus(os.environ.get("DB_PASS", "@Bismark6"))  # URL encode the password
DB_NAME = os.environ.get("DB_NAME", "kinoplex-db")

# Database connection string for MySQL
DB_URL = f"mysql+pymysql://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}"

# Global engine and session maker objects
engine = None
Session = None
metadata = MetaData()

# Query cache for performance
QUERY_CACHE = {}
CACHE_EXPIRY = 3600  # Cache expiry in seconds (1 hour)
CACHE_TIMESTAMPS = {}  # Track when items were added to cache

def init_db():
    """Initialize database connection."""
    global engine, Session, metadata
    try:
        engine = create_engine(
            DB_URL,
            poolclass=QueuePool,
            pool_size=10,
            max_overflow=20,
            pool_timeout=30,
            pool_recycle=1800
        )
        # Create session maker
        Session = sessionmaker(bind=engine)
        
        # Test connection
        with engine.connect() as conn:
            conn.execute(text("SELECT 1"))
            
        # Set up metadata for reflection
        metadata.reflect(engine)
        
        logger.info("Database connection initialized successfully")
        return True
    except Exception as e:
        logger.error(f"Error initializing database: {e}")
        return False

@contextmanager
def get_db_session():
    """
    Context manager for database sessions.
    Ensures connections are properly closed after use.
    
    Yields:
        SQLAlchemy session
    """
    global Session
    if Session is None:
        init_db()
    
    if Session is None:
        raise ValueError("Failed to initialize database session")
    
    session = Session()
    try:
        yield session
        session.commit()
    except Exception as e:
        session.rollback()
        raise e
    finally:
        session.close()

def execute_query(query: str, params: dict = None) -> pd.DataFrame:
    """
    Execute a SQL query and return results as a DataFrame.
    Handles MySQL's parameter binding quirks.
    
    Args:
        query: SQL query string
        params: Parameters for the query
        
    Returns:
        Pandas DataFrame with query results
    """
    global engine
    
    if engine is None:
        init_db()
        
    if engine is None:
        logger.error("Database engine not initialized")
        return pd.DataFrame()
    
    try:
        # Process parameters for MySQL
        processed_params = {}
        if params:
            for key, value in params.items():
                # Handle tuple parameters for IN clauses
                if isinstance(value, tuple) or isinstance(value, list):
                    # For IN clauses, convert tuple/list to comma-separated string
                    if len(value) > 0:
                        # For string values, add quotes
                        if isinstance(value[0], str):
                            formatted_values = ", ".join([f"'{v}'" for v in value])
                        else:
                            formatted_values = ", ".join([str(v) for v in value])
                        # Replace the parameter in the query directly
                        query = query.replace(f":{key}", f"({formatted_values})")
                    else:
                        # Empty tuple/list - replace with NULL to avoid SQL error
                        query = query.replace(f":{key}", "(NULL)")
                else:
                    # Regular parameters
                    processed_params[key] = value
        
        with engine.connect() as conn:
            # Use text() to ensure proper SQL parameter handling
            result = conn.execute(text(query), processed_params or {})
            df = pd.DataFrame(result.fetchall())
            if not df.empty:
                df.columns = result.keys()
            return df
    except SQLAlchemyError as e:
        logger.error(f"Database query error: {e}")
        logger.error(f"Query was: {query}")
        logger.error(f"Parameters were: {params}")
        return pd.DataFrame()
    except Exception as e:
        logger.error(f"Error executing query: {e}")
        return pd.DataFrame()
    
def execute_batch_query(query_template: str, id_list: List[str], batch_size: int = 500, 
                      id_param_name: str = "site_ids", extra_params: dict = None) -> pd.DataFrame:
    """
    Execute a batch query by breaking up a large list into smaller batches.
    Addresses MySQL's limitations with large IN clauses.
    
    Args:
        query_template: SQL query template with :id_param_name placeholder
        id_list: List of IDs to process in batches
        batch_size: Maximum number of IDs per batch
        id_param_name: Name of the ID parameter in the query
        extra_params: Additional query parameters
        
    Returns:
        Combined DataFrame with results from all batches
    """
    if not id_list:
        return pd.DataFrame()
    
    # Initialize an empty result DataFrame
    combined_results = pd.DataFrame()
    
    # Process in batches
    for i in range(0, len(id_list), batch_size):
        batch = id_list[i:i + batch_size]
        logger.info(f"Processing batch {i//batch_size + 1} with {len(batch)} items")
        
        # Prepare parameters
        params = {id_param_name: tuple(batch)}
        if extra_params:
            params.update(extra_params)
        
        # Execute query for this batch
        batch_results = execute_query(query_template, params)
        
        # Append results
        if not batch_results.empty:
            if combined_results.empty:
                combined_results = batch_results
            else:
                combined_results = pd.concat([combined_results, batch_results], ignore_index=True)
    
    return combined_results
    

def is_cache_valid(cache_key: str) -> bool:
    """Check if a cache entry is still valid based on timestamp."""
    timestamp = CACHE_TIMESTAMPS.get(cache_key)
    if timestamp is None:
        return False
    
    current_time = time.time()
    return (current_time - timestamp) < CACHE_EXPIRY

def get_phosphosite_data(site_id: str) -> Optional[Dict]:
    """
    Get phosphosite data from the database.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary with phosphosite data or None if not found
    """
    # Check cache first
    cache_key = f"phosphosite_data_{site_id}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        return QUERY_CACHE[cache_key]
    
    try:
        query = QUERY_TEMPLATES["get_phosphosite_data"]
        df = execute_query(query, {"site_id": site_id})
        
        if df.empty:
            logger.warning(f"No phosphosite data found for {site_id}")
            # Cache the negative result too
            QUERY_CACHE[cache_key] = None
            CACHE_TIMESTAMPS[cache_key] = time.time()
            return None
            
        # Convert first row to dictionary
        result = df.iloc[0].to_dict()
        
        # Cache the result
        QUERY_CACHE[cache_key] = result
        CACHE_TIMESTAMPS[cache_key] = time.time()
        
        return result
    except Exception as e:
        logger.error(f"Error getting phosphosite data: {e}")
        return None

def get_phosphosites_batch(site_ids: List[str]) -> Dict[str, Dict]:
    """
    Get phosphosite data for multiple sites in a batch.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary mapping site IDs to phosphosite data dictionaries
    """
    if not site_ids:
        return {}
        
    try:
        # Use the batch helper function
        query = QUERY_TEMPLATES["get_phosphosites_batch"]
        df = execute_batch_query(query, site_ids, batch_size=500, id_param_name="site_ids")
        
        if df.empty:
            logger.warning(f"No phosphosite data found for batch of {len(site_ids)} sites")
            return {}
            
        # Convert to dictionary of dictionaries
        result = {}
        for _, row in df.iterrows():
            row_dict = row.to_dict()
            if "PhosphositeID" in row_dict:
                site_id = row_dict["PhosphositeID"]
                result[site_id] = row_dict
                
        logger.info(f"Retrieved data for {len(result)} out of {len(site_ids)} phosphosites")
        return result
    except Exception as e:
        logger.error(f"Error getting batch phosphosite data: {e}")
        return {}

def find_structural_matches(site_id: str, rmsd_threshold: float = 5.0) -> List[Dict]:
    """
    Find structural matches for a site.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        rmsd_threshold: Maximum RMSD value for matches
        
    Returns:
        List of dictionaries with match information
    """
    # Check cache first
    cache_key = f"structural_matches_{site_id}_{rmsd_threshold}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        return QUERY_CACHE[cache_key]
    
    try:
        query = QUERY_TEMPLATES["find_structural_matches"]
        df = execute_query(query, {
            "site_id": site_id,
            "rmsd_threshold": rmsd_threshold
        })
        
        matches = []
        if not df.empty:
            # Process matches
            for _, row in df.iterrows():
                row_dict = row.to_dict()
                
                # Parse target info
                target_id = row_dict.get("Target", "")
                target_parts = target_id.split('_')
                
                if len(target_parts) >= 2:
                    target_uniprot = target_parts[0]
                    target_site_raw = target_parts[1]
                    
                    # Extract residue type and number
                    match = re.match(r'([STY]?)(\d+)', target_site_raw)
                    if match:
                        residue_type = match.group(1) or "S"  # Default to S if missing
                        residue_number = match.group(2)
                        target_site = f"{residue_type}{residue_number}"
                    else:
                        target_site = target_site_raw
                        residue_type = ""
                    
                    # Create match dictionary
                    match_dict = {
                        "query_id": site_id,
                        "target_id": target_id,
                        "target_uniprot": target_uniprot,
                        "target_site": target_site,
                        "rmsd": float(row_dict.get("RMSD", 0)),
                        "max_residues": int(row_dict.get("Max_Residues", 0)),
                        "number_mapped": int(row_dict.get("NumberMapped", 0)),
                        "map_ratio": float(row_dict.get("Map_Ratio", 0))
                    }
                    
                    matches.append(match_dict)
        
        # Cache the result
        QUERY_CACHE[cache_key] = matches
        CACHE_TIMESTAMPS[cache_key] = time.time()
        
        return matches
    except Exception as e:
        logger.error(f"Error finding structural matches: {e}")
        return []

def find_structural_matches_batch(site_ids: List[str], rmsd_threshold: float = 5.0) -> Dict[str, List[Dict]]:
    """
    Find structural matches for multiple sites in a batch.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        rmsd_threshold: Maximum RMSD value for matches
        
    Returns:
        Dictionary mapping site IDs to lists of match dictionaries
    """
    if not site_ids:
        return {}
        
    try:
        # Use the batch helper function
        query = QUERY_TEMPLATES["find_structural_matches_batch"]
        df = execute_batch_query(
            query, 
            site_ids, 
            batch_size=200,  # Smaller batch size for potentially larger result sets
            id_param_name="site_ids", 
            extra_params={"rmsd_threshold": rmsd_threshold}
        )
        
        if df.empty:
            logger.warning(f"No structural matches found for batch of {len(site_ids)} sites")
            return {}
            
        # Organize matches by query site
        results = {}
        for _, row in df.iterrows():
            row_dict = row.to_dict()
            query_id = row_dict.get("Query")
            
            if query_id not in results:
                results[query_id] = []
                
            # Process match data as before...
            # (code from previous implementation)
            
            # Parse target info
            target_id = row_dict.get("Target", "")
            target_parts = target_id.split('_')
            
            if len(target_parts) >= 2:
                target_uniprot = target_parts[0]
                target_site_raw = target_parts[1]
                
                # Extract residue type and number
                match = re.match(r'([STY]?)(\d+)', target_site_raw)
                if match:
                    residue_type = match.group(1) or "S"  # Default to S if missing
                    residue_number = match.group(2)
                    target_site = f"{residue_type}{residue_number}"
                else:
                    target_site = target_site_raw
                    residue_type = ""
                
                # Create match dictionary
                match_dict = {
                    "query_id": query_id,
                    "target_id": target_id,
                    "target_uniprot": target_uniprot,
                    "target_site": target_site,
                    "rmsd": float(row_dict.get("RMSD", 0)),
                    "max_residues": int(row_dict.get("Max_Residues", 0)),
                    "number_mapped": int(row_dict.get("NumberMapped", 0)),
                    "map_ratio": float(row_dict.get("Map_Ratio", 0))
                }
                
                results[query_id].append(match_dict)
                
        logger.info(f"Retrieved structural matches for {len(results)} out of {len(site_ids)} sites")
        return results
    except Exception as e:
        logger.error(f"Error finding batch structural matches: {e}")
        return {}


def find_sequence_matches(site_id: str, min_similarity: float = 0.4) -> List[Dict]:
    """
    Find sequence similarity matches for a site.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        min_similarity: Minimum similarity score to include (0-1)
        
    Returns:
        List of dictionaries with match information
    """
    # Check cache first
    cache_key = f"sequence_matches_{site_id}_{min_similarity}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        return QUERY_CACHE[cache_key]
    
    try:
        query = QUERY_TEMPLATES["find_sequence_matches"]
        df = execute_query(query, {
            "site_id": site_id,
            "min_similarity": min_similarity
        })
        
        matches = []
        if not df.empty:
            # Process matches
            for _, row in df.iterrows():
                row_dict = row.to_dict()
                id1 = row_dict.get("ID1")
                id2 = row_dict.get("ID2")
                similarity = float(row_dict.get("Similarity", 0))
                
                # Determine which ID is the target (the one that isn't the query)
                if id1 == site_id:
                    target_id = id2
                else:
                    target_id = id1
                    
                # Parse target info
                target_parts = target_id.split('_')
                if len(target_parts) >= 2:
                    target_uniprot = target_parts[0]
                    target_site_raw = target_parts[1]
                    
                    # Extract residue type and number if possible
                    match = re.match(r'([STY]?)(\d+)', target_site_raw)
                    if match:
                        residue_type = match.group(1) or None
                        residue_number = match.group(2)
                        target_site = f"{residue_type or ''}{residue_number}"
                    else:
                        target_site = target_site_raw
                        residue_type = None
                    
                    # Create match dictionary
                    match_dict = {
                        "query_id": site_id,
                        "target_id": target_id,
                        "target_uniprot": target_uniprot,
                        "target_site": target_site,
                        "site_type": residue_type,
                        "similarity": similarity
                    }
                    
                    matches.append(match_dict)
        
        # Cache the result
        QUERY_CACHE[cache_key] = matches
        CACHE_TIMESTAMPS[cache_key] = time.time()
        
        return matches
    except Exception as e:
        logger.error(f"Error finding sequence matches: {e}")
        return []

def find_sequence_matches_batch(site_ids: List[str], min_similarity: float = 0.4) -> Dict[str, List[Dict]]:
    """
    Find sequence similarity matches for multiple sites in a batch.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        min_similarity: Minimum similarity score to include (0-1)
        
    Returns:
        Dictionary mapping site IDs to lists of match dictionaries
    """
    if not site_ids:
        return {}
        
    try:
        # Use the batch helper function
        query = QUERY_TEMPLATES["find_sequence_matches_batch"]
        df = execute_batch_query(
            query, 
            site_ids, 
            batch_size=100,  # Even smaller batch size for potentially very large result sets
            id_param_name="site_ids", 
            extra_params={"min_similarity": min_similarity}
        )
        
        if df.empty:
            logger.warning(f"No sequence matches found for batch of {len(site_ids)} sites")
            return {}
            
        # Process results as before...
        # (remaining implementation)
        
        # Organize matches by query site
        results = {}
        for _, row in df.iterrows():
            row_dict = row.to_dict()
            id1 = row_dict.get("ID1")
            id2 = row_dict.get("ID2")
            similarity = float(row_dict.get("Similarity", 0))
            
            # Determine which ID is the query and which is the target
            if id1 in site_ids:
                query_id = id1
                target_id = id2
            else:
                query_id = id2
                target_id = id1
                
            # Initialize results for this query if needed
            if query_id not in results:
                results[query_id] = []
                
            # Parse target info
            target_parts = target_id.split('_')
            if len(target_parts) >= 2:
                target_uniprot = target_parts[0]
                target_site_raw = target_parts[1]
                
                # Extract residue type and number if possible
                match = re.match(r'([STY]?)(\d+)', target_site_raw)
                if match:
                    residue_type = match.group(1) or None
                    residue_number = match.group(2)
                    target_site = f"{residue_type or ''}{residue_number}"
                else:
                    target_site = target_site_raw
                    residue_type = None
                
                # Create match dictionary
                match_dict = {
                    "query_id": query_id,
                    "target_id": target_id,
                    "target_uniprot": target_uniprot,
                    "target_site": target_site,
                    "site_type": residue_type,
                    "similarity": similarity
                }
                
                results[query_id].append(match_dict)
                
        logger.info(f"Retrieved sequence matches for {len(results)} out of {len(site_ids)} sites")
        return results
    except Exception as e:
        logger.error(f"Error finding batch sequence matches: {e}")
        return {}

def get_kinase_scores(site_id: str, score_type: str = 'structure') -> Dict:
    """
    Get kinase scores for a specific site.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary with kinase names as keys and scores as values
    """
    # Check cache first
    cache_key = f"kinase_scores_{score_type}_{site_id}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        return QUERY_CACHE[cache_key]
    
    try:
        # Select appropriate table based on score type
        if score_type.lower() == 'structure':
            query = QUERY_TEMPLATES["get_structure_kinase_scores"]
        else:
            query = QUERY_TEMPLATES["get_sequence_kinase_scores"]
            
        df = execute_query(query, {"site_id": site_id})
        
        if df.empty:
            logger.warning(f"No {score_type} kinase scores found for {site_id}")
            # Cache empty result
            empty_result = {'site_id': site_id, 'scores': {}, 'known_kinase': None}
            QUERY_CACHE[cache_key] = empty_result
            CACHE_TIMESTAMPS[cache_key] = time.time()
            return empty_result
            
        # Extract scores from first row
        row = df.iloc[0]
        
        # Extract all kinase columns and their scores
        # Skip non-kinase columns (node and label)
        scores = {}
        known_kinase = None
        
        for column in df.columns:
            if column == 'label':
                known_kinase = row[column] if pd.notna(row[column]) else None
            elif column != 'node':
                try:
                    value = row[column]
                    if pd.notna(value):  # Check if not NaN
                        scores[column] = float(value)
                    else:
                        scores[column] = 0.0
                except (ValueError, TypeError):
                    logger.warning(f"Non-numeric value for kinase {column}: {row[column]}")
        
        # Create result
        result = {
            'site_id': site_id,
            'scores': scores,
            'known_kinase': known_kinase
        }
        
        # Cache result
        QUERY_CACHE[cache_key] = result
        CACHE_TIMESTAMPS[cache_key] = time.time()
        
        return result
    except Exception as e:
        logger.error(f"Error getting kinase scores: {e}")
        return {'site_id': site_id, 'scores': {}, 'known_kinase': None}

def get_kinase_scores_batch(site_ids: List[str], score_type: str = 'structure') -> Dict[str, Dict]:
    """
    Get kinase scores for multiple sites in a batch.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary mapping site IDs to dictionaries of kinase scores
    """
    if not site_ids:
        return {}
        
    try:
        # Select appropriate table based on score type
        if score_type.lower() == 'structure':
            query = QUERY_TEMPLATES["get_structure_kinase_scores_batch"]
        else:
            query = QUERY_TEMPLATES["get_sequence_kinase_scores_batch"]
            
        # Use the batch helper function
        df = execute_batch_query(query, site_ids)
        
        if df.empty:
            logger.warning(f"No {score_type} kinase scores found for batch of {len(site_ids)} sites")
            return {}
            
        # Process results
        results = {}
        for _, row in df.iterrows():
            site_id = row['node']
            
            # Extract all kinase columns and their scores
            # Skip non-kinase columns (node and label)
            scores = {}
            for column in df.columns:
                if column not in ['node', 'label']:
                    scores[column] = float(row[column])
                    
            # Check if any known kinases are available in the database
            # (For future enhancement - check Phosphosite_Supplementary_Data)
            known_kinase = None
            
            results[site_id] = {
                'site_id': site_id,
                'scores': scores,
                'known_kinase': known_kinase
            }
            
        logger.info(f"Retrieved kinase scores for {len(results)} out of {len(site_ids)} sites")
        return results
    except Exception as e:
        logger.error(f"Error getting batch kinase scores: {e}")
        return {}

def clear_cache():
    """Clear the query cache."""
    global QUERY_CACHE, CACHE_TIMESTAMPS
    QUERY_CACHE = {}
    CACHE_TIMESTAMPS = {}
    logger.info("Query cache cleared")

def get_db_stats():
    """
    Get database statistics.
    
    Returns:
        Dictionary with database statistics
    """
    stats = {
        'table_counts': {},
        'cache_size': len(QUERY_CACHE),
        'engine_status': engine is not None,
        'tables_available': list(TABLES.keys())
    }
    
    try:
        # Execute count queries for each table
        for table_name in TABLES:
            query = f'SELECT COUNT(*) FROM "{table_name}"'
            df = execute_query(query)
            if not df.empty:
                stats['table_counts'][table_name] = int(df.iloc[0, 0])
            else:
                stats['table_counts'][table_name] = 0
    except Exception as e:
        logger.error(f"Error getting database stats: {e}")
    
    return stats

# Initialize database on module import
init_db()