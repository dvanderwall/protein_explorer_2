"""
Database access module for KinoPlex.

This module provides functions to connect to the database and execute queries against
the phosphosite tables. It includes optimized batch query functionality.
"""

import os
import logging
import pandas as pd
from typing import Dict, List, Optional, Union, Tuple
import sqlalchemy
from sqlalchemy import create_engine, text, Table, MetaData, select, and_, or_
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import QueuePool
from contextlib import contextmanager
import numpy as np

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Global database engine and connection
db_engine = None
Session = None
metadata = MetaData()

# Cache for query results to minimize database calls
QUERY_CACHE = {}

def init_db(connection_string=None):
    """
    Initialize database connection.
    
    Args:
        connection_string: Database connection string. If None, it will be read from environment variables.
    
    Returns:
        sqlalchemy.engine: The database engine
    """
    global db_engine, Session, metadata
    
    if db_engine is not None:
        logger.info("Database already initialized")
        return db_engine
    
    # Get connection string from environment if not provided
    if connection_string is None:
        connection_string = os.environ.get('DATABASE_URL')
        
    if connection_string is None:
        raise ValueError("No database connection string provided. Set DATABASE_URL environment variable or pass connection_string parameter.")
    
    logger.info(f"Initializing database connection")
    
    # Create engine with connection pooling
    db_engine = create_engine(
        connection_string, 
        poolclass=QueuePool,
        pool_size=10,
        max_overflow=20,
        pool_timeout=30,
        pool_recycle=3600,
        pool_pre_ping=True
    )
    
    # Create session factory
    Session = sessionmaker(bind=db_engine)
    
    # Reflect the database tables
    try:
        logger.info("Reflecting database tables")
        metadata.reflect(bind=db_engine)
        logger.info(f"Successfully reflected {len(metadata.tables)} tables")
    except Exception as e:
        logger.error(f"Error reflecting database tables: {e}")
        raise
    
    return db_engine

@contextmanager
def get_db_session():
    """
    Context manager for database sessions.
    Ensures connections are properly closed after use.
    
    Yields:
        SQLAlchemy session
    """
    if Session is None:
        init_db()
    
    session = Session()
    try:
        yield session
        session.commit()
    except Exception as e:
        session.rollback()
        raise e
    finally:
        session.close()

def get_phosphosite_data(site_id: str) -> Optional[Dict]:
    """
    Get supplementary data for a specific phosphosite from the database.
    
    Args:
        site_id: The site ID in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary with supplementary data or None if not found
    """
    # Check cache first
    cache_key = f"phosphosite_data_{site_id}"
    if cache_key in QUERY_CACHE:
        return QUERY_CACHE[cache_key]
    
    try:
        with get_db_session() as session:
            # Get the Phosphosite_Supplementary_Data table
            phosphosite_table = metadata.tables.get('Phosphosite_Supplementary_Data')
            
            if not phosphosite_table:
                logger.error("Phosphosite_Supplementary_Data table not found in database")
                return None
            
            # Create query
            query = select([phosphosite_table]).where(phosphosite_table.c.PhosphositeID == site_id)
            
            # Execute query
            result = session.execute(query).fetchone()
            
            if result:
                # Convert result to dictionary
                result_dict = {column: getattr(result, column) for column in result._fields}
                
                # Cache result
                QUERY_CACHE[cache_key] = result_dict
                return result_dict
            
            return None
    except Exception as e:
        logger.error(f"Error retrieving phosphosite data for {site_id}: {e}")
        return None

def get_phosphosites_batch(site_ids: List[str]) -> Dict[str, Dict]:
    """
    Get supplementary data for multiple phosphosites in a single query.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary mapping site IDs to their data dictionaries
    """
    result_dict = {}
    
    # Return empty dict if no site_ids
    if not site_ids:
        return result_dict
    
    # Check cache first for all site_ids
    uncached_site_ids = []
    for site_id in site_ids:
        cache_key = f"phosphosite_data_{site_id}"
        if cache_key in QUERY_CACHE:
            result_dict[site_id] = QUERY_CACHE[cache_key]
        else:
            uncached_site_ids.append(site_id)
    
    # If all site_ids were in cache, return directly
    if not uncached_site_ids:
        return result_dict
    
    try:
        with get_db_session() as session:
            # Get the Phosphosite_Supplementary_Data table
            phosphosite_table = metadata.tables.get('Phosphosite_Supplementary_Data')
            
            if not phosphosite_table:
                logger.error("Phosphosite_Supplementary_Data table not found in database")
                return result_dict
            
            # Create query for all uncached site_ids
            query = select([phosphosite_table]).where(phosphosite_table.c.PhosphositeID.in_(uncached_site_ids))
            
            # Execute query
            results = session.execute(query).fetchall()
            
            # Process results
            for result in results:
                # Convert result to dictionary
                site_data = {column: getattr(result, column) for column in result._fields}
                site_id = site_data.get('PhosphositeID')
                
                if site_id:
                    # Add to result dictionary
                    result_dict[site_id] = site_data
                    
                    # Cache the result
                    cache_key = f"phosphosite_data_{site_id}"
                    QUERY_CACHE[cache_key] = site_data
            
            # For any site_ids that weren't found in the database, cache as None
            for site_id in uncached_site_ids:
                if site_id not in result_dict:
                    cache_key = f"phosphosite_data_{site_id}"
                    QUERY_CACHE[cache_key] = None
    
    except Exception as e:
        logger.error(f"Error retrieving batch phosphosite data: {e}")
    
    return result_dict

def get_kinase_scores(site_id: str, score_type: str = 'structure') -> Optional[Dict]:
    """
    Get kinase scores for a specific site from the database.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary with kinase scores or None if not found
    """
    # Check cache first
    cache_key = f"{score_type}_kinase_scores_{site_id}"
    if cache_key in QUERY_CACHE:
        return QUERY_CACHE[cache_key]
    
    try:
        with get_db_session() as session:
            # Select the appropriate table based on score_type
            table_name = 'Structure_Kinase_Scores' if score_type.lower() == 'structure' else 'Sequence_Kinase_Scores'
            kinase_table = metadata.tables.get(table_name)
            
            if not kinase_table:
                logger.error(f"{table_name} table not found in database")
                return None
            
            # Create query
            query = select([kinase_table]).where(kinase_table.c.node == site_id)
            
            # Execute query
            result = session.execute(query).fetchone()
            
            if result:
                # Convert result to dictionary
                # This includes 'label' which might be a known_kinase field, and all the kinase columns
                result_dict = {column: getattr(result, column) for column in result._fields}
                
                # Extract the scores into a separate dictionary
                scores = {}
                for column in result._fields:
                    if column not in ['node', 'label']:
                        # Convert any string scores to float
                        try:
                            value = getattr(result, column)
                            if pd.notna(value):  # Check if not NaN
                                scores[column] = float(value)
                            else:
                                scores[column] = 0.0
                        except (ValueError, TypeError):
                            logger.warning(f"Skipping non-numeric value for kinase {column}")
                
                # Structure the return value like the original function
                return_dict = {
                    'known_kinase': result_dict.get('label', 'unlabeled'),
                    'scores': scores
                }
                
                # Cache result
                QUERY_CACHE[cache_key] = return_dict
                return return_dict
            
            return None
    except Exception as e:
        logger.error(f"Error retrieving {score_type} kinase scores for {site_id}: {e}")
        return None

def get_kinase_scores_batch(site_ids: List[str], score_type: str = 'structure') -> Dict[str, Dict]:
    """
    Get kinase scores for multiple sites in a single query.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary mapping site IDs to their kinase score dictionaries
    """
    result_dict = {}
    
    # Return empty dict if no site_ids
    if not site_ids:
        return result_dict
    
    # Check cache first for all site_ids
    uncached_site_ids = []
    for site_id in site_ids:
        cache_key = f"{score_type}_kinase_scores_{site_id}"
        if cache_key in QUERY_CACHE:
            result_dict[site_id] = QUERY_CACHE[cache_key]
        else:
            uncached_site_ids.append(site_id)
    
    # If all site_ids were in cache, return directly
    if not uncached_site_ids:
        return result_dict
    
    try:
        with get_db_session() as session:
            # Select the appropriate table based on score_type
            table_name = 'Structure_Kinase_Scores' if score_type.lower() == 'structure' else 'Sequence_Kinase_Scores'
            kinase_table = metadata.tables.get(table_name)
            
            if not kinase_table:
                logger.error(f"{table_name} table not found in database")
                return result_dict
            
            # Create query for all uncached site_ids
            query = select([kinase_table]).where(kinase_table.c.node.in_(uncached_site_ids))
            
            # Execute query
            results = session.execute(query).fetchall()
            
            # Process results
            for result in results:
                # Get the site_id
                site_id = getattr(result, 'node')
                
                # Extract the scores into a separate dictionary
                scores = {}
                for column in result._fields:
                    if column not in ['node', 'label']:
                        # Convert any string scores to float
                        try:
                            value = getattr(result, column)
                            if pd.notna(value):  # Check if not NaN
                                scores[column] = float(value)
                            else:
                                scores[column] = 0.0
                        except (ValueError, TypeError):
                            logger.warning(f"Skipping non-numeric value for kinase {column}")
                
                # Structure the return value like the original function
                return_dict = {
                    'known_kinase': getattr(result, 'label', 'unlabeled'),
                    'scores': scores
                }
                
                # Add to result dictionary
                result_dict[site_id] = return_dict
                
                # Cache the result
                cache_key = f"{score_type}_kinase_scores_{site_id}"
                QUERY_CACHE[cache_key] = return_dict
            
            # For any site_ids that weren't found in the database, cache as None
            for site_id in uncached_site_ids:
                if site_id not in result_dict:
                    cache_key = f"{score_type}_kinase_scores_{site_id}"
                    QUERY_CACHE[cache_key] = None
    except Exception as e:
        logger.error(f"Error retrieving batch {score_type} kinase scores: {e}")
    
    return result_dict

def find_structural_matches(site_id: str, rmsd_threshold: float = 5.0, top_n: Optional[int] = None) -> List[Dict]:
    """
    Find structural matches for a site from the database.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        rmsd_threshold: Maximum RMSD value to include
        top_n: Maximum number of matches to return
        
    Returns:
        List of match dictionaries sorted by RMSD (closest first)
    """
    # Check cache first
    cache_key = f"structural_matches_{site_id}_{rmsd_threshold}_{top_n}"
    if cache_key in QUERY_CACHE:
        return QUERY_CACHE[cache_key]
    
    matches = []
    
    try:
        with get_db_session() as session:
            # Get the Structural_Similarity_Edges table
            edges_table = metadata.tables.get('Structural_Similarity_Edges')
            
            if not edges_table:
                logger.error("Structural_Similarity_Edges table not found in database")
                return matches
            
            # Create query for matches where the site is either Query or Target
            query = select([edges_table]).where(
                or_(
                    edges_table.c.Query == site_id,
                    edges_table.c.Target == site_id
                )
            ).where(
                edges_table.c.RMSD <= rmsd_threshold
            ).order_by(edges_table.c.RMSD)
            
            # Apply top_n limit if specified
            if top_n is not None:
                query = query.limit(top_n)
            
            # Execute query
            results = session.execute(query).fetchall()
            
            # Process results
            for result in results:
                # Convert result to dictionary
                match_data = {column: getattr(result, column) for column in result._fields}
                
                # Ensure Query is the input site_id by swapping if needed
                if match_data['Query'] != site_id:
                    # Swap Query and Target
                    match_data['Query'], match_data['Target'] = match_data['Target'], match_data['Query']
                
                # Parse Query and Target into components
                query_parts = match_data['Query'].split('_')
                target_parts = match_data['Target'].split('_')
                
                if len(query_parts) > 1 and len(target_parts) > 1:
                    # Create match dictionary in the expected format
                    match = {
                        'query_uniprot': query_parts[0],
                        'query_site': match_data.get('Mapped_Site', ''),  # Use Mapped_Site if available
                        'target_uniprot': target_parts[0],
                        'target_site': match_data.get('Mapped_Site', target_parts[1]),
                        'rmsd': float(match_data['RMSD']),
                        'fold_disco_score': float(match_data.get('FoldDisco_Score', 0)),
                        'mapped_residues': int(match_data.get('Mapped_Residues', 0)),
                        'number_mapped': int(match_data.get('NumberMapped', 0)),
                    }
                    matches.append(match)
            
            # Cache result
            QUERY_CACHE[cache_key] = matches
    except Exception as e:
        logger.error(f"Error finding structural matches for {site_id}: {e}")
    
    return matches

def find_structural_matches_batch(site_ids: List[str], rmsd_threshold: float = 5.0) -> Dict[str, List[Dict]]:
    """
    Find structural matches for multiple sites in a single query.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        rmsd_threshold: Maximum RMSD value to include
        
    Returns:
        Dictionary mapping site IDs to their match lists
    """
    result_dict = {}
    
    # Return empty dict if no site_ids
    if not site_ids:
        return result_dict
    
    # Check cache first for all site_ids
    uncached_site_ids = []
    for site_id in site_ids:
        cache_key = f"structural_matches_{site_id}_{rmsd_threshold}_None"
        if cache_key in QUERY_CACHE:
            result_dict[site_id] = QUERY_CACHE[cache_key]
        else:
            uncached_site_ids.append(site_id)
    
    # If all site_ids were in cache, return directly
    if not uncached_site_ids:
        return result_dict
    
    try:
        with get_db_session() as session:
            # Get the Structural_Similarity_Edges table
            edges_table = metadata.tables.get('Structural_Similarity_Edges')
            
            if not edges_table:
                logger.error("Structural_Similarity_Edges table not found in database")
                return result_dict
            
            # Create query for matches where any site is either Query or Target
            query = select([edges_table]).where(
                or_(
                    edges_table.c.Query.in_(uncached_site_ids),
                    edges_table.c.Target.in_(uncached_site_ids)
                )
            ).where(
                edges_table.c.RMSD <= rmsd_threshold
            ).order_by(edges_table.c.RMSD)
            
            # Execute query
            results = session.execute(query).fetchall()
            
            # First pass: collect matches by site_id
            site_matches = {site_id: [] for site_id in uncached_site_ids}
            
            # Process results
            for result in results:
                # Convert result to dictionary
                match_data = {column: getattr(result, column) for column in result._fields}
                
                query = match_data['Query']
                target = match_data['Target']
                
                # Process for query site
                if query in uncached_site_ids:
                    # Parse Query and Target into components
                    query_parts = query.split('_')
                    target_parts = target.split('_')
                    
                    if len(query_parts) > 1 and len(target_parts) > 1:
                        # Create match dictionary in the expected format
                        match = {
                            'query_uniprot': query_parts[0],
                            'query_site': match_data.get('Mapped_Site', ''),
                            'target_uniprot': target_parts[0],
                            'target_site': match_data.get('Mapped_Site', target_parts[1]),
                            'rmsd': float(match_data['RMSD']),
                            'fold_disco_score': float(match_data.get('FoldDisco_Score', 0)),
                            'mapped_residues': int(match_data.get('Mapped_Residues', 0)),
                            'number_mapped': int(match_data.get('NumberMapped', 0)),
                        }
                        site_matches[query].append(match)
                
                # Process for target site (swapping query and target)
                if target in uncached_site_ids:
                    # Parse Query and Target into components
                    query_parts = query.split('_')
                    target_parts = target.split('_')
                    
                    if len(query_parts) > 1 and len(target_parts) > 1:
                        # Create match dictionary with swapped query and target
                        match = {
                            'query_uniprot': target_parts[0],
                            'query_site': match_data.get('Mapped_Site', ''),
                            'target_uniprot': query_parts[0],
                            'target_site': match_data.get('Mapped_Site', query_parts[1]),
                            'rmsd': float(match_data['RMSD']),
                            'fold_disco_score': float(match_data.get('FoldDisco_Score', 0)),
                            'mapped_residues': int(match_data.get('Mapped_Residues', 0)),
                            'number_mapped': int(match_data.get('NumberMapped', 0)),
                        }
                        site_matches[target].append(match)
            
            # Add to result dictionary and cache
            for site_id, matches in site_matches.items():
                result_dict[site_id] = matches
                
                # Cache the result
                cache_key = f"structural_matches_{site_id}_{rmsd_threshold}_None"
                QUERY_CACHE[cache_key] = matches
    except Exception as e:
        logger.error(f"Error finding batch structural matches: {e}")
    
    return result_dict

def find_sequence_matches(site_id: str, min_similarity: float = 0.4, top_n: Optional[int] = None) -> List[Dict]:
    """
    Find sequence similarity matches for a site from the database.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        min_similarity: Minimum similarity score to include (0-1)
        top_n: Maximum number of results to return
        
    Returns:
        List of dictionaries with match information
    """
    # Check cache first
    cache_key = f"sequence_matches_{site_id}_{min_similarity}_{top_n}"
    if cache_key in QUERY_CACHE:
        return QUERY_CACHE[cache_key]
    
    matches = []
    
    try:
        with get_db_session() as session:
            # Get the Sequence_Similarity_Edges_7 table
            edges_table = metadata.tables.get('Sequence_Similarity_Edges_7')
            
            if not edges_table:
                logger.error("Sequence_Similarity_Edges_7 table not found in database")
                return matches
            
            # Create query for matches where the site is either ID1 or ID2
            query = select([edges_table]).where(
                or_(
                    edges_table.c.ID1 == site_id,
                    edges_table.c.ID2 == site_id
                )
            ).where(
                edges_table.c.Similarity >= min_similarity
            ).order_by(edges_table.c.Similarity.desc())
            
            # Apply top_n limit if specified
            if top_n is not None:
                query = query.limit(top_n)
            
            # Execute query
            results = session.execute(query).fetchall()
            
            # Process results
            for result in results:
                # Convert result to dictionary
                match_data = {column: getattr(result, column) for column in result._fields}
                
                # Ensure ID1 is the input site_id by swapping if needed
                if match_data['ID1'] != site_id:
                    # Swap ID1 and ID2
                    match_data['ID1'], match_data['ID2'] = match_data['ID2'], match_data['ID1']
                
                # Parse target_id into components
                target_id = match_data['ID2']
                target_parts = target_id.split('_')
                
                if len(target_parts) >= 2:
                    target_uniprot = target_parts[0]
                    target_site = target_parts[1]
                    
                    # Create match dictionary
                    match = {
                        'query_id': site_id,
                        'target_id': target_id,
                        'target_uniprot': target_uniprot,
                        'target_site': target_site,
                        'similarity': float(match_data['Similarity'])
                    }
                    
                    # Try to extract site type
                    site_match = re.match(r'([STY])(\d+)', target_site)
                    if site_match:
                        match['site_type'] = site_match.group(1)
                    
                    matches.append(match)
            
            # Cache result
            QUERY_CACHE[cache_key] = matches
    except Exception as e:
        logger.error(f"Error finding sequence matches for {site_id}: {e}")
    
    return matches

def find_sequence_matches_batch(site_ids: List[str], min_similarity: float = 0.4) -> Dict[str, List[Dict]]:
    """
    Find sequence similarity matches for multiple sites in a single query.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        min_similarity: Minimum similarity score to include (0-1)
        
    Returns:
        Dictionary mapping site IDs to their match lists
    """
    result_dict = {}
    
    # Return empty dict if no site_ids
    if not site_ids:
        return result_dict
    
    # Check cache first for all site_ids
    uncached_site_ids = []
    for site_id in site_ids:
        cache_key = f"sequence_matches_{site_id}_{min_similarity}_None"
        if cache_key in QUERY_CACHE:
            result_dict[site_id] = QUERY_CACHE[cache_key]
        else:
            uncached_site_ids.append(site_id)
    
    # If all site_ids were in cache, return directly
    if not uncached_site_ids:
        return result_dict
    
    try:
        with get_db_session() as session:
            # Get the Sequence_Similarity_Edges_7 table
            edges_table = metadata.tables.get('Sequence_Similarity_Edges_7')
            
            if not edges_table:
                logger.error("Sequence_Similarity_Edges_7 table not found in database")
                return result_dict
            
            # Create query for matches where any site is either ID1 or ID2
            query = select([edges_table]).where(
                or_(
                    edges_table.c.ID1.in_(uncached_site_ids),
                    edges_table.c.ID2.in_(uncached_site_ids)
                )
            ).where(
                edges_table.c.Similarity >= min_similarity
            ).order_by(edges_table.c.Similarity.desc())
            
            # Execute query
            results = session.execute(query).fetchall()
            
            # First pass: collect matches by site_id
            site_matches = {site_id: [] for site_id in uncached_site_ids}
            
            # Process results
            for result in results:
                # Convert result to dictionary
                match_data = {column: getattr(result, column) for column in result._fields}
                
                id1 = match_data['ID1']
                id2 = match_data['ID2']
                similarity = float(match_data['Similarity'])
                
                # Process for ID1 site
                if id1 in uncached_site_ids:
                    # Parse ID2 into components
                    target_parts = id2.split('_')
                    
                    if len(target_parts) >= 2:
                        target_uniprot = target_parts[0]
                        target_site = target_parts[1]
                        
                        # Create match dictionary
                        match = {
                            'query_id': id1,
                            'target_id': id2,
                            'target_uniprot': target_uniprot,
                            'target_site': target_site,
                            'similarity': similarity
                        }
                        
                        # Try to extract site type
                        site_match = re.match(r'([STY])(\d+)', target_site)
                        if site_match:
                            match['site_type'] = site_match.group(1)
                        
                        site_matches[id1].append(match)
                
                # Process for ID2 site (swapping ID1 and ID2)
                if id2 in uncached_site_ids:
                    # Parse ID1 into components
                    target_parts = id1.split('_')
                    
                    if len(target_parts) >= 2:
                        target_uniprot = target_parts[0]
                        target_site = target_parts[1]
                        
                        # Create match dictionary with swapped ID1 and ID2
                        match = {
                            'query_id': id2,
                            'target_id': id1,
                            'target_uniprot': target_uniprot,
                            'target_site': target_site,
                            'similarity': similarity
                        }
                        
                        # Try to extract site type
                        site_match = re.match(r'([STY])(\d+)', target_site)
                        if site_match:
                            match['site_type'] = site_match.group(1)
                        
                        site_matches[id2].append(match)
            
            # Add to result dictionary and cache
            for site_id, matches in site_matches.items():
                result_dict[site_id] = matches
                
                # Cache the result
                cache_key = f"sequence_matches_{site_id}_{min_similarity}_None"
                QUERY_CACHE[cache_key] = matches
    except Exception as e:
        logger.error(f"Error finding batch sequence matches: {e}")
    
    return result_dict

def clear_cache():
    """Clear the query cache."""
    global QUERY_CACHE
    QUERY_CACHE = {}
    logger.info("Query cache cleared")

def get_db_stats():
    """
    Get database statistics.
    
    Returns:
        Dictionary with database statistics
    """
    stats = {
        'table_counts': {},
        'cache_size': len(QUERY_CACHE)
    }
    
    try:
        with get_db_session() as session:
            # Get counts for each table
            for table_name in metadata.tables:
                table = metadata.tables[table_name]
                query = select([sqlalchemy.func.count()]).select_from(table)
                count = session.execute(query).scalar()
                stats['table_counts'][table_name] = count
    except Exception as e:
        logger.error(f"Error getting database stats: {e}")
    
    return stats