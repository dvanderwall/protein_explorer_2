
"""
Database package for KinoPlex application.
Re-exports all functions from db.py for consistent imports.
"""

# Import and re-export all functions from db.py
from protein_explorer.db.db import (
    init_db, get_phosphosite_data, get_phosphosites_batch,
    find_structural_matches, find_structural_matches_batch,
    find_sequence_matches, find_sequence_matches_batch,
    get_kinase_scores, get_kinase_scores_batch
)

# Make these functions available when importing from protein_explorer.db
__all__ = [
    'init_db', 'get_phosphosite_data', 'get_phosphosites_batch',
    'find_structural_matches', 'find_structural_matches_batch',
    'find_sequence_matches', 'find_sequence_matches_batch',
    'get_kinase_scores', 'get_kinase_scores_batch'
]