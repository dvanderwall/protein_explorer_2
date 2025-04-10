"""
Database access module for KinoPlex.

This module provides functions to connect to the database and execute queries against
the phosphosite tables. It includes optimized batch query functionality that leverages
database indexes for efficient data retrieval.

Indexed tables and columns:
- Phosphosite_Supplementary_Data: PhosphositeID, uniprot_id, GENE
- Structural_Similarity_Edges: Query+RMSD, Target
- Sequence_Similarity_Edges_7: ID1+Similarity, ID2+Similarity
- Structure_Kinase_Scores: node
- Sequence_Kinase_Scores: node
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
        ],
        "indexes": {
            "PhosphositeID": "idx_phosphositeID",
            "uniprot_id": "idx_uniprot_id",
            "GENE": "idx_gene"
        }
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
        ],
        "indexes": {
            "node": "idx_node_sequence"
        }
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
        ],
        "indexes": {
            "node": "idx_node_structure"
        }
    },
    "Structural_Similarity_Edges": {
        "columns": [
            "Query", "Target", "GlobalScore", "QueryNeighborhood", 
            "TargetNeighborhood", "QueryNeighborhoodSize", "TargetNeighborhoodSize", 
            "MappedQueryResidues", "MappedTargetResidues", "MappingRatio", "RMSD", 
            "QueryPercentHydrophobicity", "TargetPercentHydrophobicity", "QueryPercentHydrophilicity",  
            "TargetPercentHydrophilicity", "QueryPercentAcidity", "TargetPercentAcidity", 
            "QueryPercentBasicity", "TargetPercentBasicity", "QueryPercentAromatic", 
            "TargetPercentAromatic", "QueryPocket_pI", "TargetPocket_pI", "QueryPocket_pKa", 
            "TargetPocket_pKa"
        ],
        "indexes": {
            "Query": "idx_query_rmsd",
            "RMSD": "idx_query_rmsd",
            "Target": "idx_target"
        }
    },
    "Sequence_Similarity_Edges_Full": {
        "columns": [
            "ID1", "ID2", "Similarity"
        ],
        "indexes": {
            "ID1": "idx_id1_sim",
            "ID2": "idx_id2_sim",
            "Similarity": ["idx_id1_sim", "idx_id2_sim"]
        }
    },
    "Sequence_Similarity_Edges_4": {
        "columns": [
            "ID1", "ID2", "Similarity"
        ],
        "indexes": {
            "ID1": "idx_id1_sim",
            "ID2": "idx_id2_sim",
            "Similarity": ["idx_id1_sim", "idx_id2_sim"]
        }
    },
    "Sequence_Similarity_Edges_7": {
        "columns": [
            "ID1", "ID2", "Similarity"
        ],
        "indexes": {
            "ID1": "idx_id1_sim",
            "ID2": "idx_id2_sim",
            "Similarity": ["idx_id1_sim", "idx_id2_sim"]
        }
    }, 
    "STY_Structural_Annotations": {
        "id_column": "PhosphositeID",
        "columns": [
            "UniProtID", "ResidueNumber", "Site", "ResidueType", "Motif", 
            "pLDDT", "NeighborCount", "ChainID", "SecondaryStructure", 
            "SecondaryStructureDesc", "HydroxylExposure", "BackboneContacts", 
            "HydrogenBonds", "MotifPLDDT", "SeqDistToHelix", "SpatialDistToHelix", 
            "SeqDistToStrand", "SpatialDistToStrand", "SeqDistToTurn", 
            "SpatialDistToTurn", "hydrophobic_count", "polar_count", 
            "positive_count", "negative_count", "net_charge", "small_count", 
            "medium_count", "large_count", "hydrophobic_ratio", "charged_ratio", 
            "HSE_CA_U", "HSE_CA_D", "HSE_CA_RATIO", "HSE_CB_U", "HSE_CB_D", 
            "HSE_CB_RATIO", "PhosphositeID"
        ],
        "indexes": {
            "PhosphositeID": "PRIMARY",
            "UniProtID": "idx_uniprot_residue",
            "ResidueNumber": "idx_uniprot_residue",
            "ResidueType": "idx_residue_type",
            "pLDDT": "idx_plddt",
            "SecondaryStructure": "idx_secondary_structure"
        }
    }, 
    "Cantley_Kinome_Scores_All_STs": {
        "id_column": "SiteID",
        "columns": [
            "SiteID", "Motif", "AAK1", "ACVR2A", "ACVR2B", "AKT1", "AKT2", "AKT3", "ALK2", "ALK4", 
            "ALPHAK3", "AMPKA1", "AMPKA2", "ANKRD3", "ASK1", "ATM", "ATR", "AURA", "AURB", "AURC",
            "BCKDK", "BIKE", "BMPR1A", "BMPR1B", "BMPR2", "BRAF", "BRSK1", "BRSK2", "BUB1", "CAMK1A",
            "CAMK1B", "CAMK1D", "CAMK1G", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G", "CAMK4", "CAMKK1", "CAMKK2",
            "CAMLCK", "CDC7", "CDK1", "CDK10", "CDK12", "CDK13", "CDK14", "CDK16", "CDK17", "CDK18",
            "CDK19", "CDK2", "CDK3", "CDK4", "CDK5", "CDK6", "CDK7", "CDK8", "CDK9", "CDKL1",
            "CDKL5", "CHAK1", "CHAK2", "CHK1", "CHK2", "CK1A", "CK1A2", "CK1D", "CK1E", "CK1G1",
            "CK1G2", "CK1G3", "CK2A1", "CK2A2", "CLK1", "CLK2", "CLK3", "CLK4", "COT", "CRIK",
            "DAPK1", "DAPK2", "DAPK3", "DCAMKL1", "DCAMKL2", "DLK", "DMPK1", "DNAPK", "DRAK1", "DSTYK",
            "DYRK1A", "DYRK1B", "DYRK2", "DYRK3", "DYRK4", "EEF2K", "ERK1", "ERK2", "ERK5", "ERK7",
            "FAM20C", "GAK", "GCK", "GCN2", "GRK1", "GRK2", "GRK3", "GRK4", "GRK5", "GRK6",
            "GRK7", "GSK3A", "GSK3B", "HASPIN", "HGK", "HIPK1", "HIPK2", "HIPK3", "HIPK4", "HPK1",
            "HRI", "HUNK", "ICK", "IKKA", "IKKB", "IKKE", "IRAK1", "IRAK4", "IRE1", "IRE2",
            "JNK1", "JNK2", "JNK3", "KHS1", "KHS2", "KIS", "LATS1", "LATS2", "LKB1", "LOK",
            "LRRK2", "MAK", "MAP3K15", "MAPKAPK2", "MAPKAPK3", "MAPKAPK5", "MARK1", "MARK2", "MARK3", "MARK4",
            "MASTL", "MEK1", "MEK2", "MEK5", "MEKK1", "MEKK2", "MEKK3", "MEKK6", "MELK", "MINK",
            "MLK1", "MLK2", "MLK3", "MLK4", "MNK1", "MNK2", "MOK", "MOS", "MPSK1", "MRCKA",
            "MRCKB", "MSK1", "MSK2", "MST1", "MST2", "MST3", "MST4", "MTOR", "MYLK4", "MYO3A",
            "MYO3B", "NDR1", "NDR2", "NEK1", "NEK11", "NEK2", "NEK3", "NEK4", "NEK5", "NEK6",
            "NEK7", "NEK8", "NEK9", "NIK", "NIM1", "NLK", "NUAK1", "NUAK2", "OSR1", "P38A",
            "P38B", "P38D", "P38G", "P70S6K", "P70S6KB", "P90RSK", "PAK1", "PAK2", "PAK3", "PAK4",
            "PAK5", "PAK6", "PASK", "PBK", "PDHK1", "PDHK4", "PDK1", "PERK", "PHKG1", "PHKG2",
            "PIM1", "PIM2", "PIM3", "PINK1", "PKACA", "PKACB", "PKACG", "PKCA", "PKCB", "PKCD",
            "PKCE", "PKCG", "PKCH", "PKCI", "PKCT", "PKCZ", "PKG1", "PKG2", "PKN1", "PKN2",
            "PKN3", "PKR", "PLK1", "PLK2", "PLK3", "PLK4", "PRKD1", "PRKD2", "PRKD3", "PRKX",
            "PRP4", "PRPK", "QIK", "QSK", "RAF1", "RIPK1", "RIPK2", "RIPK3", "ROCK1", "ROCK2",
            "RSK2", "RSK3", "RSK4", "SBK", "SGK1", "SGK3", "SIK", "SKMLCK", "SLK", "SMG1",
            "SMMLCK", "SNRK", "SRPK1", "SRPK2", "SRPK3", "SSTK", "STK33", "STLK3", "TAK1", "TAO1",
            "TAO2", "TAO3", "TBK1", "TGFBR1", "TGFBR2", "TLK1", "TLK2", "TNIK", "TSSK1", "TSSK2",
            "TTBK1", "TTBK2", "TTK", "ULK1", "ULK2", "VRK1", "VRK2", "WNK1", "WNK3", "WNK4",
            "YANK2", "YANK3", "YSK1", "YSK4", "ZAK"
        ],
        "indexes": {
            "SiteID": "PRIMARY"
        }
    },

    "Cantley_Kinome_Scores_All_Ys": {
        "id_column": "SiteID",
        "columns": [
            "SiteID", "Motif", "ABL", "ACK", "ALK", "ARG", "AXL", "BLK", "BMPR2_TYR", "BRK", 
            "BTK", "CSFR", "CSK", "CTK", "DDR1", "DDR2", "EGFR", "EPHA1", "EPHA2", "EPHA3", 
            "EPHA4", "EPHA5", "EPHA6", "EPHA7", "EPHA8", "EPHB1", "EPHB2", "EPHB3", "EPHB4", 
            "ETK", "FAK", "FER", "FES", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FGR", "FLT3", 
            "FRK", "FYN", "HCK", "HER2", "HER4", "IGF1R", "INSR", "IRR", "ITK", "JAK1", 
            "JAK2", "JAK3", "KIT", "LCK", "LIMK1_TYR", "LIMK2_TYR", "LTK", "LYN", "MER", 
            "MET", "MKK4_TYR", "MKK6_TYR", "MKK7_TYR", "MST1R", "MUSK", "MYT1_TYR", "NEK10_TYR", 
            "PDGFRA", "PDGFRB", "PDHK1_TYR", "PDHK3_TYR", "PDHK4_TYR", "PINK1_TYR", "PYK2", 
            "RET", "ROS", "SRC", "SRMS", "SYK", "TEC", "TESK1_TYR", "TIE2", "TNK1", 
            "TNNI3K_TYR", "TRKA", "TRKB", "TRKC", "TXK", "TYK2", "TYRO3", "VEGFR1", "VEGFR2", 
            "VEGFR3", "WEE1_TYR", "YES", "ZAP70"
        ],
        "indexes": {
            "SiteID": "PRIMARY"
        }
    }
}

# Batch sizes optimized for indexed queries
# These sizes are tuned based on the indexes and table sizes
BATCH_SIZES = {
    "phosphosites": 1000,        # Increased due to idx_phosphositeID index
    "structural_matches": 300,    # Optimized for idx_query_rmsd
    "sequence_matches": 200,      # Optimized for idx_id1_sim and idx_id2_sim
    "kinase_scores": 800         # Optimized for idx_node_structure/sequence
}

# Query templates for MySQL - optimized for the new indexes
QUERY_TEMPLATES = {
    "get_phosphosite_data": """
        SELECT * FROM `Phosphosite_Supplementary_Data`
        WHERE `PhosphositeID` = :site_id
        /* Uses index: idx_phosphositeID */
    """,
    
    "get_phosphosites_batch": """
        SELECT * FROM `Phosphosite_Supplementary_Data`
        WHERE `PhosphositeID` IN :site_ids
        /* Uses index: idx_phosphositeID */
    """,
    
    "find_structural_matches": """
        SELECT * FROM `Structural_Similarity_Edges`
        WHERE `Query` = :site_id AND `RMSD` <= :rmsd_threshold
        ORDER BY `RMSD` ASC
        /* Uses index: idx_query_rmsd (composite) */
    """,
    
    "find_structural_matches_batch": """
        SELECT * FROM `Structural_Similarity_Edges`
        WHERE `Query` IN :site_ids AND `RMSD` <= :rmsd_threshold
        ORDER BY `RMSD` ASC
        /* Uses index: idx_query_rmsd (composite) */
    """,

    "find_structural_matches_batch_by_target": """
        SELECT * FROM `Structural_Similarity_Edges`
        WHERE `Target` IN :site_ids AND `RMSD` <= :rmsd_threshold
        ORDER BY `RMSD` ASC
        /* Uses index: idx_query_rmsd (composite) */
    """,
    
    "find_sequence_matches": """
        SELECT * FROM `Sequence_Similarity_Edges_Full`
        WHERE (`ID1` = :site_id OR `ID2` = :site_id) AND `Similarity` >= :min_similarity
        ORDER BY `Similarity` DESC
        /* Uses indexes: idx_id1_sim and idx_id2_sim (composite) */
    """,
    
    "find_sequence_matches_batch": """
        SELECT * FROM `Sequence_Similarity_Edges_Full`
        WHERE (`ID1` IN :site_ids OR `ID2` IN :site_ids) AND `Similarity` >= :min_similarity
        ORDER BY `Similarity` DESC
        /* Uses indexes: idx_id1_sim and idx_id2_sim (composite) */
    """,
    
    "get_structure_kinase_scores": """
        SELECT * FROM `Structure_Kinase_Scores`
        WHERE `node` = :site_id
        /* Uses index: idx_node_structure */
    """,
    
    "get_sequence_kinase_scores": """
        SELECT * FROM `Sequence_Kinase_Scores`
        WHERE `node` = :site_id
        /* Uses index: idx_node_sequence */
    """,
    
    "get_structure_kinase_scores_batch": """
        SELECT * FROM `Structure_Kinase_Scores`
        WHERE `node` IN :site_ids
        /* Uses index: idx_node_structure */
    """,
    
    "get_sequence_kinase_scores_batch": """
        SELECT * FROM `Sequence_Kinase_Scores`
        WHERE `node` IN :site_ids
        /* Uses index: idx_node_sequence */
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
        SELECT COUNT(*) as count FROM `Sequence_Similarity_Edges_Full`
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
        /* Uses index: idx_phosphositeID */
    """,
    
    "find_sites_by_protein": """
        SELECT * FROM `Phosphosite_Supplementary_Data`
        WHERE `uniprot_id` = :uniprot_id
        /* Uses index: idx_uniprot_id */
    """,
    
    "find_sites_by_gene": """
        SELECT * FROM `Phosphosite_Supplementary_Data`
        WHERE `GENE` = :gene_symbol
        /* Uses index: idx_gene */
    """,
    
    "get_top_similar_sites": """
        SELECT `ID1`, `ID2`, `Similarity` 
        FROM `Sequence_Similarity_Edges_Full`
        WHERE (`ID1` = :site_id OR `ID2` = :site_id)
        AND `Similarity` >= :min_similarity
        ORDER BY `Similarity` DESC
        LIMIT :limit
        /* Uses indexes: idx_id1_sim and idx_id2_sim (composite) */
    """,
    
    "get_top_structural_matches": """
        SELECT * 
        FROM `Structural_Similarity_Edges`
        WHERE `Query` = :site_id
        ORDER BY `RMSD` ASC
        LIMIT :limit
        /* Uses index: idx_query_rmsd */
    """,
    
    "get_top_kinases": """
        SELECT :kinase_columns
        FROM `:score_table`
        WHERE `node` = :site_id
        /* Uses indexes: idx_node_structure or idx_node_sequence */
    """,
    
    # New optimized batch queries leveraging indexes
    "find_sites_by_protein_batch": """
        SELECT * FROM `Phosphosite_Supplementary_Data`
        WHERE `uniprot_id` IN :uniprot_ids
        /* Uses index: idx_uniprot_id */
    """,
    
    "find_sites_by_gene_batch": """
        SELECT * FROM `Phosphosite_Supplementary_Data`
        WHERE `GENE` IN :gene_symbols
        /* Uses index: idx_gene */
    """,
    
    "find_structural_matches_by_target_batch": """
        SELECT * FROM `Structural_Similarity_Edges`
        WHERE `Target` IN :target_ids AND `RMSD` <= :rmsd_threshold
        ORDER BY `RMSD` ASC
        /* Uses index: idx_target */
    """, 
    "get_structural_annotations": """
        SELECT * FROM `STY_Structural_Annotations`
        WHERE `Site` = :site_id
        /* Uses primary key */
    """,

    "get_structural_annotations_batch": """
        SELECT * FROM `STY_Structural_Annotations`
        WHERE `Site` IN :site_ids
        /* Uses primary key */
    """,

    "find_structural_annotations_by_uniprot": """
        SELECT * FROM `STY_Structural_Annotations`
        WHERE `UniProtID` = :uniprot_id
        /* Uses index: idx_uniprot_residue */
    """,

    "find_structural_annotations_by_residue_type": """
        SELECT * FROM `STY_Structural_Annotations`
        WHERE `ResidueType` = :residue_type
        /* Uses index: idx_residue_type */
    """,

    "find_high_quality_sites": """
        SELECT * FROM `STY_Structural_Annotations`
        WHERE `pLDDT` >= :min_plddt
        ORDER BY `pLDDT` DESC
        /* Uses index: idx_plddt */
    """,

    "find_sites_by_secondary_structure": """
        SELECT * FROM `STY_Structural_Annotations`
        WHERE `SecondaryStructure` = :structure_type
        /* Uses index: idx_secondary_structure */
    """,

    "find_exposed_sites": """
        SELECT * FROM `STY_Structural_Annotations`
        WHERE `HydroxylExposure` >= :min_exposure
        ORDER BY `HydroxylExposure` DESC
    """,

    "count_structural_annotations": """
        SELECT COUNT(*) as count FROM `STY_Structural_Annotations`
    """, 
    "get_cantley_st_kinase_scores": """
        SELECT * FROM `Cantley_Kinome_Scores_All_STs`
        WHERE `SiteID` = :site_id
        /* Uses primary key */
    """,

    "get_cantley_st_kinase_scores_batch": """
        SELECT * FROM `Cantley_Kinome_Scores_All_STs`
        WHERE `SiteID` IN :site_ids
        /* Uses primary key */
    """,

    "get_cantley_y_kinase_scores": """
        SELECT * FROM `Cantley_Kinome_Scores_All_Ys`
        WHERE `SiteID` = :site_id
        /* Uses primary key */
    """,

    "get_cantley_y_kinase_scores_batch": """
        SELECT * FROM `Cantley_Kinome_Scores_All_Ys`
        WHERE `SiteID` IN :site_ids
        /* Uses primary key */
    """,

    "get_cantley_st_top_kinases": """
        SELECT `SiteID`, :kinase_names, `Motif`
        FROM `Cantley_Kinome_Scores_All_STs`
        WHERE `SiteID` = :site_id
        /* Uses primary key */
    """,

    "get_cantley_y_top_kinases": """
        SELECT `SiteID`, :kinase_names, `Motif`
        FROM `Cantley_Kinome_Scores_All_Ys`
        WHERE `SiteID` = :site_id
        /* Uses primary key */
    """,

    "get_sites_by_cantley_st_kinase": """
        SELECT `SiteID`, `Motif`, `:kinase_name`
        FROM `Cantley_Kinome_Scores_All_STs`
        WHERE `:kinase_name` >= :score_threshold
        ORDER BY `:kinase_name` DESC
        LIMIT :limit
        /* Uses full table scan with filter */
    """,

    "get_sites_by_cantley_y_kinase": """
        SELECT `SiteID`, `Motif`, `:kinase_name`
        FROM `Cantley_Kinome_Scores_All_Ys`
        WHERE `:kinase_name` >= :score_threshold
        ORDER BY `:kinase_name` DESC
        LIMIT :limit
        /* Uses full table scan with filter */
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

# Performance metrics
PERFORMANCE_METRICS = {
    "query_times": {},
    "cache_hits": 0,
    "cache_misses": 0,
    "db_errors": 0
}

def init_db():
    """Initialize database connection with optimized connection pooling settings."""
    global engine, Session, metadata
    try:
        engine = create_engine(
            DB_URL,
            poolclass=QueuePool,
            pool_size=15,  # Increased pool size for better concurrent query performance
            max_overflow=30,  # Increased overflow connections
            pool_timeout=30,
            pool_recycle=1800,  # Ensure connections are recycled every 30 minutes
            pool_pre_ping=True  # Add connection health check
        )
        # Create session maker
        Session = sessionmaker(bind=engine)
        
        # Test connection
        with engine.connect() as conn:
            conn.execute(text("SELECT 1"))
            
        # Set up metadata for reflection
        metadata.reflect(engine)
        
        logger.info("Database connection initialized successfully")

        print_structural_edges_sample(100)

        return True
    except Exception as e:
        logger.error(f"Error initializing database: {e}")
        PERFORMANCE_METRICS["db_errors"] += 1
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
    Handles MySQL's parameter binding quirks and tracks query performance.
    
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
        PERFORMANCE_METRICS["db_errors"] += 1
        return pd.DataFrame()
    
    # Strip comment lines with index usage info - for logging only
    query_key = re.sub(r'/\*.*?\*/', '', query).strip()
    
    # Track query execution time
    start_time = time.time()
    
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
            
            # Record performance metrics
            query_time = time.time() - start_time
            if query_key in PERFORMANCE_METRICS["query_times"]:
                PERFORMANCE_METRICS["query_times"][query_key] = (
                    PERFORMANCE_METRICS["query_times"][query_key] * 0.9 + query_time * 0.1  # Weighted average
                )
            else:
                PERFORMANCE_METRICS["query_times"][query_key] = query_time
                
            if query_time > 1.0:  # Log slow queries (> 1 second)
                logger.warning(f"Slow query detected ({query_time:.2f}s): {query_key[:100]}...")
                
            return df
    except SQLAlchemyError as e:
        PERFORMANCE_METRICS["db_errors"] += 1
        logger.error(f"Database query error: {e}")
        logger.error(f"Query was: {query}")
        logger.error(f"Parameters were: {params}")
        return pd.DataFrame()
    except Exception as e:
        PERFORMANCE_METRICS["db_errors"] += 1
        logger.error(f"Error executing query: {e}")
        return pd.DataFrame()
    
def execute_batch_query(query_template: str, id_list: List[str], batch_size: int = None, 
                      id_param_name: str = "site_ids", extra_params: dict = None) -> pd.DataFrame:
    """
    Execute a batch query by breaking up a large list into smaller batches.
    Uses appropriate batch sizes based on the table and indexes.
    
    Args:
        query_template: SQL query template with :id_param_name placeholder
        id_list: List of IDs to process in batches
        batch_size: Maximum number of IDs per batch (if None, uses optimized default from BATCH_SIZES)
        id_param_name: Name of the ID parameter in the query
        extra_params: Additional query parameters
        
    Returns:
        Combined DataFrame with results from all batches
    """
    if not id_list:
        return pd.DataFrame()
    
    # Determine appropriate batch size based on query type
    if batch_size is None:
        if "Phosphosite_Supplementary_Data" in query_template:
            batch_size = BATCH_SIZES["phosphosites"]
        elif "Structural_Similarity_Edges" in query_template:
            batch_size = BATCH_SIZES["structural_matches"]
        elif "Sequence_Similarity_Edges_Full" in query_template:
            batch_size = BATCH_SIZES["sequence_matches"]
        elif "Kinase_Scores" in query_template:
            batch_size = BATCH_SIZES["kinase_scores"]
        else:
            batch_size = 500  # Default batch size
    
    # Initialize an empty result DataFrame
    combined_results = pd.DataFrame()
    
    # Calculate total number of batches for progress tracking
    total_batches = (len(id_list) + batch_size - 1) // batch_size
    
    # Process in batches
    for i in range(0, len(id_list), batch_size):
        batch = id_list[i:i + batch_size]
        batch_num = i // batch_size + 1
        logger.info(f"Processing batch {batch_num}/{total_batches} with {len(batch)} items")
        
        # Prepare parameters
        params = {id_param_name: tuple(batch)}
        if extra_params:
            params.update(extra_params)
        
        # Execute query for this batch with timing
        start_time = time.time()
        batch_results = execute_query(query_template, params)
        query_time = time.time() - start_time
        
        # Log batch performance
        if query_time > 1.0:  # Only log slow batches
            logger.info(f"Batch {batch_num}/{total_batches} completed in {query_time:.2f}s")
        
        # Append results
        if not batch_results.empty:
            if combined_results.empty:
                combined_results = batch_results
            else:
                combined_results = pd.concat([combined_results, batch_results], ignore_index=True)
    
    logger.info(f"Completed processing {total_batches} batches with {len(id_list)} total items")
    
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
    Uses idx_phosphositeID index for optimized lookup.
    Ensures that motif sequences are returned in uppercase.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary with phosphosite data or None if not found
    """
    # Check cache first
    cache_key = f"phosphosite_data_{site_id}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        PERFORMANCE_METRICS["cache_hits"] += 1
        result = QUERY_CACHE[cache_key]
        # Ensure motif is uppercase if present
        if result and ('SITE_+/-7_AA' in result) and result['SITE_+/-7_AA']:
            result['SITE_+/-7_AA'] = result['SITE_+/-7_AA'].upper()
        if result and ('motif' in result) and result['motif']:
            result['motif'] = result['motif'].upper()
        return result
    
    PERFORMANCE_METRICS["cache_misses"] += 1
    
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
        
        # Ensure motif sequences are uppercase
        if 'SITE_+/-7_AA' in result and result['SITE_+/-7_AA']:
            result['SITE_+/-7_AA'] = result['SITE_+/-7_AA'].upper()
        if 'motif' in result and result['motif']:
            result['motif'] = result['motif'].upper()
        
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
    Uses idx_phosphositeID index for optimized lookup.
    Ensures that motif sequences are returned in uppercase.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary mapping site IDs to phosphosite data dictionaries
    """
    if not site_ids:
        return {}
    
    # Check cache first for all sites
    cached_results = {}
    missing_sites = []
    
    for site_id in site_ids:
        cache_key = f"phosphosite_data_{site_id}"
        if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
            PERFORMANCE_METRICS["cache_hits"] += 1
            cached_value = QUERY_CACHE[cache_key]
            if cached_value is not None:  # Don't include None values
                # Ensure motif sequences are uppercase
                if 'SITE_+/-7_AA' in cached_value and cached_value['SITE_+/-7_AA']:
                    cached_value['SITE_+/-7_AA'] = cached_value['SITE_+/-7_AA'].upper()
                if 'motif' in cached_value and cached_value['motif']:
                    cached_value['motif'] = cached_value['motif'].upper()
                cached_results[site_id] = cached_value
        else:
            PERFORMANCE_METRICS["cache_misses"] += 1
            missing_sites.append(site_id)
    
    # If all results were in cache, return immediately
    if not missing_sites:
        return cached_results
        
    try:
        # Use the batch helper function for missing sites
        query = QUERY_TEMPLATES["get_phosphosites_batch"]
        # Use optimized batch size from BATCH_SIZES
        df = execute_batch_query(
            query, 
            missing_sites, 
            batch_size=BATCH_SIZES["phosphosites"],
            id_param_name="site_ids"
        )
        
        if not df.empty:
            # Convert to dictionary of dictionaries
            for _, row in df.iterrows():
                row_dict = row.to_dict()
                if "PhosphositeID" in row_dict:
                    site_id = row_dict["PhosphositeID"]
                    
                    # Ensure motif sequences are uppercase
                    if 'SITE_+/-7_AA' in row_dict and row_dict['SITE_+/-7_AA']:
                        row_dict['SITE_+/-7_AA'] = row_dict['SITE_+/-7_AA'].upper()
                    if 'motif' in row_dict and row_dict['motif']:
                        row_dict['motif'] = row_dict['motif'].upper()
                    
                    # Add to results
                    cached_results[site_id] = row_dict
                    # Also update cache
                    cache_key = f"phosphosite_data_{site_id}"
                    QUERY_CACHE[cache_key] = row_dict
                    CACHE_TIMESTAMPS[cache_key] = time.time()
        
        # Cache negative results for missing sites
        for site_id in missing_sites:
            if site_id not in cached_results:
                cache_key = f"phosphosite_data_{site_id}"
                QUERY_CACHE[cache_key] = None
                CACHE_TIMESTAMPS[cache_key] = time.time()
                
        logger.info(f"Retrieved data for {len(cached_results)} out of {len(site_ids)} phosphosites")
        return cached_results
    except Exception as e:
        logger.error(f"Error getting batch phosphosite data: {e}")
        return cached_results  # Return whatever was found in cache even if the rest failed

def find_structural_matches(site_id: str, rmsd_threshold: float = 5.0) -> List[Dict]:
    """
    Find structural matches for a site.
    Uses composite idx_query_rmsd index for optimized lookup.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        rmsd_threshold: Maximum RMSD value for matches
        
    Returns:
        List of dictionaries with match information
    """
    # Check cache first
    cache_key = f"structural_matches_{site_id}_{rmsd_threshold}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        PERFORMANCE_METRICS["cache_hits"] += 1
        return QUERY_CACHE[cache_key]
    
    PERFORMANCE_METRICS["cache_misses"] += 1
    
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
                        "map_ratio": float(row_dict.get("MappingRatio", 0))
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
    Uses composite idx_query_rmsd index for optimized lookup.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        rmsd_threshold: Maximum RMSD value for matches
        
    Returns:
        Dictionary mapping site IDs to lists of match dictionaries
    """
    if not site_ids:
        return {}
    
    # Check cache first for all sites
    cached_results = {}
    missing_sites = []
    
    for site_id in site_ids:
        cache_key = f"structural_matches_{site_id}_{rmsd_threshold}"
        if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
            PERFORMANCE_METRICS["cache_hits"] += 1
            cached_results[site_id] = QUERY_CACHE[cache_key]
        else:
            PERFORMANCE_METRICS["cache_misses"] += 1
            missing_sites.append(site_id)
    
    # If all results were in cache, return immediately
    print("MISSING SITES")
    print(missing_sites) 
    if not missing_sites:
        return cached_results
    print("MISSING SITES")
    print(missing_sites)    
    try:
        # Use the batch helper function for missing sites
        query = QUERY_TEMPLATES["find_structural_matches_batch"]
        df = execute_batch_query(
            query, 
            missing_sites, 
            batch_size=BATCH_SIZES["structural_matches"],  # Use optimized batch size
            id_param_name="site_ids", 
            extra_params={"rmsd_threshold": rmsd_threshold}
        )
        query_2 = QUERY_TEMPLATES["find_structural_matches_batch_by_target"]
        df_2 = execute_batch_query(
            query_2,
            missing_sites,  
            batch_size=BATCH_SIZES["structural_matches"],  # Use optimized batch size
            id_param_name="site_ids",
            extra_params={"rmsd_threshold": rmsd_threshold}
        )
        # Identify the correct column names in df_2
        query_col = "Query"
        target_col = "Target"

        # Create a copy of df_2 to avoid modifying the original
        df_2_swapped = df_2.copy()

        # Swap the Query and Target columns in df_2
        df_2_swapped[query_col] = df_2[target_col]
        df_2_swapped[target_col] = df_2[query_col]

        # Merge df and df_2_swapped
        merged_df = pd.concat([df, df_2_swapped], ignore_index=True)

        # Remove duplicate rows based on Query and Target columns
        merged_df = merged_df.drop_duplicates(subset=[query_col, target_col], keep="first")

        # Reset index
        df = merged_df.reset_index(drop=True)

        print('df1')
        print(df)
        print('df2')
        print(df_2_swapped)
        print('Here is Pre-Structural Data Structure')
        #print(merged_df)
        if not df.empty:
            # Organize matches by query site
            for _, row in df.iterrows():
                row_dict = row.to_dict()
                query_id = row_dict.get("Query")
                
                if query_id not in cached_results:
                    cached_results[query_id] = []
                    
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
                        "map_ratio": float(row_dict.get("MappingRatio", 0))
                    }
                    
                    cached_results[query_id].append(match_dict)
        
        # Update cache with new results
        for site_id in missing_sites:
            cache_key = f"structural_matches_{site_id}_{rmsd_threshold}"
            matches = cached_results.get(site_id, [])
            QUERY_CACHE[cache_key] = matches
            CACHE_TIMESTAMPS[cache_key] = time.time()
                
        logger.info(f"Retrieved structural matches for {len(cached_results)} out of {len(site_ids)} sites")
        return cached_results
    except Exception as e:
        logger.error(f"Error finding batch structural matches: {e}")
        return cached_results  # Return whatever was found in cache even if the rest failed


def find_sequence_matches(site_id: str, min_similarity: float = 0.4) -> List[Dict]:
    """
    Find sequence similarity matches for a site.
    Uses composite idx_id1_sim and idx_id2_sim indexes for optimized lookup.
    Ensures that motif sequences are returned in uppercase.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        min_similarity: Minimum similarity score to include (0-1)
        
    Returns:
        List of dictionaries with match information
    """
    # Check cache first
    cache_key = f"sequence_matches_{site_id}_{min_similarity}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        PERFORMANCE_METRICS["cache_hits"] += 1
        matches = QUERY_CACHE[cache_key]
        # Ensure any motifs in cached matches are uppercase
        for match in matches:
            if 'motif' in match and match['motif']:
                match['motif'] = match['motif'].upper()
        return matches
    
    PERFORMANCE_METRICS["cache_misses"] += 1
    
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
                    
                    # Get motif if possible and ensure it's uppercase
                    try:
                        target_data = get_phosphosite_data(target_id)
                        if target_data and 'SITE_+/-7_AA' in target_data and target_data['SITE_+/-7_AA']:
                            match_dict['motif'] = target_data['SITE_+/-7_AA'].upper()
                        elif target_data and 'motif' in target_data and target_data['motif']:
                            match_dict['motif'] = target_data['motif'].upper()
                    except Exception as e:
                        # Continue without motif if we can't get it
                        logger.debug(f"Could not get motif for {target_id}: {e}")
                    
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
    Uses composite idx_id1_sim and idx_id2_sim indexes for optimized lookup.
    Ensures that motif sequences are returned in uppercase.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        min_similarity: Minimum similarity score to include (0-1)
        
    Returns:
        Dictionary mapping site IDs to lists of match dictionaries
    """
    if not site_ids:
        return {}
    
    # Check cache first for all sites
    cached_results = {}
    missing_sites = []
    
    for site_id in site_ids:
        cache_key = f"sequence_matches_{site_id}_{min_similarity}"
        if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
            PERFORMANCE_METRICS["cache_hits"] += 1
            matches = QUERY_CACHE[cache_key]
            # Ensure any motifs in cached matches are uppercase
            for match in matches:
                if 'motif' in match and match['motif']:
                    match['motif'] = match['motif'].upper()
            cached_results[site_id] = matches
        else:
            PERFORMANCE_METRICS["cache_misses"] += 1
            missing_sites.append(site_id)
    
    # If all results were in cache, return immediately
    if not missing_sites:
        return cached_results
        
    try:
        # Use the batch helper function
        query = QUERY_TEMPLATES["find_sequence_matches_batch"]
        df = execute_batch_query(
            query, 
            missing_sites, 
            batch_size=BATCH_SIZES["sequence_matches"],
            id_param_name="site_ids", 
            extra_params={"min_similarity": min_similarity}
        )
        
        if not df.empty:
            # First, build all the match dictionaries without motif information
            target_ids_to_fetch = set()  # Set to collect all target IDs
            match_data_by_query = {}  # Temporary storage for match data

            # Process results
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
                    
                # Add target_id to the set of IDs we need to fetch
                target_ids_to_fetch.add(target_id)
                    
                # Initialize storage for this query if needed
                if query_id not in match_data_by_query:
                    match_data_by_query[query_id] = []
                    
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
                    
                    # Add to the list for this query
                    match_data_by_query[query_id].append(match_dict)
            
            # Fetch all motif data in one batch query
            if target_ids_to_fetch:
                target_data_batch = get_phosphosites_batch(list(target_ids_to_fetch))
                
                # Now add motif information to each match and build final results
                for query_id, matches in match_data_by_query.items():
                    if query_id not in cached_results:
                        cached_results[query_id] = []
                        
                    for match_dict in matches:
                        target_id = match_dict["target_id"]
                        target_data = target_data_batch.get(target_id, {})
                        
                        # Add motif information if available
                        if target_data:
                            if 'SITE_+/-7_AA' in target_data and target_data['SITE_+/-7_AA']:
                                match_dict['motif'] = target_data['SITE_+/-7_AA'].upper()
                            elif 'motif' in target_data and target_data['motif']:
                                match_dict['motif'] = target_data['motif'].upper()
                        
                        # Add to the final results
                        cached_results[query_id].append(match_dict)
        
        # Update cache with new results
        for site_id in missing_sites:
            cache_key = f"sequence_matches_{site_id}_{min_similarity}"
            matches = cached_results.get(site_id, [])
            QUERY_CACHE[cache_key] = matches
            CACHE_TIMESTAMPS[cache_key] = time.time()
                
        logger.info(f"Retrieved sequence matches for {len(cached_results)} out of {len(site_ids)} sites")
        return cached_results
    except Exception as e:
        logger.error(f"Error finding batch sequence matches: {e}")
        return cached_results  # Return whatever was found in cache even if the rest failed

def get_kinase_scores(site_id: str, score_type: str = 'structure') -> Dict:
    """
    Get kinase scores for a specific site.
    Uses idx_node_structure or idx_node_sequence index for optimized lookup.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary with kinase names as keys and scores as values
    """
    # Check cache first
    cache_key = f"kinase_scores_{score_type}_{site_id}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        PERFORMANCE_METRICS["cache_hits"] += 1
        return QUERY_CACHE[cache_key]
    
    PERFORMANCE_METRICS["cache_misses"] += 1
    
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
    Uses idx_node_structure or idx_node_sequence index for optimized lookup.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary mapping site IDs to dictionaries of kinase scores
    """
    if not site_ids:
        return {}
    
    # Check cache first for all sites
    cached_results = {}
    missing_sites = []
    
    for site_id in site_ids:
        cache_key = f"kinase_scores_{score_type}_{site_id}"
        if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
            PERFORMANCE_METRICS["cache_hits"] += 1
            cached_results[site_id] = QUERY_CACHE[cache_key]
        else:
            PERFORMANCE_METRICS["cache_misses"] += 1
            missing_sites.append(site_id)
    
    # If all results were in cache, return immediately
    if not missing_sites:
        return cached_results
        
    try:
        # Select appropriate table and query based on score type
        if score_type.lower() == 'structure':
            query = QUERY_TEMPLATES["get_structure_kinase_scores_batch"]
        else:
            query = QUERY_TEMPLATES["get_sequence_kinase_scores_batch"]
            
        # Use the batch helper function
        df = execute_batch_query(
            query,
            missing_sites,
            batch_size=BATCH_SIZES["kinase_scores"],  # Use optimized batch size
            id_param_name="site_ids"
        )
        
        if not df.empty:
            # Process results
            for _, row in df.iterrows():
                site_id = str(row['node'])
                
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
                
                # Add to results and update cache
                cached_results[site_id] = result
                cache_key = f"kinase_scores_{score_type}_{site_id}"
                QUERY_CACHE[cache_key] = result
                CACHE_TIMESTAMPS[cache_key] = time.time()
            
        # Add empty results for missing sites
        for site_id in missing_sites:
            if site_id not in cached_results:
                empty_result = {'site_id': site_id, 'scores': {}, 'known_kinase': None}
                cached_results[site_id] = empty_result
                cache_key = f"kinase_scores_{score_type}_{site_id}"
                QUERY_CACHE[cache_key] = empty_result
                CACHE_TIMESTAMPS[cache_key] = time.time()
                
        logger.info(f"Retrieved kinase scores for {len(cached_results)} out of {len(site_ids)} sites")
        return cached_results
    except Exception as e:
        logger.error(f"Error getting batch kinase scores: {e}")
        return cached_results  # Return whatever was found in cache even if the rest failed

def find_sites_by_protein(uniprot_id: str) -> List[Dict]:
    """
    Find all phosphosites for a given protein by UniProt ID.
    Uses idx_uniprot_id index for optimized lookup.
    
    Args:
        uniprot_id: UniProt ID for the protein
        
    Returns:
        List of dictionaries with phosphosite data
    """
    # Check cache first
    cache_key = f"sites_by_protein_{uniprot_id}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        PERFORMANCE_METRICS["cache_hits"] += 1
        return QUERY_CACHE[cache_key]
    
    PERFORMANCE_METRICS["cache_misses"] += 1
    
    try:
        query = QUERY_TEMPLATES["find_sites_by_protein"]
        df = execute_query(query, {"uniprot_id": uniprot_id})
        
        sites = []
        if not df.empty:
            for _, row in df.iterrows():
                sites.append(row.to_dict())
        
        # Cache the result
        QUERY_CACHE[cache_key] = sites
        CACHE_TIMESTAMPS[cache_key] = time.time()
        
        return sites
    except Exception as e:
        logger.error(f"Error finding sites for protein {uniprot_id}: {e}")
        return []

def find_sites_by_protein_batch(uniprot_ids: List[str]) -> Dict[str, List[Dict]]:
    """
    Find all phosphosites for multiple proteins by UniProt ID.
    Uses idx_uniprot_id index for optimized lookup.
    
    Args:
        uniprot_ids: List of UniProt IDs
        
    Returns:
        Dictionary mapping UniProt IDs to lists of phosphosite data
    """
    if not uniprot_ids:
        return {}
    
    # Check cache first for all proteins
    cached_results = {}
    missing_proteins = []
    
    for uniprot_id in uniprot_ids:
        cache_key = f"sites_by_protein_{uniprot_id}"
        if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
            PERFORMANCE_METRICS["cache_hits"] += 1
            cached_results[uniprot_id] = QUERY_CACHE[cache_key]
        else:
            PERFORMANCE_METRICS["cache_misses"] += 1
            missing_proteins.append(uniprot_id)
    
    # If all results were in cache, return immediately
    if not missing_proteins:
        return cached_results
        
    try:
        # Use the batch helper function
        query = QUERY_TEMPLATES["find_sites_by_protein_batch"]
        df = execute_batch_query(
            query, 
            missing_proteins, 
            batch_size=BATCH_SIZES["phosphosites"],  # Use optimized batch size
            id_param_name="uniprot_ids"
        )
        
        if not df.empty:
            # Organize by protein
            for _, row in df.iterrows():
                row_dict = row.to_dict()
                uniprot_id = row_dict.get("uniprot_id")
                
                if uniprot_id not in cached_results:
                    cached_results[uniprot_id] = []
                
                cached_results[uniprot_id].append(row_dict)
        
        # Update cache with new results and add empty entries for proteins with no sites
        for uniprot_id in missing_proteins:
            cache_key = f"sites_by_protein_{uniprot_id}"
            sites = cached_results.get(uniprot_id, [])
            QUERY_CACHE[cache_key] = sites
            CACHE_TIMESTAMPS[cache_key] = time.time()
            
            # Ensure all proteins have an entry in results
            if uniprot_id not in cached_results:
                cached_results[uniprot_id] = []
                
        logger.info(f"Retrieved phosphosites for {len(cached_results)} proteins")
        return cached_results
    except Exception as e:
        logger.error(f"Error finding sites for protein batch: {e}")
        return cached_results  # Return whatever was found in cache even if the rest failed

def find_sites_by_gene(gene_symbol: str) -> List[Dict]:
    """
    Find all phosphosites for a given gene symbol.
    Uses idx_gene index for optimized lookup.
    
    Args:
        gene_symbol: Gene symbol
        
    Returns:
        List of dictionaries with phosphosite data
    """
    # Check cache first
    cache_key = f"sites_by_gene_{gene_symbol}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        PERFORMANCE_METRICS["cache_hits"] += 1
        return QUERY_CACHE[cache_key]
    
    PERFORMANCE_METRICS["cache_misses"] += 1
    
    try:
        query = QUERY_TEMPLATES["find_sites_by_gene"]
        df = execute_query(query, {"gene_symbol": gene_symbol})
        
        sites = []
        if not df.empty:
            for _, row in df.iterrows():
                sites.append(row.to_dict())
        
        # Cache the result
        QUERY_CACHE[cache_key] = sites
        CACHE_TIMESTAMPS[cache_key] = time.time()
        
        return sites
    except Exception as e:
        logger.error(f"Error finding sites for gene {gene_symbol}: {e}")
        return []

def find_sites_by_gene_batch(gene_symbols: List[str]) -> Dict[str, List[Dict]]:
    """
    Find all phosphosites for multiple gene symbols.
    Uses idx_gene index for optimized lookup.
    
    Args:
        gene_symbols: List of gene symbols
        
    Returns:
        Dictionary mapping gene symbols to lists of phosphosite data
    """
    if not gene_symbols:
        return {}
    
    # Check cache first for all genes
    cached_results = {}
    missing_genes = []
    
    for gene_symbol in gene_symbols:
        cache_key = f"sites_by_gene_{gene_symbol}"
        if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
            PERFORMANCE_METRICS["cache_hits"] += 1
            cached_results[gene_symbol] = QUERY_CACHE[cache_key]
        else:
            PERFORMANCE_METRICS["cache_misses"] += 1
            missing_genes.append(gene_symbol)
    
    # If all results were in cache, return immediately
    if not missing_genes:
        return cached_results
        
    try:
        # Use the batch helper function
        query = QUERY_TEMPLATES["find_sites_by_gene_batch"]
        df = execute_batch_query(
            query, 
            missing_genes, 
            batch_size=BATCH_SIZES["phosphosites"],  # Use optimized batch size
            id_param_name="gene_symbols"
        )
        
        if not df.empty:
            # Organize by gene
            for _, row in df.iterrows():
                row_dict = row.to_dict()
                gene_symbol = row_dict.get("GENE")
                
                if gene_symbol not in cached_results:
                    cached_results[gene_symbol] = []
                
                cached_results[gene_symbol].append(row_dict)
        
        # Update cache with new results and add empty entries for genes with no sites
        for gene_symbol in missing_genes:
            cache_key = f"sites_by_gene_{gene_symbol}"
            sites = cached_results.get(gene_symbol, [])
            QUERY_CACHE[cache_key] = sites
            CACHE_TIMESTAMPS[cache_key] = time.time()
            
            # Ensure all genes have an entry in results
            if gene_symbol not in cached_results:
                cached_results[gene_symbol] = []
                
        logger.info(f"Retrieved phosphosites for {len(cached_results)} genes")
        return cached_results
    except Exception as e:
        logger.error(f"Error finding sites for gene batch: {e}")
        return cached_results  # Return whatever was found in cache even if the rest failed

def find_structural_matches_by_target(target_id: str, rmsd_threshold: float = 5.0) -> List[Dict]:
    """
    Find structural matches where the given site is the target.
    Uses idx_target index for optimized lookup.
    
    Args:
        target_id: Target site ID in format 'UniProtID_ResidueNumber'
        rmsd_threshold: Maximum RMSD value for matches
        
    Returns:
        List of dictionaries with match information
    """
    # Check cache first
    cache_key = f"structural_matches_by_target_{target_id}_{rmsd_threshold}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        PERFORMANCE_METRICS["cache_hits"] += 1
        return QUERY_CACHE[cache_key]
    
    PERFORMANCE_METRICS["cache_misses"] += 1
    
    try:
        query = """
            SELECT * FROM `Structural_Similarity_Edges`
            WHERE `Target` = :target_id AND `RMSD` <= :rmsd_threshold
            ORDER BY `RMSD` ASC
            /* Uses index: idx_target */
        """
        df = execute_query(query, {
            "target_id": target_id,
            "rmsd_threshold": rmsd_threshold
        })
        
        matches = []
        if not df.empty:
            # Process matches
            for _, row in df.iterrows():
                row_dict = row.to_dict()
                
                # Parse query info
                query_id = row_dict.get("Query", "")
                query_parts = query_id.split('_')
                
                if len(query_parts) >= 2:
                    query_uniprot = query_parts[0]
                    query_site_raw = query_parts[1]
                    
                    # Create match dictionary
                    match_dict = {
                        "target_id": target_id,
                        "query_id": query_id,
                        "query_uniprot": query_uniprot,
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
        logger.error(f"Error finding structural matches by target: {e}")
        return []
    

def get_structural_annotations(site_id: str) -> Optional[Dict]:
    """
    Get structural annotations for a phosphosite.
    Uses primary key for optimized lookup.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary with structural annotations or None if not found
    """
    # Check cache first
    cache_key = f"structural_annotations_{site_id}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        PERFORMANCE_METRICS["cache_hits"] += 1
        return QUERY_CACHE[cache_key]
    
    PERFORMANCE_METRICS["cache_misses"] += 1
    
    try:
        query = QUERY_TEMPLATES["get_structural_annotations"]
        df = execute_query(query, {"site_id": site_id})
        
        if df.empty:
            logger.warning(f"No structural annotations found for {site_id}")
            # Cache the negative result too
            QUERY_CACHE[cache_key] = None
            CACHE_TIMESTAMPS[cache_key] = time.time()
            return None
            
        # Convert first row to dictionary
        result = df.iloc[0].to_dict()
        
        # Ensure motif is uppercase if present
        if 'Motif' in result and result['Motif']:
            result['Motif'] = str(result['Motif']).upper()
        
        # Cache the result
        QUERY_CACHE[cache_key] = result
        CACHE_TIMESTAMPS[cache_key] = time.time()
        
        return result
    except Exception as e:
        logger.error(f"Error getting structural annotations: {e}")
        return None

def get_structural_annotations_batch(site_ids: List[str]) -> Dict[str, Dict]:
    """
    Get structural annotations for multiple sites in a batch.
    Uses primary key for optimized lookup.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary mapping site IDs to structural annotation dictionaries
    """
    if not site_ids:
        return {}
    
    # Check cache first for all sites
    cached_results = {}
    missing_sites = []
    
    for site_id in site_ids:
        cache_key = f"structural_annotations_{site_id}"
        if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
            PERFORMANCE_METRICS["cache_hits"] += 1
            cached_value = QUERY_CACHE[cache_key]
            if cached_value is not None:  # Don't include None values
                # Ensure motif is uppercase if present
                if 'Motif' in cached_value and cached_value['Motif']:
                    cached_value['Motif'] = str(cached_value['Motif']).upper()
                cached_results[site_id] = cached_value
        else:
            PERFORMANCE_METRICS["cache_misses"] += 1
            missing_sites.append(site_id)
    
    # If all results were in cache, return immediately
    if not missing_sites:
        return cached_results
        
    try:
        # Use the batch helper function for missing sites
        query = QUERY_TEMPLATES["get_structural_annotations_batch"]
        df = execute_batch_query(
            query, 
            missing_sites, 
            batch_size=BATCH_SIZES["phosphosites"],  # Use optimized batch size
            id_param_name="site_ids"
        )
        
        if not df.empty:
            # Convert to dictionary of dictionaries
            for _, row in df.iterrows():
                row_dict = row.to_dict()
                if "PhosphositeID" in row_dict:
                    site_id = row_dict["PhosphositeID"]
                    
                    # Ensure motif is uppercase if present
                    if 'Motif' in row_dict and row_dict['Motif']:
                        row_dict['Motif'] = str(row_dict['Motif']).upper()
                    
                    # Add to results
                    cached_results[site_id] = row_dict
                    # Also update cache
                    cache_key = f"structural_annotations_{site_id}"
                    QUERY_CACHE[cache_key] = row_dict
                    CACHE_TIMESTAMPS[cache_key] = time.time()
        
        # Cache negative results for missing sites
        for site_id in missing_sites:
            if site_id not in cached_results:
                cache_key = f"structural_annotations_{site_id}"
                QUERY_CACHE[cache_key] = None
                CACHE_TIMESTAMPS[cache_key] = time.time()
                
        logger.info(f"Retrieved structural annotations for {len(cached_results)} out of {len(site_ids)} phosphosites")
        return cached_results
    except Exception as e:
        logger.error(f"Error getting batch structural annotations: {e}")
        return cached_results  # Return whatever was found in cache even if the rest failed

def find_structural_annotations_by_uniprot(uniprot_id: str) -> List[Dict]:
    """
    Find all structural annotations for a given protein by UniProt ID.
    Uses idx_uniprot_residue index for optimized lookup.
    
    Args:
        uniprot_id: UniProt ID for the protein
        
    Returns:
        List of dictionaries with structural annotation data
    """
    # Check cache first
    cache_key = f"structural_annotations_by_uniprot_{uniprot_id}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        PERFORMANCE_METRICS["cache_hits"] += 1
        return QUERY_CACHE[cache_key]
    
    PERFORMANCE_METRICS["cache_misses"] += 1
    
    try:
        query = QUERY_TEMPLATES["find_structural_annotations_by_uniprot"]
        df = execute_query(query, {"uniprot_id": uniprot_id})
        
        annotations = []
        if not df.empty:
            for _, row in df.iterrows():
                row_dict = row.to_dict()
                # Ensure motif is uppercase if present
                if 'Motif' in row_dict and row_dict['Motif']:
                    row_dict['Motif'] = str(row_dict['Motif']).upper()
                annotations.append(row_dict)
        
        # Cache the result
        QUERY_CACHE[cache_key] = annotations
        CACHE_TIMESTAMPS[cache_key] = time.time()
        
        return annotations
    except Exception as e:
        logger.error(f"Error finding structural annotations for protein {uniprot_id}: {e}")
        return []

def get_comprehensive_site_data(site_id: str) -> Dict:
    """
    Get comprehensive data for a specific site by combining STY_Structural_Annotations 
    with Phosphosite_Supplementary_Data.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        
    Returns:
        Integrated dictionary with all site data
    """
    # Get structural annotations
    structural_data = get_structural_annotations(site_id)
    
    # Get phosphosite supplementary data
    supplementary_data = get_phosphosite_data(site_id)
    
    # Create integrated data dictionary
    integrated_data = {}
    
    # Add structural data if available
    if structural_data:
        integrated_data.update(structural_data)
    
    # Add supplementary data if available (will overwrite duplicate keys)
    if supplementary_data:
        for key, value in supplementary_data.items():
            # Skip None values when supplementary data has them
            if value is not None:
                integrated_data[key] = value
    
    # Ensure site ID is included
    if 'PhosphositeID' not in integrated_data:
        integrated_data['PhosphositeID'] = site_id
    
    return integrated_data

def get_comprehensive_site_data_batch(site_ids: List[str]) -> Dict[str, Dict]:
    """
    Get comprehensive data for multiple sites by combining STY_Structural_Annotations
    with Phosphosite_Supplementary_Data in an efficient batch operation.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary mapping site IDs to comprehensive data dictionaries
    """
    if not site_ids:
        return {}
    
    result = {}
    
    # Get structural annotations in batch
    structural_batch = get_structural_annotations_batch(site_ids)
    
    # Get phosphosite supplementary data in batch
    supplementary_batch = get_phosphosites_batch(site_ids)
    
    # Combine the data for each site
    for site_id in site_ids:
        integrated_data = {}
        
        # Add structural data if available
        if site_id in structural_batch:
            integrated_data.update(structural_batch[site_id])
        
        # Add supplementary data if available
        if site_id in supplementary_batch:
            for key, value in supplementary_batch[site_id].items():
                # Skip None values
                if value is not None:
                    integrated_data[key] = value
        
        # Only add to results if we have some data
        if integrated_data:
            # Ensure site ID is included
            if 'PhosphositeID' not in integrated_data:
                integrated_data['PhosphositeID'] = site_id
                
            result[site_id] = integrated_data
    
    return result



def get_cantley_st_kinase_scores(site_id: str) -> Optional[Dict]:
    """
    Get Cantley lab S/T kinase scores for a specific phosphosite.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary with kinase scores or None if not found
    """
    # Check cache first
    cache_key = f"cantley_st_scores_{site_id}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        PERFORMANCE_METRICS["cache_hits"] += 1
        return QUERY_CACHE[cache_key]
    
    PERFORMANCE_METRICS["cache_misses"] += 1
    
    try:
        query = QUERY_TEMPLATES["get_cantley_st_kinase_scores"]
        df = execute_query(query, {"site_id": site_id})
        
        if df.empty:
            logger.warning(f"No Cantley S/T kinase scores found for {site_id}")
            # Cache the negative result too
            QUERY_CACHE[cache_key] = None
            CACHE_TIMESTAMPS[cache_key] = time.time()
            return None
            
        # Convert first row to dictionary
        result = df.iloc[0].to_dict()
        
        # Convert scores to a separate dictionary
        scores = {}
        for key, value in result.items():
            if key not in ['SiteID', 'Motif'] and pd.notna(value):
                scores[key] = float(value)
        
        # Create final result with organized structure
        final_result = {
            'site_id': site_id,
            'motif': result.get('Motif', ''),
            'scores': scores
        }
        
        # Cache the result
        QUERY_CACHE[cache_key] = final_result
        CACHE_TIMESTAMPS[cache_key] = time.time()
        
        return final_result
    except Exception as e:
        logger.error(f"Error getting Cantley S/T kinase scores: {e}")
        return None

def get_cantley_st_kinase_scores_batch(site_ids: List[str]) -> Dict[str, Dict]:
    """
    Get Cantley lab S/T kinase scores for multiple phosphosites in a batch.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary mapping site IDs to dictionaries of scores
    """
    if not site_ids:
        return {}
    
    # Check cache first for all sites
    cached_results = {}
    missing_sites = []
    
    for site_id in site_ids:
        cache_key = f"cantley_st_scores_{site_id}"
        if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
            PERFORMANCE_METRICS["cache_hits"] += 1
            cached_value = QUERY_CACHE[cache_key]
            if cached_value is not None:  # Don't include None values
                cached_results[site_id] = cached_value
        else:
            PERFORMANCE_METRICS["cache_misses"] += 1
            missing_sites.append(site_id)
    
    # If all results were in cache, return immediately
    if not missing_sites:
        return cached_results
        
    try:
        # Use the batch helper function for missing sites
        query = QUERY_TEMPLATES["get_cantley_st_kinase_scores_batch"]
        df = execute_batch_query(
            query, 
            missing_sites, 
            batch_size=BATCH_SIZES["phosphosites"],  # Use the same batch size as phosphosites
            id_param_name="site_ids"
        )
        
        if not df.empty:
            # Convert to dictionary of dictionaries
            for _, row in df.iterrows():
                row_dict = row.to_dict()
                if "SiteID" in row_dict:
                    site_id = row_dict["SiteID"]
                    
                    # Convert scores to a separate dictionary
                    scores = {}
                    for key, value in row_dict.items():
                        if key not in ['SiteID', 'Motif'] and pd.notna(value):
                            scores[key] = float(value)
                    
                    # Create final result with organized structure
                    final_result = {
                        'site_id': site_id,
                        'motif': row_dict.get('Motif', ''),
                        'scores': scores
                    }
                    
                    # Add to results
                    cached_results[site_id] = final_result
                    # Also update cache
                    cache_key = f"cantley_st_scores_{site_id}"
                    QUERY_CACHE[cache_key] = final_result
                    CACHE_TIMESTAMPS[cache_key] = time.time()
        
        # Cache negative results for missing sites
        for site_id in missing_sites:
            if site_id not in cached_results:
                cache_key = f"cantley_st_scores_{site_id}"
                QUERY_CACHE[cache_key] = None
                CACHE_TIMESTAMPS[cache_key] = time.time()
                
        logger.info(f"Retrieved Cantley S/T kinase scores for {len(cached_results)} out of {len(site_ids)} phosphosites")
        return cached_results
    except Exception as e:
        logger.error(f"Error getting batch Cantley S/T kinase scores: {e}")
        return cached_results  # Return whatever was found in cache even if the rest failed

def get_cantley_y_kinase_scores(site_id: str) -> Optional[Dict]:
    """
    Get Cantley lab Y kinase scores for a specific phosphosite.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary with kinase scores or None if not found
    """
    # Check cache first
    cache_key = f"cantley_y_scores_{site_id}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        PERFORMANCE_METRICS["cache_hits"] += 1
        return QUERY_CACHE[cache_key]
    
    PERFORMANCE_METRICS["cache_misses"] += 1
    
    try:
        query = QUERY_TEMPLATES["get_cantley_y_kinase_scores"]
        df = execute_query(query, {"site_id": site_id})
        
        if df.empty:
            logger.warning(f"No Cantley Y kinase scores found for {site_id}")
            # Cache the negative result too
            QUERY_CACHE[cache_key] = None
            CACHE_TIMESTAMPS[cache_key] = time.time()
            return None
            
        # Convert first row to dictionary
        result = df.iloc[0].to_dict()
        
        # Convert scores to a separate dictionary
        scores = {}
        for key, value in result.items():
            if key not in ['SiteID', 'Motif'] and pd.notna(value):
                scores[key] = float(value)
        
        # Create final result with organized structure
        final_result = {
            'site_id': site_id,
            'motif': result.get('Motif', ''),
            'scores': scores
        }
        
        # Cache the result
        QUERY_CACHE[cache_key] = final_result
        CACHE_TIMESTAMPS[cache_key] = time.time()
        
        return final_result
    except Exception as e:
        logger.error(f"Error getting Cantley Y kinase scores: {e}")
        return None

def get_cantley_y_kinase_scores_batch(site_ids: List[str]) -> Dict[str, Dict]:
    """
    Get Cantley lab Y kinase scores for multiple phosphosites in a batch.
    
    Args:
        site_ids: List of site IDs in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary mapping site IDs to dictionaries of scores
    """
    if not site_ids:
        return {}
    
    # Check cache first for all sites
    cached_results = {}
    missing_sites = []
    
    for site_id in site_ids:
        cache_key = f"cantley_y_scores_{site_id}"
        if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
            PERFORMANCE_METRICS["cache_hits"] += 1
            cached_value = QUERY_CACHE[cache_key]
            if cached_value is not None:  # Don't include None values
                cached_results[site_id] = cached_value
        else:
            PERFORMANCE_METRICS["cache_misses"] += 1
            missing_sites.append(site_id)
    
    # If all results were in cache, return immediately
    if not missing_sites:
        return cached_results
        
    try:
        # Use the batch helper function for missing sites
        query = QUERY_TEMPLATES["get_cantley_y_kinase_scores_batch"]
        df = execute_batch_query(
            query, 
            missing_sites, 
            batch_size=BATCH_SIZES["phosphosites"],  # Use the same batch size as phosphosites
            id_param_name="site_ids"
        )
        
        if not df.empty:
            # Convert to dictionary of dictionaries
            for _, row in df.iterrows():
                row_dict = row.to_dict()
                if "SiteID" in row_dict:
                    site_id = row_dict["SiteID"]
                    
                    # Convert scores to a separate dictionary
                    scores = {}
                    for key, value in row_dict.items():
                        if key not in ['SiteID', 'Motif'] and pd.notna(value):
                            scores[key] = float(value)
                    
                    # Create final result with organized structure
                    final_result = {
                        'site_id': site_id,
                        'motif': row_dict.get('Motif', ''),
                        'scores': scores
                    }
                    
                    # Add to results
                    cached_results[site_id] = final_result
                    # Also update cache
                    cache_key = f"cantley_y_scores_{site_id}"
                    QUERY_CACHE[cache_key] = final_result
                    CACHE_TIMESTAMPS[cache_key] = time.time()
        
        # Cache negative results for missing sites
        for site_id in missing_sites:
            if site_id not in cached_results:
                cache_key = f"cantley_y_scores_{site_id}"
                QUERY_CACHE[cache_key] = None
                CACHE_TIMESTAMPS[cache_key] = time.time()
                
        logger.info(f"Retrieved Cantley Y kinase scores for {len(cached_results)} out of {len(site_ids)} phosphosites")
        return cached_results
    except Exception as e:
        logger.error(f"Error getting batch Cantley Y kinase scores: {e}")
        return cached_results  # Return whatever was found in cache even if the rest failed

def get_cantley_top_kinases(site_id: str, residue_type: str = None, top_n: int = 5) -> Dict:
    """
    Get top predicted kinases for a phosphosite from Cantley lab scores.
    Automatically determines whether to use S/T or Y kinases based on the residue type.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        residue_type: Optional residue type ('S', 'T', or 'Y'). If None, determined from site_id.
        top_n: Number of top kinases to return
        
    Returns:
        Dictionary with top kinases and their scores
    """
    # Determine residue type if not provided
    if residue_type is None:
        # Try to extract from site_id
        parts = site_id.split('_')
        if len(parts) >= 2:
            # Check for residue type in second part
            match = re.match(r'([STY])(\d+)', parts[1])
            if match:
                residue_type = match.group(1)
            else:
                # If not found, check if we have phosphosite data
                site_data = get_phosphosite_data(site_id)
                if site_data and 'Residue' in site_data:
                    residue_type = site_data['Residue']
                else:
                    # Default to S/T if we can't determine
                    residue_type = 'S'
        else:
            # Default to S/T if malformed site_id
            residue_type = 'S'
    
    # Get appropriate kinase scores
    if residue_type.upper() == 'Y':
        scores_data = get_cantley_y_kinase_scores(site_id)
    else:  # S or T
        scores_data = get_cantley_st_kinase_scores(site_id)
    
    if not scores_data or 'scores' not in scores_data:
        return {
            'site_id': site_id,
            'residue_type': residue_type.upper(),
            'top_kinases': [],
            'motif': ''
        }
    
    # Sort kinases by score and get top N
    sorted_kinases = sorted(
        [(kinase, score) for kinase, score in scores_data['scores'].items()],
        key=lambda x: x[1],
        reverse=True
    )
    
    top_kinases = sorted_kinases[:top_n]
    
    return {
        'site_id': site_id,
        'residue_type': residue_type.upper(),
        'top_kinases': [{'kinase': k, 'score': s} for k, s in top_kinases],
        'motif': scores_data.get('motif', '')
    }

def get_sites_by_cantley_kinase(kinase_name: str, score_threshold: float = 0.9, limit: int = 100) -> List[Dict]:
    """
    Find sites with high prediction scores for a specific kinase using Cantley lab scores.
    Automatically determines whether to look in S/T or Y table based on kinase name.
    
    Args:
        kinase_name: Name of the kinase to look for
        score_threshold: Minimum score to include (0-1)
        limit: Maximum number of sites to return
        
    Returns:
        List of sites with high scores for the specified kinase
    """
    # Check cache first
    cache_key = f"sites_by_cantley_{kinase_name}_{score_threshold}_{limit}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        PERFORMANCE_METRICS["cache_hits"] += 1
        return QUERY_CACHE[cache_key]
    
    PERFORMANCE_METRICS["cache_misses"] += 1
    
    # Determine which table to query based on kinase name
    # List of kinases that exist only in the Y table
    y_specific_kinases = [
        "ABL", "ACK", "ALK", "ARG", "AXL", "BLK", "BMPR2_TYR", "BRK", "BTK", "CSFR", 
        "CSK", "CTK", "DDR1", "DDR2", "EGFR", "EPHA1", "EPHA2", "EPHA3", "EPHA4", "EPHA5", 
        "EPHA6", "EPHA7", "EPHA8", "EPHB1", "EPHB2", "EPHB3", "EPHB4", "ETK", "FAK", "FER", 
        "FES", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FGR", "FLT3", "FRK", "FYN", "HCK", 
        "HER2", "HER4", "IGF1R", "INSR", "IRR", "ITK", "JAK1", "JAK2", "JAK3", "KIT", 
        "LCK", "LIMK1_TYR", "LIMK2_TYR", "LTK", "LYN", "MER", "MET", "MKK4_TYR", "MKK6_TYR", 
        "MKK7_TYR", "MST1R", "MUSK", "MYT1_TYR", "NEK10_TYR", "PDGFRA", "PDGFRB", "PDHK1_TYR", 
        "PDHK3_TYR", "PDHK4_TYR", "PINK1_TYR", "PYK2", "RET", "ROS", "SRC", "SRMS", "SYK", 
        "TEC", "TESK1_TYR", "TIE2", "TNK1", "TNNI3K_TYR", "TRKA", "TRKB", "TRKC", "TXK", 
        "TYK2", "TYRO3", "VEGFR1", "VEGFR2", "VEGFR3", "WEE1_TYR", "YES", "ZAP70"
    ]
    
    # Check if this kinase is Y-specific
    is_y_kinase = kinase_name in y_specific_kinases or kinase_name.endswith("_TYR")
    
    try:
        if is_y_kinase:
            query = QUERY_TEMPLATES["get_sites_by_cantley_y_kinase"]
        else:
            query = QUERY_TEMPLATES["get_sites_by_cantley_st_kinase"]
        
        # Replace the kinase placeholder in the query
        modified_query = query.replace(':kinase_name', kinase_name)
        
        df = execute_query(modified_query, {
            "score_threshold": score_threshold,
            "limit": limit
        })
        
        results = []
        if not df.empty:
            for _, row in df.iterrows():
                site_id = row['SiteID']
                motif = row.get('Motif', '')
                score = float(row.get(kinase_name, 0.0))
                
                results.append({
                    'site_id': site_id,
                    'motif': motif,
                    'score': score
                })
        
        # Cache the result
        QUERY_CACHE[cache_key] = results
        CACHE_TIMESTAMPS[cache_key] = time.time()
        
        return results
    except Exception as e:
        logger.error(f"Error finding sites by Cantley kinase {kinase_name}: {e}")
        return []

def get_cantley_kinase_names(residue_type: str = 'all') -> Dict[str, List[str]]:
    """
    Get list of available kinase names in the Cantley datasets.
    
    Args:
        residue_type: Type of kinases to return - 'S', 'T', 'Y', or 'all'
        
    Returns:
        Dictionary with kinase names by type
    """
    # Check cache first
    cache_key = f"cantley_kinase_names_{residue_type}"
    if cache_key in QUERY_CACHE and is_cache_valid(cache_key):
        PERFORMANCE_METRICS["cache_hits"] += 1
        return QUERY_CACHE[cache_key]
    
    PERFORMANCE_METRICS["cache_misses"] += 1
    
    try:
        st_kinases = []
        y_kinases = []
        
        # Get S/T kinases if needed
        if residue_type.lower() in ['s', 't', 'all']:
            query = "SHOW COLUMNS FROM `Cantley_Kinome_Scores_All_STs` WHERE Field != 'SiteID' AND Field != 'Motif'"
            df = execute_query(query)
            if not df.empty:
                st_kinases = df['Field'].tolist()
        
        # Get Y kinases if needed
        if residue_type.lower() in ['y', 'all']:
            query = "SHOW COLUMNS FROM `Cantley_Kinome_Scores_All_Ys` WHERE Field != 'SiteID' AND Field != 'Motif'"
            df = execute_query(query)
            if not df.empty:
                y_kinases = df['Field'].tolist()
        
        result = {
            'S/T_kinases': st_kinases,
            'Y_kinases': y_kinases
        }
        
        # Cache the result
        QUERY_CACHE[cache_key] = result
        CACHE_TIMESTAMPS[cache_key] = time.time()
        
        return result
    except Exception as e:
        logger.error(f"Error getting Cantley kinase names: {e}")
        return {'S/T_kinases': [], 'Y_kinases': []}

def clear_cache():
    """Clear the query cache."""
    global QUERY_CACHE, CACHE_TIMESTAMPS
    QUERY_CACHE = {}
    CACHE_TIMESTAMPS = {}
    logger.info("Query cache cleared")

def reset_performance_metrics():
    """Reset performance tracking metrics."""
    global PERFORMANCE_METRICS
    PERFORMANCE_METRICS = {
        "query_times": {},
        "cache_hits": 0,
        "cache_misses": 0,
        "db_errors": 0
    }
    logger.info("Performance metrics reset")

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
        'tables_available': list(TABLES.keys()),
        'performance': {
            'cache_hit_rate': 0,
            'average_query_time': 0,
            'errors': PERFORMANCE_METRICS["db_errors"]
        }
    }
    
    # Calculate cache hit rate
    total_cache_requests = PERFORMANCE_METRICS["cache_hits"] + PERFORMANCE_METRICS["cache_misses"]
    if total_cache_requests > 0:
        stats['performance']['cache_hit_rate'] = PERFORMANCE_METRICS["cache_hits"] / total_cache_requests
    
    # Calculate average query time
    if PERFORMANCE_METRICS["query_times"]:
        avg_time = sum(PERFORMANCE_METRICS["query_times"].values()) / len(PERFORMANCE_METRICS["query_times"])
        stats['performance']['average_query_time'] = avg_time
    
    try:
        # Execute count queries for each table
        for table_name in TABLES:
            query = f'SELECT COUNT(*) FROM `{table_name}`'
            df = execute_query(query)
            if not df.empty:
                stats['table_counts'][table_name] = int(df.iloc[0, 0])
            else:
                stats['table_counts'][table_name] = 0
    except Exception as e:
        logger.error(f"Error getting database stats: {e}")
    
    return stats

def print_structural_edges_sample(sample_size=100):
    """
    Print a sample of the first rows from the Structural_Similarity_Edges table.
    Useful for debugging and verification after table updates.
    
    Args:
        sample_size: Number of rows to sample (default: 100)
    """
    try:
        # Ensure DB connection is initialized
        if engine is None:
            init_db()
            
        # Query to get the first N rows
        sample_query = f"""
        SELECT * FROM Structural_Similarity_Edges
        LIMIT {sample_size}
        """
        
        df = execute_query(sample_query)
        
        if df.empty:
            logger.warning("No data found in Structural_Similarity_Edges table")
            return
            
        # Log basic stats
        logger.info(f"Structural_Similarity_Edges table sample ({len(df)} rows):")
        logger.info(f"Columns: {list(df.columns)}")
        
        # Get a count of total rows
        count_query = "SELECT COUNT(*) as total FROM Structural_Similarity_Edges"
        count_df = execute_query(count_query)
        total_rows = count_df.iloc[0]['total'] if not count_df.empty else "unknown"
        logger.info(f"Total rows in table: {total_rows}")
        
        # For each row in the sample, print a summary
        for i, row in df.iterrows():
            # Limit to first 10 rows of detailed output to avoid log spam
            if i < 10:
                logger.info(f"Row {i+1}:")
                logger.info(f"  Query: {row.get('Query', 'N/A')}")
                logger.info(f"  Target: {row.get('Target', 'N/A')}")
                logger.info(f"  RMSD: {row.get('RMSD', 'N/A')}")
                
                # Print a few of the new columns if they exist
                for col in ['GlobalScore', 'QueryNeighborhood', 'MappingRatio']:
                    if col in row:
                        value = row[col]
                        # Truncate long string values
                        if isinstance(value, str) and len(value) > 50:
                            value = value[:50] + "..."
                        logger.info(f"  {col}: {value}")
                logger.info("  ---")
                
        # Print a simple summary for the remaining rows
        if len(df) > 10:
            logger.info(f"... and {len(df) - 10} more rows in sample")
            
    except Exception as e:
        logger.error(f"Error printing Structural_Similarity_Edges sample: {e}")


# Initialize database on module import
clear_cache()
init_db()
#clear_cache()
