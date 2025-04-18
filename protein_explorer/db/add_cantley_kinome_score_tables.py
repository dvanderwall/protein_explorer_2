#!/usr/bin/env python
"""
Script to bulk load Cantley kinome score data into the MySQL database.
Uses LOAD DATA INFILE for optimal performance.

Usage:
  python add_cantley_kinome_tables.py
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
        logging.FileHandler('cantley_load.log')
    ]
)

logger = logging.getLogger(__name__)

# Database connection parameters
DB_USER = os.environ.get('DB_USER', 'root')
DB_PASS = os.environ.get('DB_PASS', '@Bismark6')
DB_HOST = os.environ.get('DB_HOST', '35.245.113.195')
DB_PORT = os.environ.get('DB_PORT', '3306')
DB_NAME = os.environ.get('DB_NAME', 'kinoplex-db')

# Path to the CSV files
ST_FILE_PATH = "F:/Kinome/motif_scores_percentile_STs.csv"
Y_FILE_PATH = "F:/Kinome/motif_scores_percentile_Ys.csv"

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

def create_st_table(cursor, connection):
    """Create the S/T kinase scores table with proper structure and indexes."""
    logger.info("Creating S/T kinase scores table...")
    
    create_table_query = """
    CREATE TABLE IF NOT EXISTS `Cantley_Kinome_Scores_All_STs` (
        `SiteID` VARCHAR(50) NOT NULL,
        `Motif` VARCHAR(50),
        `AAK1` FLOAT, `ACVR2A` FLOAT, `ACVR2B` FLOAT, `AKT1` FLOAT, `AKT2` FLOAT, 
        `AKT3` FLOAT, `ALK2` FLOAT, `ALK4` FLOAT, `ALPHAK3` FLOAT, `AMPKA1` FLOAT, 
        `AMPKA2` FLOAT, `ANKRD3` FLOAT, `ASK1` FLOAT, `ATM` FLOAT, `ATR` FLOAT, 
        `AURA` FLOAT, `AURB` FLOAT, `AURC` FLOAT, `BCKDK` FLOAT, `BIKE` FLOAT, 
        `BMPR1A` FLOAT, `BMPR1B` FLOAT, `BMPR2` FLOAT, `BRAF` FLOAT, `BRSK1` FLOAT, 
        `BRSK2` FLOAT, `BUB1` FLOAT, `CAMK1A` FLOAT, `CAMK1B` FLOAT, `CAMK1D` FLOAT, 
        `CAMK1G` FLOAT, `CAMK2A` FLOAT, `CAMK2B` FLOAT, `CAMK2D` FLOAT, `CAMK2G` FLOAT, 
        `CAMK4` FLOAT, `CAMKK1` FLOAT, `CAMKK2` FLOAT, `CAMLCK` FLOAT, `CDC7` FLOAT, 
        `CDK1` FLOAT, `CDK10` FLOAT, `CDK12` FLOAT, `CDK13` FLOAT, `CDK14` FLOAT, 
        `CDK16` FLOAT, `CDK17` FLOAT, `CDK18` FLOAT, `CDK19` FLOAT, `CDK2` FLOAT, 
        `CDK3` FLOAT, `CDK4` FLOAT, `CDK5` FLOAT, `CDK6` FLOAT, `CDK7` FLOAT, 
        `CDK8` FLOAT, `CDK9` FLOAT, `CDKL1` FLOAT, `CDKL5` FLOAT, `CHAK1` FLOAT, 
        `CHAK2` FLOAT, `CHK1` FLOAT, `CHK2` FLOAT, `CK1A` FLOAT, `CK1A2` FLOAT, 
        `CK1D` FLOAT, `CK1E` FLOAT, `CK1G1` FLOAT, `CK1G2` FLOAT, `CK1G3` FLOAT, 
        `CK2A1` FLOAT, `CK2A2` FLOAT, `CLK1` FLOAT, `CLK2` FLOAT, `CLK3` FLOAT, 
        `CLK4` FLOAT, `COT` FLOAT, `CRIK` FLOAT, `DAPK1` FLOAT, `DAPK2` FLOAT, 
        `DAPK3` FLOAT, `DCAMKL1` FLOAT, `DCAMKL2` FLOAT, `DLK` FLOAT, `DMPK1` FLOAT, 
        `DNAPK` FLOAT, `DRAK1` FLOAT, `DSTYK` FLOAT, `DYRK1A` FLOAT, `DYRK1B` FLOAT, 
        `DYRK2` FLOAT, `DYRK3` FLOAT, `DYRK4` FLOAT, `EEF2K` FLOAT, `ERK1` FLOAT, 
        `ERK2` FLOAT, `ERK5` FLOAT, `ERK7` FLOAT, `FAM20C` FLOAT, `GAK` FLOAT, 
        `GCK` FLOAT, `GCN2` FLOAT, `GRK1` FLOAT, `GRK2` FLOAT, `GRK3` FLOAT, 
        `GRK4` FLOAT, `GRK5` FLOAT, `GRK6` FLOAT, `GRK7` FLOAT, `GSK3A` FLOAT, 
        `GSK3B` FLOAT, `HASPIN` FLOAT, `HGK` FLOAT, `HIPK1` FLOAT, `HIPK2` FLOAT, 
        `HIPK3` FLOAT, `HIPK4` FLOAT, `HPK1` FLOAT, `HRI` FLOAT, `HUNK` FLOAT, 
        `ICK` FLOAT, `IKKA` FLOAT, `IKKB` FLOAT, `IKKE` FLOAT, `IRAK1` FLOAT, 
        `IRAK4` FLOAT, `IRE1` FLOAT, `IRE2` FLOAT, `JNK1` FLOAT, `JNK2` FLOAT, 
        `JNK3` FLOAT, `KHS1` FLOAT, `KHS2` FLOAT, `KIS` FLOAT, `LATS1` FLOAT, 
        `LATS2` FLOAT, `LKB1` FLOAT, `LOK` FLOAT, `LRRK2` FLOAT, `MAK` FLOAT, 
        `MAP3K15` FLOAT, `MAPKAPK2` FLOAT, `MAPKAPK3` FLOAT, `MAPKAPK5` FLOAT, 
        `MARK1` FLOAT, `MARK2` FLOAT, `MARK3` FLOAT, `MARK4` FLOAT, `MASTL` FLOAT, 
        `MEK1` FLOAT, `MEK2` FLOAT, `MEK5` FLOAT, `MEKK1` FLOAT, `MEKK2` FLOAT, 
        `MEKK3` FLOAT, `MEKK6` FLOAT, `MELK` FLOAT, `MINK` FLOAT, `MLK1` FLOAT, 
        `MLK2` FLOAT, `MLK3` FLOAT, `MLK4` FLOAT, `MNK1` FLOAT, `MNK2` FLOAT, 
        `MOK` FLOAT, `MOS` FLOAT, `MPSK1` FLOAT, `MRCKA` FLOAT, `MRCKB` FLOAT, 
        `MSK1` FLOAT, `MSK2` FLOAT, `MST1` FLOAT, `MST2` FLOAT, `MST3` FLOAT, 
        `MST4` FLOAT, `MTOR` FLOAT, `MYLK4` FLOAT, `MYO3A` FLOAT, `MYO3B` FLOAT, 
        `NDR1` FLOAT, `NDR2` FLOAT, `NEK1` FLOAT, `NEK11` FLOAT, `NEK2` FLOAT, 
        `NEK3` FLOAT, `NEK4` FLOAT, `NEK5` FLOAT, `NEK6` FLOAT, `NEK7` FLOAT, 
        `NEK8` FLOAT, `NEK9` FLOAT, `NIK` FLOAT, `NIM1` FLOAT, `NLK` FLOAT, 
        `NUAK1` FLOAT, `NUAK2` FLOAT, `OSR1` FLOAT, `P38A` FLOAT, `P38B` FLOAT, 
        `P38D` FLOAT, `P38G` FLOAT, `P70S6K` FLOAT, `P70S6KB` FLOAT, `P90RSK` FLOAT, 
        `PAK1` FLOAT, `PAK2` FLOAT, `PAK3` FLOAT, `PAK4` FLOAT, `PAK5` FLOAT, 
        `PAK6` FLOAT, `PASK` FLOAT, `PBK` FLOAT, `PDHK1` FLOAT, `PDHK4` FLOAT, 
        `PDK1` FLOAT, `PERK` FLOAT, `PHKG1` FLOAT, `PHKG2` FLOAT, `PIM1` FLOAT, 
        `PIM2` FLOAT, `PIM3` FLOAT, `PINK1` FLOAT, `PKACA` FLOAT, `PKACB` FLOAT, 
        `PKACG` FLOAT, `PKCA` FLOAT, `PKCB` FLOAT, `PKCD` FLOAT, `PKCE` FLOAT, 
        `PKCG` FLOAT, `PKCH` FLOAT, `PKCI` FLOAT, `PKCT` FLOAT, `PKCZ` FLOAT, 
        `PKG1` FLOAT, `PKG2` FLOAT, `PKN1` FLOAT, `PKN2` FLOAT, `PKN3` FLOAT, 
        `PKR` FLOAT, `PLK1` FLOAT, `PLK2` FLOAT, `PLK3` FLOAT, `PLK4` FLOAT, 
        `PRKD1` FLOAT, `PRKD2` FLOAT, `PRKD3` FLOAT, `PRKX` FLOAT, `PRP4` FLOAT, 
        `PRPK` FLOAT, `QIK` FLOAT, `QSK` FLOAT, `RAF1` FLOAT, `RIPK1` FLOAT, 
        `RIPK2` FLOAT, `RIPK3` FLOAT, `ROCK1` FLOAT, `ROCK2` FLOAT, `RSK2` FLOAT, 
        `RSK3` FLOAT, `RSK4` FLOAT, `SBK` FLOAT, `SGK1` FLOAT, `SGK3` FLOAT, 
        `SIK` FLOAT, `SKMLCK` FLOAT, `SLK` FLOAT, `SMG1` FLOAT, `SMMLCK` FLOAT, 
        `SNRK` FLOAT, `SRPK1` FLOAT, `SRPK2` FLOAT, `SRPK3` FLOAT, `SSTK` FLOAT, 
        `STK33` FLOAT, `STLK3` FLOAT, `TAK1` FLOAT, `TAO1` FLOAT, `TAO2` FLOAT, 
        `TAO3` FLOAT, `TBK1` FLOAT, `TGFBR1` FLOAT, `TGFBR2` FLOAT, `TLK1` FLOAT, 
        `TLK2` FLOAT, `TNIK` FLOAT, `TSSK1` FLOAT, `TSSK2` FLOAT, `TTBK1` FLOAT, 
        `TTBK2` FLOAT, `TTK` FLOAT, `ULK1` FLOAT, `ULK2` FLOAT, `VRK1` FLOAT, 
        `VRK2` FLOAT, `WNK1` FLOAT, `WNK3` FLOAT, `WNK4` FLOAT, `YANK2` FLOAT, 
        `YANK3` FLOAT, `YSK1` FLOAT, `YSK4` FLOAT, `ZAK` FLOAT,
        PRIMARY KEY (`SiteID`),
        INDEX `idx_Cantley_Kinome_Scores_All_STs_siteid` (`SiteID`)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci
    """
    
    cursor.execute(create_table_query)
    connection.commit()  # Use connection directly instead of cursor.connection
    logger.info("S/T kinase scores table created successfully")

def create_y_table(cursor, connection):
    """Create the Y kinase scores table with proper structure and indexes."""
    logger.info("Creating Y kinase scores table...")
    
    create_table_query = """
    CREATE TABLE IF NOT EXISTS `Cantley_Kinome_Scores_All_Ys` (
        `SiteID` VARCHAR(50) NOT NULL,
        `Motif` VARCHAR(50),
        `ABL` FLOAT, `ACK` FLOAT, `ALK` FLOAT, `ARG` FLOAT, `AXL` FLOAT, 
        `BLK` FLOAT, `BMPR2_TYR` FLOAT, `BRK` FLOAT, `BTK` FLOAT, `CSFR` FLOAT, 
        `CSK` FLOAT, `CTK` FLOAT, `DDR1` FLOAT, `DDR2` FLOAT, `EGFR` FLOAT, 
        `EPHA1` FLOAT, `EPHA2` FLOAT, `EPHA3` FLOAT, `EPHA4` FLOAT, `EPHA5` FLOAT, 
        `EPHA6` FLOAT, `EPHA7` FLOAT, `EPHA8` FLOAT, `EPHB1` FLOAT, `EPHB2` FLOAT, 
        `EPHB3` FLOAT, `EPHB4` FLOAT, `ETK` FLOAT, `FAK` FLOAT, `FER` FLOAT, 
        `FES` FLOAT, `FGFR1` FLOAT, `FGFR2` FLOAT, `FGFR3` FLOAT, `FGFR4` FLOAT, 
        `FGR` FLOAT, `FLT3` FLOAT, `FRK` FLOAT, `FYN` FLOAT, `HCK` FLOAT, 
        `HER2` FLOAT, `HER4` FLOAT, `IGF1R` FLOAT, `INSR` FLOAT, `IRR` FLOAT, 
        `ITK` FLOAT, `JAK1` FLOAT, `JAK2` FLOAT, `JAK3` FLOAT, `KIT` FLOAT, 
        `LCK` FLOAT, `LIMK1_TYR` FLOAT, `LIMK2_TYR` FLOAT, `LTK` FLOAT, `LYN` FLOAT, 
        `MER` FLOAT, `MET` FLOAT, `MKK4_TYR` FLOAT, `MKK6_TYR` FLOAT, `MKK7_TYR` FLOAT, 
        `MST1R` FLOAT, `MUSK` FLOAT, `MYT1_TYR` FLOAT, `NEK10_TYR` FLOAT, `PDGFRA` FLOAT, 
        `PDGFRB` FLOAT, `PDHK1_TYR` FLOAT, `PDHK3_TYR` FLOAT, `PDHK4_TYR` FLOAT, 
        `PINK1_TYR` FLOAT, `PYK2` FLOAT, `RET` FLOAT, `ROS` FLOAT, `SRC` FLOAT, 
        `SRMS` FLOAT, `SYK` FLOAT, `TEC` FLOAT, `TESK1_TYR` FLOAT, `TIE2` FLOAT, 
        `TNK1` FLOAT, `TNNI3K_TYR` FLOAT, `TRKA` FLOAT, `TRKB` FLOAT, `TRKC` FLOAT, 
        `TXK` FLOAT, `TYK2` FLOAT, `TYRO3` FLOAT, `VEGFR1` FLOAT, `VEGFR2` FLOAT, 
        `VEGFR3` FLOAT, `WEE1_TYR` FLOAT, `YES` FLOAT, `ZAP70` FLOAT,
        PRIMARY KEY (`SiteID`),
        INDEX `idx_Cantley_Kinome_Scores_All_Ys_siteid` (`SiteID`)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci
    """
    
    cursor.execute(create_table_query)
    connection.commit()  # Use connection directly instead of cursor.connection
    logger.info("Y kinase scores table created successfully")

def bulk_load_st_data(connection, driver_type):
    """Bulk load S/T kinase scores from CSV file into MySQL table."""
    try:
        cursor = connection.cursor()
        
        # Check if table exists
        try:
            cursor.execute("SELECT COUNT(*) FROM `Cantley_Kinome_Scores_All_STs`")
            count = cursor.fetchone()[0]
            logger.info(f"S/T kinase scores table exists with {count} rows")
            
            # Ask user what to do
            choice = input("S/T kinase scores table already exists. Choose an option:\n"
                           "1. Drop and recreate table\n"
                           "2. Truncate table (faster but keeps structure)\n"
                           "3. Append to existing table\n"
                           "4. Exit without changes\n"
                           "Enter choice (1-4): ")
            
            if choice == '1':
                logger.info("Dropping table...")
                cursor.execute("DROP TABLE IF EXISTS `Cantley_Kinome_Scores_All_STs`")
                connection.commit()
                # Recreate table
                create_st_table(cursor, connection)
            elif choice == '2':
                logger.info("Truncating table...")
                cursor.execute("TRUNCATE TABLE `Cantley_Kinome_Scores_All_STs`")
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
            logger.info(f"S/T kinase scores table doesn't exist yet: {e}")
            # Create the table
            create_st_table(cursor, connection)
        
        # Start the bulk load process
        logger.info("Starting S/T kinase scores bulk load process...")
        
        # Disable keys/indexes for faster loading
        logger.info("Disabling keys for faster loading...")
        cursor.execute("ALTER TABLE `Cantley_Kinome_Scores_All_STs` DISABLE KEYS")
        connection.commit()
        
        # Set local_infile to 1
        cursor.execute("SET GLOBAL local_infile = 1")
        
        # Prepare the LOAD DATA INFILE statement with proper path formatting
        csv_path = ST_FILE_PATH.replace('\\', '/')
        
        # Load data, skipping the first column (Unnamed: 0)
        load_query = """
        LOAD DATA LOCAL INFILE '{}'
        INTO TABLE `Cantley_Kinome_Scores_All_STs`
        FIELDS TERMINATED BY ','
        ENCLOSED BY '"'
        LINES TERMINATED BY '\\n'
        IGNORE 1 LINES
        (@dummy, SiteID, Motif, AAK1, ACVR2A, ACVR2B, AKT1, AKT2, AKT3, ALK2, ALK4, ALPHAK3, AMPKA1, 
        AMPKA2, ANKRD3, ASK1, ATM, ATR, AURA, AURB, AURC, BCKDK, BIKE, BMPR1A, BMPR1B, BMPR2, BRAF, 
        BRSK1, BRSK2, BUB1, CAMK1A, CAMK1B, CAMK1D, CAMK1G, CAMK2A, CAMK2B, CAMK2D, CAMK2G, CAMK4, 
        CAMKK1, CAMKK2, CAMLCK, CDC7, CDK1, CDK10, CDK12, CDK13, CDK14, CDK16, CDK17, CDK18, CDK19, 
        CDK2, CDK3, CDK4, CDK5, CDK6, CDK7, CDK8, CDK9, CDKL1, CDKL5, CHAK1, CHAK2, CHK1, CHK2, CK1A, 
        CK1A2, CK1D, CK1E, CK1G1, CK1G2, CK1G3, CK2A1, CK2A2, CLK1, CLK2, CLK3, CLK4, COT, CRIK, DAPK1, 
        DAPK2, DAPK3, DCAMKL1, DCAMKL2, DLK, DMPK1, DNAPK, DRAK1, DSTYK, DYRK1A, DYRK1B, DYRK2, DYRK3, 
        DYRK4, EEF2K, ERK1, ERK2, ERK5, ERK7, FAM20C, GAK, GCK, GCN2, GRK1, GRK2, GRK3, GRK4, GRK5, 
        GRK6, GRK7, GSK3A, GSK3B, HASPIN, HGK, HIPK1, HIPK2, HIPK3, HIPK4, HPK1, HRI, HUNK, ICK, IKKA, 
        IKKB, IKKE, IRAK1, IRAK4, IRE1, IRE2, JNK1, JNK2, JNK3, KHS1, KHS2, KIS, LATS1, LATS2, LKB1, 
        LOK, LRRK2, MAK, MAP3K15, MAPKAPK2, MAPKAPK3, MAPKAPK5, MARK1, MARK2, MARK3, MARK4, MASTL, 
        MEK1, MEK2, MEK5, MEKK1, MEKK2, MEKK3, MEKK6, MELK, MINK, MLK1, MLK2, MLK3, MLK4, MNK1, MNK2, 
        MOK, MOS, MPSK1, MRCKA, MRCKB, MSK1, MSK2, MST1, MST2, MST3, MST4, MTOR, MYLK4, MYO3A, MYO3B, 
        NDR1, NDR2, NEK1, NEK11, NEK2, NEK3, NEK4, NEK5, NEK6, NEK7, NEK8, NEK9, NIK, NIM1, NLK, NUAK1, 
        NUAK2, OSR1, P38A, P38B, P38D, P38G, P70S6K, P70S6KB, P90RSK, PAK1, PAK2, PAK3, PAK4, PAK5, 
        PAK6, PASK, PBK, PDHK1, PDHK4, PDK1, PERK, PHKG1, PHKG2, PIM1, PIM2, PIM3, PINK1, PKACA, PKACB, 
        PKACG, PKCA, PKCB, PKCD, PKCE, PKCG, PKCH, PKCI, PKCT, PKCZ, PKG1, PKG2, PKN1, PKN2, PKN3, PKR, 
        PLK1, PLK2, PLK3, PLK4, PRKD1, PRKD2, PRKD3, PRKX, PRP4, PRPK, QIK, QSK, RAF1, RIPK1, RIPK2, 
        RIPK3, ROCK1, ROCK2, RSK2, RSK3, RSK4, SBK, SGK1, SGK3, SIK, SKMLCK, SLK, SMG1, SMMLCK, SNRK, 
        SRPK1, SRPK2, SRPK3, SSTK, STK33, STLK3, TAK1, TAO1, TAO2, TAO3, TBK1, TGFBR1, TGFBR2, TLK1, 
        TLK2, TNIK, TSSK1, TSSK2, TTBK1, TTBK2, TTK, ULK1, ULK2, VRK1, VRK2, WNK1, WNK3, WNK4, YANK2, 
        YANK3, YSK1, YSK4, ZAK);
        """.format(csv_path)
        
        # Execute the bulk load
        start_time = time.time()
        logger.info("Executing LOAD DATA INFILE for S/T kinase scores...")
        cursor.execute(load_query)
        connection.commit()
        
        # Enable keys/indexes
        logger.info("Re-enabling keys and building indexes...")
        cursor.execute("ALTER TABLE `Cantley_Kinome_Scores_All_STs` ENABLE KEYS")
        connection.commit()
        
        # Get final row count
        cursor.execute("SELECT COUNT(*) FROM `Cantley_Kinome_Scores_All_STs`")
        final_count = cursor.fetchone()[0]
        
        # Log completion
        elapsed_time = time.time() - start_time
        logger.info(f"S/T kinase scores bulk load completed in {elapsed_time:.2f} seconds")
        logger.info(f"Loaded {final_count} rows into Cantley_Kinome_Scores_All_STs")
        
        return True
    except Exception as e:
        logger.error(f"Error during S/T kinase scores bulk load: {e}")
        return False

def bulk_load_y_data(connection, driver_type):
    """Bulk load Y kinase scores from CSV file into MySQL table."""
    try:
        cursor = connection.cursor()
        
        # Check if table exists
        try:
            cursor.execute("SELECT COUNT(*) FROM `Cantley_Kinome_Scores_All_Ys`")
            count = cursor.fetchone()[0]
            logger.info(f"Y kinase scores table exists with {count} rows")
            
            # Ask user what to do
            choice = input("Y kinase scores table already exists. Choose an option:\n"
                           "1. Drop and recreate table\n"
                           "2. Truncate table (faster but keeps structure)\n"
                           "3. Append to existing table\n"
                           "4. Exit without changes\n"
                           "Enter choice (1-4): ")
            
            if choice == '1':
                logger.info("Dropping table...")
                cursor.execute("DROP TABLE IF EXISTS `Cantley_Kinome_Scores_All_Ys`")
                connection.commit()
                # Recreate table
                create_y_table(cursor, connection)
            elif choice == '2':
                logger.info("Truncating table...")
                cursor.execute("TRUNCATE TABLE `Cantley_Kinome_Scores_All_Ys`")
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
            logger.info(f"Y kinase scores table doesn't exist yet: {e}")
            # Create the table
            create_y_table(cursor, connection)
        
        # Start the bulk load process
        logger.info("Starting Y kinase scores bulk load process...")
        
        # Disable keys/indexes for faster loading
        logger.info("Disabling keys for faster loading...")
        cursor.execute("ALTER TABLE `Cantley_Kinome_Scores_All_Ys` DISABLE KEYS")
        connection.commit()
        
        # Set local_infile to 1
        cursor.execute("SET GLOBAL local_infile = 1")
        
        # Prepare the LOAD DATA INFILE statement with proper path formatting
        csv_path = Y_FILE_PATH.replace('\\', '/')
        
        # Load data, skipping the first column (Unnamed: 0)
        load_query = """
        LOAD DATA LOCAL INFILE '{}'
        INTO TABLE `Cantley_Kinome_Scores_All_Ys`
        FIELDS TERMINATED BY ','
        ENCLOSED BY '"'
        LINES TERMINATED BY '\\n'
        IGNORE 1 LINES
        (@dummy, SiteID, Motif, ABL, ACK, ALK, ARG, AXL, BLK, BMPR2_TYR, BRK, BTK, CSFR, CSK, 
        CTK, DDR1, DDR2, EGFR, EPHA1, EPHA2, EPHA3, EPHA4, EPHA5, EPHA6, EPHA7, EPHA8, EPHB1, 
        EPHB2, EPHB3, EPHB4, ETK, FAK, FER, FES, FGFR1, FGFR2, FGFR3, FGFR4, FGR, FLT3, FRK, 
        FYN, HCK, HER2, HER4, IGF1R, INSR, IRR, ITK, JAK1, JAK2, JAK3, KIT, LCK, LIMK1_TYR, 
        LIMK2_TYR, LTK, LYN, MER, MET, MKK4_TYR, MKK6_TYR, MKK7_TYR, MST1R, MUSK, MYT1_TYR, 
        NEK10_TYR, PDGFRA, PDGFRB, PDHK1_TYR, PDHK3_TYR, PDHK4_TYR, PINK1_TYR, PYK2, RET, ROS, 
        SRC, SRMS, SYK, TEC, TESK1_TYR, TIE2, TNK1, TNNI3K_TYR, TRKA, TRKB, TRKC, TXK, TYK2, 
        TYRO3, VEGFR1, VEGFR2, VEGFR3, WEE1_TYR, YES, ZAP70);
        """.format(csv_path)
        
        # Execute the bulk load
        start_time = time.time()
        logger.info("Executing LOAD DATA INFILE for Y kinase scores...")
        cursor.execute(load_query)
        connection.commit()
        
        # Enable keys/indexes
        logger.info("Re-enabling keys and building indexes...")
        cursor.execute("ALTER TABLE `Cantley_Kinome_Scores_All_Ys` ENABLE KEYS")
        connection.commit()
        
        # Get final row count
        cursor.execute("SELECT COUNT(*) FROM `Cantley_Kinome_Scores_All_Ys`")
        final_count = cursor.fetchone()[0]
        
        # Log completion
        elapsed_time = time.time() - start_time
        logger.info(f"Y kinase scores bulk load completed in {elapsed_time:.2f} seconds")
        logger.info(f"Loaded {final_count} rows into Cantley_Kinome_Scores_All_Ys")
        
        return True
    except Exception as e:
        logger.error(f"Error during Y kinase scores bulk load: {e}")
        return False

def try_direct_mysql_command_st():
    """Try using the mysql command-line client directly for S/T kinase scores."""
    try:
        logger.info("Attempting to use mysql command-line client for S/T kinase scores...")
        
        # Create SQL file with commands
        sql_file_path = "load_st_data.sql"
        csv_path = ST_FILE_PATH.replace('\\', '/')
        
        with open(sql_file_path, 'w') as f:
            f.write("""
            -- Create or replace the table
            DROP TABLE IF EXISTS `Cantley_Kinome_Scores_All_STs`;
            
            CREATE TABLE `Cantley_Kinome_Scores_All_STs` (
                `SiteID` VARCHAR(50) NOT NULL,
                `Motif` VARCHAR(50),
                `AAK1` FLOAT, `ACVR2A` FLOAT, `ACVR2B` FLOAT, `AKT1` FLOAT, `AKT2` FLOAT, 
                `AKT3` FLOAT, `ALK2` FLOAT, `ALK4` FLOAT, `ALPHAK3` FLOAT, `AMPKA1` FLOAT, 
                `AMPKA2` FLOAT, `ANKRD3` FLOAT, `ASK1` FLOAT, `ATM` FLOAT, `ATR` FLOAT, 
                `AURA` FLOAT, `AURB` FLOAT, `AURC` FLOAT, `BCKDK` FLOAT, `BIKE` FLOAT, 
                `BMPR1A` FLOAT, `BMPR1B` FLOAT, `BMPR2` FLOAT, `BRAF` FLOAT, `BRSK1` FLOAT, 
                `BRSK2` FLOAT, `BUB1` FLOAT, `CAMK1A` FLOAT, `CAMK1B` FLOAT, `CAMK1D` FLOAT, 
                `CAMK1G` FLOAT, `CAMK2A` FLOAT, `CAMK2B` FLOAT, `CAMK2D` FLOAT, `CAMK2G` FLOAT, 
                `CAMK4` FLOAT, `CAMKK1` FLOAT, `CAMKK2` FLOAT, `CAMLCK` FLOAT, `CDC7` FLOAT, 
                `CDK1` FLOAT, `CDK10` FLOAT, `CDK12` FLOAT, `CDK13` FLOAT, `CDK14` FLOAT, 
                `CDK16` FLOAT, `CDK17` FLOAT, `CDK18` FLOAT, `CDK19` FLOAT, `CDK2` FLOAT, 
                `CDK3` FLOAT, `CDK4` FLOAT, `CDK5` FLOAT, `CDK6` FLOAT, `CDK7` FLOAT, 
                `CDK8` FLOAT, `CDK9` FLOAT, `CDKL1` FLOAT, `CDKL5` FLOAT, `CHAK1` FLOAT, 
                `CHAK2` FLOAT, `CHK1` FLOAT, `CHK2` FLOAT, `CK1A` FLOAT, `CK1A2` FLOAT, 
                `CK1D` FLOAT, `CK1E` FLOAT, `CK1G1` FLOAT, `CK1G2` FLOAT, `CK1G3` FLOAT, 
                `CK2A1` FLOAT, `CK2A2` FLOAT, `CLK1` FLOAT, `CLK2` FLOAT, `CLK3` FLOAT, 
                `CLK4` FLOAT, `COT` FLOAT, `CRIK` FLOAT, `DAPK1` FLOAT, `DAPK2` FLOAT, 
                `DAPK3` FLOAT, `DCAMKL1` FLOAT, `DCAMKL2` FLOAT, `DLK` FLOAT, `DMPK1` FLOAT, 
                `DNAPK` FLOAT, `DRAK1` FLOAT, `DSTYK` FLOAT, `DYRK1A` FLOAT, `DYRK1B` FLOAT, 
                `DYRK2` FLOAT, `DYRK3` FLOAT, `DYRK4` FLOAT, `EEF2K` FLOAT, `ERK1` FLOAT, 
                `ERK2` FLOAT, `ERK5` FLOAT, `ERK7` FLOAT, `FAM20C` FLOAT, `GAK` FLOAT, 
                `GCK` FLOAT, `GCN2` FLOAT, `GRK1` FLOAT, `GRK2` FLOAT, `GRK3` FLOAT, 
                `GRK4` FLOAT, `GRK5` FLOAT, `GRK6` FLOAT, `GRK7` FLOAT, `GSK3A` FLOAT, 
                `GSK3B` FLOAT, `HASPIN` FLOAT, `HGK` FLOAT, `HIPK1` FLOAT, `HIPK2` FLOAT, 
                `HIPK3` FLOAT, `HIPK4` FLOAT, `HPK1` FLOAT, `HRI` FLOAT, `HUNK` FLOAT, 
                `ICK` FLOAT, `IKKA` FLOAT, `IKKB` FLOAT, `IKKE` FLOAT, `IRAK1` FLOAT, 
                `IRAK4` FLOAT, `IRE1` FLOAT, `IRE2` FLOAT, `JNK1` FLOAT, `JNK2` FLOAT, 
                `JNK3` FLOAT, `KHS1` FLOAT, `KHS2` FLOAT, `KIS` FLOAT, `LATS1` FLOAT, 
                `LATS2` FLOAT, `LKB1` FLOAT, `LOK` FLOAT, `LRRK2` FLOAT, `MAK` FLOAT, 
                `MAP3K15` FLOAT, `MAPKAPK2` FLOAT, `MAPKAPK3` FLOAT, `MAPKAPK5` FLOAT, 
                `MARK1` FLOAT, `MARK2` FLOAT, `MARK3` FLOAT, `MARK4` FLOAT, `MASTL` FLOAT, 
                `MEK1` FLOAT, `MEK2` FLOAT, `MEK5` FLOAT, `MEKK1` FLOAT, `MEKK2` FLOAT, 
                `MEKK3` FLOAT, `MEKK6` FLOAT, `MELK` FLOAT, `MINK` FLOAT, `MLK1` FLOAT, 
                `MLK2` FLOAT, `MLK3` FLOAT, `MLK4` FLOAT, `MNK1` FLOAT, `MNK2` FLOAT, 
                `MOK` FLOAT, `MOS` FLOAT, `MPSK1` FLOAT, `MRCKA` FLOAT, `MRCKB` FLOAT, 
                `MSK1` FLOAT, `MSK2` FLOAT, `MST1` FLOAT, `MST2` FLOAT, `MST3` FLOAT, 
                `MST4` FLOAT, `MTOR` FLOAT, `MYLK4` FLOAT, `MYO3A` FLOAT, `MYO3B` FLOAT, 
                `NDR1` FLOAT, `NDR2` FLOAT, `NEK1` FLOAT, `NEK11` FLOAT, `NEK2` FLOAT, 
                `NEK3` FLOAT, `NEK4` FLOAT, `NEK5` FLOAT, `NEK6` FLOAT, `NEK7` FLOAT, 
                `NEK8` FLOAT, `NEK9` FLOAT, `NIK` FLOAT, `NIM1` FLOAT, `NLK` FLOAT, 
                `NUAK1` FLOAT, `NUAK2` FLOAT, `OSR1` FLOAT, `P38A` FLOAT, `P38B` FLOAT, 
                `P38D` FLOAT, `P38G` FLOAT, `P70S6K` FLOAT, `P70S6KB` FLOAT, `P90RSK` FLOAT, 
                `PAK1` FLOAT, `PAK2` FLOAT, `PAK3` FLOAT, `PAK4` FLOAT, `PAK5` FLOAT, 
                `PAK6` FLOAT, `PASK` FLOAT, `PBK` FLOAT, `PDHK1` FLOAT, `PDHK4` FLOAT, 
                `PDK1` FLOAT, `PERK` FLOAT, `PHKG1` FLOAT, `PHKG2` FLOAT, `PIM1` FLOAT, 
                `PIM2` FLOAT, `PIM3` FLOAT, `PINK1` FLOAT, `PKACA` FLOAT, `PKACB` FLOAT, 
                `PKACG` FLOAT, `PKCA` FLOAT, `PKCB` FLOAT, `PKCD` FLOAT, `PKCE` FLOAT, 
                `PKCG` FLOAT, `PKCH` FLOAT, `PKCI` FLOAT, `PKCT` FLOAT, `PKCZ` FLOAT, 
                `PKG1` FLOAT, `PKG2` FLOAT, `PKN1` FLOAT, `PKN2` FLOAT, `PKN3` FLOAT, 
                `PKR` FLOAT, `PLK1` FLOAT, `PLK2` FLOAT, `PLK3` FLOAT, `PLK4` FLOAT, 
                `PRKD1` FLOAT, `PRKD2` FLOAT, `PRKD3` FLOAT, `PRKX` FLOAT, `PRP4` FLOAT, 
                `PRPK` FLOAT, `QIK` FLOAT, `QSK` FLOAT, `RAF1` FLOAT, `RIPK1` FLOAT, 
                `RIPK2` FLOAT, `RIPK3` FLOAT, `ROCK1` FLOAT, `ROCK2` FLOAT, `RSK2` FLOAT, 
                `RSK3` FLOAT, `RSK4` FLOAT, `SBK` FLOAT, `SGK1` FLOAT, `SGK3` FLOAT, 
                `SIK` FLOAT, `SKMLCK` FLOAT, `SLK` FLOAT, `SMG1` FLOAT, `SMMLCK` FLOAT, 
                `SNRK` FLOAT, `SRPK1` FLOAT, `SRPK2` FLOAT, `SRPK3` FLOAT, `SSTK` FLOAT, 
                `STK33` FLOAT, `STLK3` FLOAT, `TAK1` FLOAT, `TAO1` FLOAT, `TAO2` FLOAT, 
                `TAO3` FLOAT, `TBK1` FLOAT, `TGFBR1` FLOAT, `TGFBR2` FLOAT, `TLK1` FLOAT, 
                `TLK2` FLOAT, `TNIK` FLOAT, `TSSK1` FLOAT, `TSSK2` FLOAT, `TTBK1` FLOAT, 
                `TTBK2` FLOAT, `TTK` FLOAT, `ULK1` FLOAT, `ULK2` FLOAT, `VRK1` FLOAT, 
                `VRK2` FLOAT, `WNK1` FLOAT, `WNK3` FLOAT, `WNK4` FLOAT, `YANK2` FLOAT, 
                `YANK3` FLOAT, `YSK1` FLOAT, `YSK4` FLOAT, `ZAK` FLOAT,
                PRIMARY KEY (`SiteID`),
                INDEX `idx_Cantley_Kinome_Scores_All_STs_siteid` (`SiteID`)
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
            
            -- Disable keys for faster loading
            ALTER TABLE `Cantley_Kinome_Scores_All_STs` DISABLE KEYS;
            
            -- Load data, skipping the first column (Unnamed: 0)
            LOAD DATA LOCAL INFILE '{}' 
            INTO TABLE `Cantley_Kinome_Scores_All_STs`
            FIELDS TERMINATED BY ','
            ENCLOSED BY '"'
            LINES TERMINATED BY '\\n'
            IGNORE 1 LINES
            (@dummy, SiteID, Motif, AAK1, ACVR2A, ACVR2B, AKT1, AKT2, AKT3, ALK2, ALK4, ALPHAK3, AMPKA1, 
            AMPKA2, ANKRD3, ASK1, ATM, ATR, AURA, AURB, AURC, BCKDK, BIKE, BMPR1A, BMPR1B, BMPR2, BRAF, 
            BRSK1, BRSK2, BUB1, CAMK1A, CAMK1B, CAMK1D, CAMK1G, CAMK2A, CAMK2B, CAMK2D, CAMK2G, CAMK4, 
            CAMKK1, CAMKK2, CAMLCK, CDC7, CDK1, CDK10, CDK12, CDK13, CDK14, CDK16, CDK17, CDK18, CDK19, 
            CDK2, CDK3, CDK4, CDK5, CDK6, CDK7, CDK8, CDK9, CDKL1, CDKL5, CHAK1, CHAK2, CHK1, CHK2, CK1A, 
            CK1A2, CK1D, CK1E, CK1G1, CK1G2, CK1G3, CK2A1, CK2A2, CLK1, CLK2, CLK3, CLK4, COT, CRIK, DAPK1, 
            DAPK2, DAPK3, DCAMKL1, DCAMKL2, DLK, DMPK1, DNAPK, DRAK1, DSTYK, DYRK1A, DYRK1B, DYRK2, DYRK3, 
            DYRK4, EEF2K, ERK1, ERK2, ERK5, ERK7, FAM20C, GAK, GCK, GCN2, GRK1, GRK2, GRK3, GRK4, GRK5, 
            GRK6, GRK7, GSK3A, GSK3B, HASPIN, HGK, HIPK1, HIPK2, HIPK3, HIPK4, HPK1, HRI, HUNK, ICK, IKKA, 
            IKKB, IKKE, IRAK1, IRAK4, IRE1, IRE2, JNK1, JNK2, JNK3, KHS1, KHS2, KIS, LATS1, LATS2, LKB1, 
            LOK, LRRK2, MAK, MAP3K15, MAPKAPK2, MAPKAPK3, MAPKAPK5, MARK1, MARK2, MARK3, MARK4, MASTL, 
            MEK1, MEK2, MEK5, MEKK1, MEKK2, MEKK3, MEKK6, MELK, MINK, MLK1, MLK2, MLK3, MLK4, MNK1, MNK2, 
            MOK, MOS, MPSK1, MRCKA, MRCKB, MSK1, MSK2, MST1, MST2, MST3, MST4, MTOR, MYLK4, MYO3A, MYO3B, 
            NDR1, NDR2, NEK1, NEK11, NEK2, NEK3, NEK4, NEK5, NEK6, NEK7, NEK8, NEK9, NIK, NIM1, NLK, NUAK1, 
            NUAK2, OSR1, P38A, P38B, P38D, P38G, P70S6K, P70S6KB, P90RSK, PAK1, PAK2, PAK3, PAK4, PAK5, 
            PAK6, PASK, PBK, PDHK1, PDHK4, PDK1, PERK, PHKG1, PHKG2, PIM1, PIM2, PIM3, PINK1, PKACA, PKACB, 
            PKACG, PKCA, PKCB, PKCD, PKCE, PKCG, PKCH, PKCI, PKCT, PKCZ, PKG1, PKG2, PKN1, PKN2, PKN3, PKR, 
            PLK1, PLK2, PLK3, PLK4, PRKD1, PRKD2, PRKD3, PRKX, PRP4, PRPK, QIK, QSK, RAF1, RIPK1, RIPK2, 
            RIPK3, ROCK1, ROCK2, RSK2, RSK3, RSK4, SBK, SGK1, SGK3, SIK, SKMLCK, SLK, SMG1, SMMLCK, SNRK, 
            SRPK1, SRPK2, SRPK3, SSTK, STK33, STLK3, TAK1, TAO1, TAO2, TAO3, TBK1, TGFBR1, TGFBR2, TLK1, 
            TLK2, TNIK, TSSK1, TSSK2, TTBK1, TTBK2, TTK, ULK1, ULK2, VRK1, VRK2, WNK1, WNK3, WNK4, YANK2, 
            YANK3, YSK1, YSK4, ZAK);
            
            -- Re-enable keys
            ALTER TABLE `Cantley_Kinome_Scores_All_STs` ENABLE KEYS;
            
            -- Show row count
            SELECT COUNT(*) FROM `Cantley_Kinome_Scores_All_STs`;
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
            logger.info("MySQL command for S/T kinase scores executed successfully")
            logger.info(f"Output: {result.stdout}")
            return True
        else:
            logger.error(f"MySQL command for S/T kinase scores failed: {result.stderr}")
            return False
    except Exception as e:
        logger.error(f"Error executing MySQL command for S/T kinase scores: {e}")
        return False

def try_direct_mysql_command_y():
    """Try using the mysql command-line client directly for Y kinase scores."""
    try:
        logger.info("Attempting to use mysql command-line client for Y kinase scores...")
        
        # Create SQL file with commands
        sql_file_path = "load_y_data.sql"
        csv_path = Y_FILE_PATH.replace('\\', '/')
        
        with open(sql_file_path, 'w') as f:
            f.write("""
            -- Create or replace the table
            DROP TABLE IF EXISTS `Cantley_Kinome_Scores_All_Ys`;
            
            CREATE TABLE `Cantley_Kinome_Scores_All_Ys` (
                `SiteID` VARCHAR(50) NOT NULL,
                `Motif` VARCHAR(50),
                `ABL` FLOAT, `ACK` FLOAT, `ALK` FLOAT, `ARG` FLOAT, `AXL` FLOAT, 
                `BLK` FLOAT, `BMPR2_TYR` FLOAT, `BRK` FLOAT, `BTK` FLOAT, `CSFR` FLOAT, 
                `CSK` FLOAT, `CTK` FLOAT, `DDR1` FLOAT, `DDR2` FLOAT, `EGFR` FLOAT, 
                `EPHA1` FLOAT, `EPHA2` FLOAT, `EPHA3` FLOAT, `EPHA4` FLOAT, `EPHA5` FLOAT, 
                `EPHA6` FLOAT, `EPHA7` FLOAT, `EPHA8` FLOAT, `EPHB1` FLOAT, `EPHB2` FLOAT, 
                `EPHB3` FLOAT, `EPHB4` FLOAT, `ETK` FLOAT, `FAK` FLOAT, `FER` FLOAT, 
                `FES` FLOAT, `FGFR1` FLOAT, `FGFR2` FLOAT, `FGFR3` FLOAT, `FGFR4` FLOAT, 
                `FGR` FLOAT, `FLT3` FLOAT, `FRK` FLOAT, `FYN` FLOAT, `HCK` FLOAT, 
                `HER2` FLOAT, `HER4` FLOAT, `IGF1R` FLOAT, `INSR` FLOAT, `IRR` FLOAT, 
                `ITK` FLOAT, `JAK1` FLOAT, `JAK2` FLOAT, `JAK3` FLOAT, `KIT` FLOAT, 
                `LCK` FLOAT, `LIMK1_TYR` FLOAT, `LIMK2_TYR` FLOAT, `LTK` FLOAT, `LYN` FLOAT, 
                `MER` FLOAT, `MET` FLOAT, `MKK4_TYR` FLOAT, `MKK6_TYR` FLOAT, `MKK7_TYR` FLOAT, 
                `MST1R` FLOAT, `MUSK` FLOAT, `MYT1_TYR` FLOAT, `NEK10_TYR` FLOAT, `PDGFRA` FLOAT, 
                `PDGFRB` FLOAT, `PDHK1_TYR` FLOAT, `PDHK3_TYR` FLOAT, `PDHK4_TYR` FLOAT, 
                `PINK1_TYR` FLOAT, `PYK2` FLOAT, `RET` FLOAT, `ROS` FLOAT, `SRC` FLOAT, 
                `SRMS` FLOAT, `SYK` FLOAT, `TEC` FLOAT, `TESK1_TYR` FLOAT, `TIE2` FLOAT, 
                `TNK1` FLOAT, `TNNI3K_TYR` FLOAT, `TRKA` FLOAT, `TRKB` FLOAT, `TRKC` FLOAT, 
                `TXK` FLOAT, `TYK2` FLOAT, `TYRO3` FLOAT, `VEGFR1` FLOAT, `VEGFR2` FLOAT, 
                `VEGFR3` FLOAT, `WEE1_TYR` FLOAT, `YES` FLOAT, `ZAP70` FLOAT,
                PRIMARY KEY (`SiteID`),
                INDEX `idx_Cantley_Kinome_Scores_All_Ys_siteid` (`SiteID`)
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
            
            -- Disable keys for faster loading
            ALTER TABLE `Cantley_Kinome_Scores_All_Ys` DISABLE KEYS;
            
            -- Load data, skipping the first column (Unnamed: 0)
            LOAD DATA LOCAL INFILE '{}' 
            INTO TABLE `Cantley_Kinome_Scores_All_Ys`
            FIELDS TERMINATED BY ','
            ENCLOSED BY '"'
            LINES TERMINATED BY '\\n'
            IGNORE 1 LINES
            (@dummy, SiteID, Motif, ABL, ACK, ALK, ARG, AXL, BLK, BMPR2_TYR, BRK, BTK, CSFR, CSK, 
            CTK, DDR1, DDR2, EGFR, EPHA1, EPHA2, EPHA3, EPHA4, EPHA5, EPHA6, EPHA7, EPHA8, EPHB1, 
            EPHB2, EPHB3, EPHB4, ETK, FAK, FER, FES, FGFR1, FGFR2, FGFR3, FGFR4, FGR, FLT3, FRK, 
            FYN, HCK, HER2, HER4, IGF1R, INSR, IRR, ITK, JAK1, JAK2, JAK3, KIT, LCK, LIMK1_TYR, 
            LIMK2_TYR, LTK, LYN, MER, MET, MKK4_TYR, MKK6_TYR, MKK7_TYR, MST1R, MUSK, MYT1_TYR, 
            NEK10_TYR, PDGFRA, PDGFRB, PDHK1_TYR, PDHK3_TYR, PDHK4_TYR, PINK1_TYR, PYK2, RET, ROS, 
            SRC, SRMS, SYK, TEC, TESK1_TYR, TIE2, TNK1, TNNI3K_TYR, TRKA, TRKB, TRKC, TXK, TYK2, 
            TYRO3, VEGFR1, VEGFR2, VEGFR3, WEE1_TYR, YES, ZAP70);
            
            -- Re-enable keys
            ALTER TABLE `Cantley_Kinome_Scores_All_Ys` ENABLE KEYS;
            
            -- Show row count
            SELECT COUNT(*) FROM `Cantley_Kinome_Scores_All_Ys`;
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
        logger.info(f"Executing: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            logger.info("MySQL command for Y kinase scores executed successfully")
            logger.info(f"Output: {result.stdout}")
            return True
        else:
            logger.error(f"MySQL command for Y kinase scores failed: {result.stderr}")
            return False
    except Exception as e:
        logger.error(f"Error executing MySQL command for Y kinase scores: {e}")
        return False
    

def use_pandas_fallback_st():
    """Use pandas as a fallback method to load S/T kinase scores data."""
    try:
        logger.info("Using pandas as fallback method for loading S/T kinase scores...")
        
        # Import SQLAlchemy
        try:
            from sqlalchemy import create_engine
        except ImportError:
            logger.error("SQLAlchemy is required for pandas fallback method. Install it with: pip install sqlalchemy")
            return False
        
        # Create engine
        import urllib.parse
        password = urllib.parse.quote_plus(DB_PASS)
        engine = create_engine(f"mysql+pymysql://{DB_USER}:{password}@{DB_HOST}:{DB_PORT}/{DB_NAME}")
        
        # Load data from CSV and drop the unwanted column
        logger.info(f"Loading S/T kinase scores CSV: {ST_FILE_PATH}")
        df = pd.read_csv(ST_FILE_PATH)
        logger.info(f"CSV columns: {df.columns.tolist()[:10]}...")
        
        # Drop the unwanted column if it exists
        if 'Unnamed: 0' in df.columns:
            logger.info("Removing 'Unnamed: 0' row number column")
            df = df.drop(columns=['Unnamed: 0'])
        
        logger.info(f"Loaded {len(df)} rows from CSV")
        
        # Create table if it doesn't exist
        with engine.connect() as conn:
            create_table_query = """
            CREATE TABLE IF NOT EXISTS `Cantley_Kinome_Scores_All_STs` (
                `SiteID` VARCHAR(50) NOT NULL,
                `Motif` VARCHAR(50),
                -- Many other columns omitted for brevity
                PRIMARY KEY (`SiteID`),
                INDEX `idx_Cantley_Kinome_Scores_All_STs_siteid` (`SiteID`)
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci
            """
            conn.execute(create_table_query)
            
            # Ask user what to do if table exists
            try:
                result = conn.execute("SELECT COUNT(*) FROM `Cantley_Kinome_Scores_All_STs`")
                count = result.fetchone()[0]
                if count > 0:
                    choice = input("S/T kinase scores table already has data. Choose an option:\n"
                                "1. Drop and recreate table\n"
                                "2. Truncate table (faster but keeps structure)\n"
                                "3. Append to existing table\n"
                                "4. Exit without changes\n"
                                "Enter choice (1-4): ")
                    
                    if choice == '1':
                        logger.info("Dropping table...")
                        conn.execute("DROP TABLE IF EXISTS `Cantley_Kinome_Scores_All_STs`")
                        # Table will be recreated by to_sql
                    elif choice == '2':
                        logger.info("Truncating table...")
                        conn.execute("TRUNCATE TABLE `Cantley_Kinome_Scores_All_STs`")
                    elif choice == '3':
                        logger.info("Will append data to existing table")
                    elif choice == '4':
                        logger.info("Exiting without changes")
                        return False
                    else:
                        logger.warning("Invalid choice. Exiting.")
                        return False
            except:
                # Table doesn't exist yet
                pass
        
        # Write to database in chunks
        start_time = time.time()
        logger.info("Writing S/T kinase scores data to MySQL...")
        
        chunksize = 1000
        for i in range(0, len(df), chunksize):
            logger.info(f"Writing chunk {i // chunksize + 1}/{(len(df) + chunksize - 1) // chunksize}...")
            chunk = df.iloc[i:i+chunksize]
            
            if i == 0:
                # First chunk replaces or creates the table
                if_exists = 'replace'
            else:
                # Subsequent chunks append to the table
                if_exists = 'append'
            
            chunk.to_sql(
                'Cantley_Kinome_Scores_All_STs',
                engine,
                if_exists=if_exists,
                index=False,  # Crucial to prevent adding an index column
                chunksize=100  # Smaller sub-chunks for each SQLAlchemy transaction
            )
        
        # Create index if using replace
        with engine.connect() as conn:
            conn.execute("ALTER TABLE `Cantley_Kinome_Scores_All_STs` ADD PRIMARY KEY (`SiteID`)")
            conn.execute("CREATE INDEX `idx_Cantley_Kinome_Scores_All_STs_siteid` ON `Cantley_Kinome_Scores_All_STs` (`SiteID`)")
            
            # Get final count
            result = conn.execute("SELECT COUNT(*) FROM `Cantley_Kinome_Scores_All_STs`")
            count = result.fetchone()[0]
        
        elapsed_time = time.time() - start_time
        logger.info(f"S/T kinase scores data loaded successfully using pandas in {elapsed_time:.2f} seconds")
        logger.info(f"Loaded {count} rows")
        
        return True
    except Exception as e:
        logger.error(f"Error in pandas fallback for S/T kinase scores: {e}")
        return False

def use_pandas_fallback_y():
    """Use pandas as a fallback method to load Y kinase scores data."""
    try:
        logger.info("Using pandas as fallback method for loading Y kinase scores...")
        
        # Import SQLAlchemy
        try:
            from sqlalchemy import create_engine
        except ImportError:
            logger.error("SQLAlchemy is required for pandas fallback method. Install it with: pip install sqlalchemy")
            return False
        
        # Create engine
        import urllib.parse
        password = urllib.parse.quote_plus(DB_PASS)
        engine = create_engine(f"mysql+pymysql://{DB_USER}:{password}@{DB_HOST}:{DB_PORT}/{DB_NAME}")
        
        # Load data from CSV and drop the unwanted column
        logger.info(f"Loading Y kinase scores CSV: {Y_FILE_PATH}")
        df = pd.read_csv(Y_FILE_PATH)
        logger.info(f"CSV columns: {df.columns.tolist()[:10]}...")
        
        # Drop the unwanted column if it exists
        if 'Unnamed: 0' in df.columns:
            logger.info("Removing 'Unnamed: 0' row number column")
            df = df.drop(columns=['Unnamed: 0'])
        
        logger.info(f"Loaded {len(df)} rows from CSV")
        
        # Create table if it doesn't exist
        with engine.connect() as conn:
            create_table_query = """
            CREATE TABLE IF NOT EXISTS `Cantley_Kinome_Scores_All_Ys` (
                `SiteID` VARCHAR(50) NOT NULL,
                `Motif` VARCHAR(50),
                -- Many other columns omitted for brevity
                PRIMARY KEY (`SiteID`),
                INDEX `idx_Cantley_Kinome_Scores_All_Ys_siteid` (`SiteID`)
            ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci
            """
            conn.execute(create_table_query)
            
            # Ask user what to do if table exists
            try:
                result = conn.execute("SELECT COUNT(*) FROM `Cantley_Kinome_Scores_All_Ys`")
                count = result.fetchone()[0]
                if count > 0:
                    choice = input("Y kinase scores table already has data. Choose an option:\n"
                                "1. Drop and recreate table\n"
                                "2. Truncate table (faster but keeps structure)\n"
                                "3. Append to existing table\n"
                                "4. Exit without changes\n"
                                "Enter choice (1-4): ")
                    
                    if choice == '1':
                        logger.info("Dropping table...")
                        conn.execute("DROP TABLE IF EXISTS `Cantley_Kinome_Scores_All_Ys`")
                        # Table will be recreated by to_sql
                    elif choice == '2':
                        logger.info("Truncating table...")
                        conn.execute("TRUNCATE TABLE `Cantley_Kinome_Scores_All_Ys`")
                    elif choice == '3':
                        logger.info("Will append data to existing table")
                    elif choice == '4':
                        logger.info("Exiting without changes")
                        return False
                    else:
                        logger.warning("Invalid choice. Exiting.")
                        return False
            except:
                # Table doesn't exist yet
                pass
        
        # Write to database in chunks
        start_time = time.time()
        logger.info("Writing Y kinase scores data to MySQL...")
        
        chunksize = 1000
        for i in range(0, len(df), chunksize):
            logger.info(f"Writing chunk {i // chunksize + 1}/{(len(df) + chunksize - 1) // chunksize}...")
            chunk = df.iloc[i:i+chunksize]
            
            if i == 0:
                # First chunk replaces or creates the table
                if_exists = 'replace'
            else:
                # Subsequent chunks append to the table
                if_exists = 'append'
            
            chunk.to_sql(
                'Cantley_Kinome_Scores_All_Ys',
                engine,
                if_exists=if_exists,
                index=False,  # Crucial to prevent adding an index column
                chunksize=100  # Smaller sub-chunks for each SQLAlchemy transaction
            )
        
        # Create index if using replace
        with engine.connect() as conn:
            conn.execute("ALTER TABLE `Cantley_Kinome_Scores_All_Ys` ADD PRIMARY KEY (`SiteID`)")
            conn.execute("CREATE INDEX `idx_Cantley_Kinome_Scores_All_Ys_siteid` ON `Cantley_Kinome_Scores_All_Ys` (`SiteID`)")
            
            # Get final count
            result = conn.execute("SELECT COUNT(*) FROM `Cantley_Kinome_Scores_All_Ys`")
            count = result.fetchone()[0]
        
        elapsed_time = time.time() - start_time
        logger.info(f"Y kinase scores data loaded successfully using pandas in {elapsed_time:.2f} seconds")
        logger.info(f"Loaded {count} rows")
        
        return True
    except Exception as e:
        logger.error(f"Error in pandas fallback for Y kinase scores: {e}")
        return False

def verify_data_load(connection):
    """Verify that data was loaded correctly by reading samples from the tables."""
    try:
        cursor = connection.cursor()
        
        # Check S/T kinase scores table
        try:
            # Get column names
            cursor.execute("SHOW COLUMNS FROM `Cantley_Kinome_Scores_All_STs`")
            st_columns = [row[0] for row in cursor.fetchall()]
            logger.info(f"S/T kinase scores table has {len(st_columns)} columns")
            
            # Check if the 'Unnamed: 0' column is present (it shouldn't be)
            if 'Unnamed: 0' in st_columns:
                logger.warning("WARNING: 'Unnamed: 0' column found in S/T table!")
            
            # Get row count
            cursor.execute("SELECT COUNT(*) FROM `Cantley_Kinome_Scores_All_STs`")
            st_count = cursor.fetchone()[0]
            logger.info(f"S/T kinase scores table has {st_count} rows")
            
            # Read a few sample rows
            cursor.execute("SELECT SiteID, Motif, AKT1, CDK1, ERK1 FROM `Cantley_Kinome_Scores_All_STs` LIMIT 5")
            st_samples = cursor.fetchall()
            logger.info("Sample S/T kinase scores data:")
            for sample in st_samples:
                logger.info(f"  SiteID: {sample[0]}, Motif: {sample[1]}, AKT1: {sample[2]}, CDK1: {sample[3]}, ERK1: {sample[4]}")
            
            # Get top kinase for a random site
            cursor.execute("""
            SELECT SiteID, Motif, 
                   GREATEST(AAK1, ACVR2A, AKT1, CDK1, ERK1, GSK3A, GSK3B, P38A, PKACA) AS max_score,
                   CASE 
                       WHEN AAK1 = GREATEST(AAK1, ACVR2A, AKT1, CDK1, ERK1, GSK3A, GSK3B, P38A, PKACA) THEN 'AAK1'
                       WHEN ACVR2A = GREATEST(AAK1, ACVR2A, AKT1, CDK1, ERK1, GSK3A, GSK3B, P38A, PKACA) THEN 'ACVR2A'
                       WHEN AKT1 = GREATEST(AAK1, ACVR2A, AKT1, CDK1, ERK1, GSK3A, GSK3B, P38A, PKACA) THEN 'AKT1'
                       WHEN CDK1 = GREATEST(AAK1, ACVR2A, AKT1, CDK1, ERK1, GSK3A, GSK3B, P38A, PKACA) THEN 'CDK1'
                       WHEN ERK1 = GREATEST(AAK1, ACVR2A, AKT1, CDK1, ERK1, GSK3A, GSK3B, P38A, PKACA) THEN 'ERK1'
                       WHEN GSK3A = GREATEST(AAK1, ACVR2A, AKT1, CDK1, ERK1, GSK3A, GSK3B, P38A, PKACA) THEN 'GSK3A'
                       WHEN GSK3B = GREATEST(AAK1, ACVR2A, AKT1, CDK1, ERK1, GSK3A, GSK3B, P38A, PKACA) THEN 'GSK3B'
                       WHEN P38A = GREATEST(AAK1, ACVR2A, AKT1, CDK1, ERK1, GSK3A, GSK3B, P38A, PKACA) THEN 'P38A'
                       WHEN PKACA = GREATEST(AAK1, ACVR2A, AKT1, CDK1, ERK1, GSK3A, GSK3B, P38A, PKACA) THEN 'PKACA'
                   END AS top_kinase
            FROM `Cantley_Kinome_Scores_All_STs`
            ORDER BY RAND()
            LIMIT 3
            """)
            st_top_kinases = cursor.fetchall()
            logger.info("Top kinases for random S/T sites:")
            for site in st_top_kinases:
                logger.info(f"  SiteID: {site[0]}, Motif: {site[1]}, Top Kinase: {site[3]} (Score: {site[2]})")
                
        except Exception as e:
            logger.warning(f"Error verifying S/T kinase scores table: {e}")
        
        # Check Y kinase scores table
        try:
            # Get column names
            cursor.execute("SHOW COLUMNS FROM `Cantley_Kinome_Scores_All_Ys`")
            y_columns = [row[0] for row in cursor.fetchall()]
            logger.info(f"Y kinase scores table has {len(y_columns)} columns")
            
            # Check if the 'Unnamed: 0' column is present (it shouldn't be)
            if 'Unnamed: 0' in y_columns:
                logger.warning("WARNING: 'Unnamed: 0' column found in Y table!")
            
            # Get row count
            cursor.execute("SELECT COUNT(*) FROM `Cantley_Kinome_Scores_All_Ys`")
            y_count = cursor.fetchone()[0]
            logger.info(f"Y kinase scores table has {y_count} rows")
            
            # Read a few sample rows
            cursor.execute("SELECT SiteID, Motif, SRC, EGFR, JAK1 FROM `Cantley_Kinome_Scores_All_Ys` LIMIT 5")
            y_samples = cursor.fetchall()
            logger.info("Sample Y kinase scores data:")
            for sample in y_samples:
                logger.info(f"  SiteID: {sample[0]}, Motif: {sample[1]}, SRC: {sample[2]}, EGFR: {sample[3]}, JAK1: {sample[4]}")
            
            # Get top kinase for a random site
            cursor.execute("""
            SELECT SiteID, Motif, 
                   GREATEST(SRC, EGFR, JAK1, ABL, LCK) AS max_score,
                   CASE 
                       WHEN SRC = GREATEST(SRC, EGFR, JAK1, ABL, LCK) THEN 'SRC'
                       WHEN EGFR = GREATEST(SRC, EGFR, JAK1, ABL, LCK) THEN 'EGFR'
                       WHEN JAK1 = GREATEST(SRC, EGFR, JAK1, ABL, LCK) THEN 'JAK1'
                       WHEN ABL = GREATEST(SRC, EGFR, JAK1, ABL, LCK) THEN 'ABL'
                       WHEN LCK = GREATEST(SRC, EGFR, JAK1, ABL, LCK) THEN 'LCK'
                   END AS top_kinase
            FROM `Cantley_Kinome_Scores_All_Ys`
            ORDER BY RAND()
            LIMIT 3
            """)
            y_top_kinases = cursor.fetchall()
            logger.info("Top kinases for random Y sites:")
            for site in y_top_kinases:
                logger.info(f"  SiteID: {site[0]}, Motif: {site[1]}, Top Kinase: {site[3]} (Score: {site[2]})")
                
        except Exception as e:
            logger.warning(f"Error verifying Y kinase scores table: {e}")
            
    except Exception as e:
        logger.error(f"Error during data verification: {e}")

def main():
    """Main function to execute the bulk load process."""
    logger.info("Starting bulk loading process for Cantley kinome scores")
    
    # Check if input files exist
    if not os.path.exists(ST_FILE_PATH):
        logger.error(f"S/T kinase scores file not found: {ST_FILE_PATH}")
        return False
    else:
        logger.info(f"Found S/T kinase scores file: {ST_FILE_PATH}")
    
    if not os.path.exists(Y_FILE_PATH):
        logger.error(f"Y kinase scores file not found: {Y_FILE_PATH}")
        return False
    else:
        logger.info(f"Found Y kinase scores file: {Y_FILE_PATH}")
    
    # Try to use Python's MySQL libraries first
    connection, driver_type = try_connect_mysql()
    
    st_loaded = False
    y_loaded = False
    
    if connection:
        # Bulk load S/T kinase scores using Python
        logger.info("Loading S/T kinase scores...")
        st_loaded = bulk_load_st_data(connection, driver_type)
        
        if st_loaded:
            logger.info("S/T kinase scores loaded successfully using Python")
        else:
            logger.warning("Failed to load S/T kinase scores using Python")
        
        # Bulk load Y kinase scores using Python
        logger.info("Loading Y kinase scores...")
        y_loaded = bulk_load_y_data(connection, driver_type)
        
        if y_loaded:
            logger.info("Y kinase scores loaded successfully using Python")
        else:
            logger.warning("Failed to load Y kinase scores using Python")
        
        # Add verification step - read some data from the tables
        if st_loaded or y_loaded:
            logger.info("Verifying data in tables...")
            verify_data_load(connection)
        
        connection.close()
    else:
        logger.warning("Could not connect to MySQL with Python drivers, trying direct mysql command...")
    
    # If Python approach failed for S/T kinase scores, try direct mysql command
    if not st_loaded:
        logger.info("Trying direct mysql command for S/T kinase scores...")
        st_loaded = try_direct_mysql_command_st()
        
        if st_loaded:
            logger.info("S/T kinase scores loaded successfully using mysql command")
        else:
            logger.warning("Failed to load S/T kinase scores using mysql command")
    
    # If Python approach failed for Y kinase scores, try direct mysql command
    if not y_loaded:
        logger.info("Trying direct mysql command for Y kinase scores...")
        y_loaded = try_direct_mysql_command_y()
        
        if y_loaded:
            logger.info("Y kinase scores loaded successfully using mysql command")
        else:
            logger.warning("Failed to load Y kinase scores using mysql command")
    
    # If both approaches failed for S/T kinase scores, try pandas fallback
    if not st_loaded:
        logger.info("Trying pandas fallback for S/T kinase scores...")
        st_loaded = use_pandas_fallback_st()
        
        if st_loaded:
            logger.info("S/T kinase scores loaded successfully using pandas fallback")
        else:
            logger.error("All methods failed to load S/T kinase scores")
    
    # If both approaches failed for Y kinase scores, try pandas fallback
    if not y_loaded:
        logger.info("Trying pandas fallback for Y kinase scores...")
        y_loaded = use_pandas_fallback_y()
        
        if y_loaded:
            logger.info("Y kinase scores loaded successfully using pandas fallback")
        else:
            logger.error("All methods failed to load Y kinase scores")
    
    # Final verification after all methods have been tried
    if st_loaded or y_loaded:
        logger.info("Performing final data verification...")
        final_conn, final_driver_type = try_connect_mysql()
        if final_conn:
            verify_data_load(final_conn)
            final_conn.close()
    
    return st_loaded and y_loaded

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