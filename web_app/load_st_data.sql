
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
            
            -- Load data
            LOAD DATA LOCAL INFILE 'F:/Kinome/motif_scores_percentile_STs.csv' 
            INTO TABLE `Cantley_Kinome_Scores_All_STs`
            FIELDS TERMINATED BY ','
            ENCLOSED BY '"'
            LINES TERMINATED BY '\n'
            IGNORE 1 LINES;
            
            -- Re-enable keys
            ALTER TABLE `Cantley_Kinome_Scores_All_STs` ENABLE KEYS;
            
            -- Show row count
            SELECT COUNT(*) FROM `Cantley_Kinome_Scores_All_STs`;
            