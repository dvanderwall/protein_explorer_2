
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
            
            -- Load data
            LOAD DATA LOCAL INFILE 'F:/Kinome/motif_scores_percentile_Ys.csv' 
            INTO TABLE `Cantley_Kinome_Scores_All_Ys`
            FIELDS TERMINATED BY ','
            ENCLOSED BY '"'
            LINES TERMINATED BY '\n'
            IGNORE 1 LINES;
            
            -- Re-enable keys
            ALTER TABLE `Cantley_Kinome_Scores_All_Ys` ENABLE KEYS;
            
            -- Show row count
            SELECT COUNT(*) FROM `Cantley_Kinome_Scores_All_Ys`;
            