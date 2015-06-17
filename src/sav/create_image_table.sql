LOAD DATA INFILE '/tmp/WGS_ALL_6_15_2015.csv' INTO TABLE wgs_img
COLUMNS TERMINATED BY ','
OPTIONALLY ENCLOSED BY '"'
ESCAPED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 LINES;-- extract human gene map by the given list of gene names

-- delete from wgs_img where img_sn <> -1
