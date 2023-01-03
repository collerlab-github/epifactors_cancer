#!/bin/bash

if [ "$#" -ne 2 ];
then
	echo "Usage: bash master_script.sh [cancer] [purity]"
	exit 1
fi

cancer=$1
purity=$2

# 01 data collection
echo 'Data Collection'
cd 01_data_collection
echo 'Extracting purity manifest'
Rscript 01_extract_purity_manifest.R $cancer $purity
echo 'Gathering raw counts'
bash 02_unzip.sh $cancer
echo 'Consolidating patients'
Rscript 03_raw_counts_dataframe.R $cancer
cd ../

# 02 filtering
echo 'Data Filtering'
cd 02_filtering
echo 'Normalizing Data'
Rscript 01_normalize_counts.R $cancer
echo 'Filtering for epigenes'
Rscript 02_epigenes_filter.R $cancer

echo 'Standard Deviation Cutoffs written'
cd ../

echo 'finished'