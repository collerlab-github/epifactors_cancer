#!/bin/bash

# NOTE: raw data for the cancer type must be downloaded and stored in the same directory.
# 		Download raw data directly from GDC portal repository

if [ "$#" -ne 2 ];
then
	echo "Usage: bash 02_unzip.sh [project] [cancer]"
	exit 1
fi

project=$1
cancer=$2

cd ../../${project}_${cancer}/01_data_collection
manifest_file=purity_gdc_manifest_$cancer.txt


# run gdc-client to download raw TCGA data from manifest file 
if [ -d raw_data/ ];
then
	rm -r raw_data/
fi

mkdir raw_data/
cd raw_data/
~/Desktop/gdc-client download -m ../${manifest_file}

# unzip raw data files
for file in */
do
	gunzip ${file}*.gz
    mv ${file}*.htseq* ./
    rm -r $file
done
