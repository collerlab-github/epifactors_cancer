#!/bin/bash

if [ "$#" -ne 2 ];
then
echo "Usage: bash create_directories.sh [TCGA/TARGET] [cancer]"
exit 1
fi

project=$1
cancer=$2
cd ../../
mkdir ${project}_${cancer}
cd ${project}_${cancer}
mkdir 01_data_collection
mkdir 02_filtering
mkdir 03_nmf
mkdir 04_clinical_analysis
mkdir 05_gene_clinical_prediction_analysis
mkdir 06_diff_gene_info
mkdir papers
mkdir supplemental
cp ../scripts/00_setup/02_purity_collection.R ./01_data_collection/
