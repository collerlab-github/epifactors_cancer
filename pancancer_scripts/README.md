# Pancancer scripts

This folder contains scripts for aggregating and analyzing all TCGA cancers as a pancancer study. Here are the following folders and scripts.

## 01_filtering
### 01_raw_counts_dataframe.R
This script combines all TCGA cancers' raw RNA-seq data into one dataframe.
### 02_normalize_data.R
This script log2 normalizes the pancancer dataset and subsets epigenes.

## 02_UMAP
### 01_umap
This script generates a umap from the pancancer dataset when using all epigenes or only top NMF genes.
