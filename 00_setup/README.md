# Setup
This folder contains scripts for collecting data for epigene analysis.

Here are the descriptions for the scripts.

## 01_create_directories.sh
This script creates all the folders necessary for the analysis.
- 01_data_collection
- 02_filtering
- 03_nmf
- 04_clinical_analysis
- 05_gene_clinical_prediction_analysis
- 06_diff_gene_info
- papers
- supplemental


## 02_purity_collection.R

This script outputs a table of cutoffs and patients for a user-defined cancer type. Before running, the patient list must be downloaded to the 'supplemental' folder from the original TCGA or TARGET manuscript.
After deciding a cutoff, the script is run again to output a data frame of patients and purity scores.


