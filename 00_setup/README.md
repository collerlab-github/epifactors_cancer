# Data Collection

This folder contains scripts for collecting data for epigene analysis.

## Purity Patient Cutoffs

For each cancer type, patients lists were collected from their original TCGA manuscript. Tumor purity values patients collected from Aran et al, 2015. Cutoffs were decided based on the number of patients who exceed the purity threshold in either CPE or CHAT purity score. The one exception is PAAD (PDAC) in which we used the 0.33 ABSOLUTE purity cutoff from the TCGA manuscript.
Cancers in the study and cutoff used

BRCA	0.7
THCA	0.7
OV	0.7
LGG	0.7
PRAD	0.7
SKCM	0.7
UCEC	0.7
KIRC	0.7
CRC	0.7
CESC	0.7
LIHC	0.7
SARC	0.7
HNSC	0.7
KIRP	0.7
GBM	0.7
PCPG	0.7
LUAD	0.7
STAD	0.7
LUSC	0.7
ESCA	0.6
TGCT	0.6
ACC	0.7
BLCA	0.7
PAAD	0.33

 
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
- 07_prognosis_ML_classifier
- papers
- supplemental


# 02_purity_collection.R

This script outputs a table of cutoffs and patients for a user-defined cancer type. Before running, the patient list must be downloaded to the 'supplemental' folder from the original TCGA manuscript.
After deciding a cutoff, the script is run again to output a data frame of patients and purity scores.


