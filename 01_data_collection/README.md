# Data Collection

This folder contains scripts for collecting data for epigene analysis.

## Purity Patient Cutoffs

For each cancer type, patients lists were collected from their original TCGA manuscript. Tumor purity values patients collected from Aran et al, 2015. Cutoffs were decided based on the number of patients who exceed the purity threshold in either CPE or CHAT purity score. The one exception is PAAD (PDAC) in which we used the 0.33 ABSOLUTE purity cutoff from the TCGA manuscript.
Cancers in the study and cutoff used

| Syntax      | Description |
| ----------- | ----------- |
| Header      | Title       |
| Paragraph   | Text        |
| BRCA	| 0.7 |
| THCA |	0.7 |
| OV	| 0.7 |
| LGG	| 0.7 |
| PRAD	| 0.7 |
| SKCM	| 0.7 |
| UCEC	| 0.7 |
| KIRC	| 0.7 |
| CRC	| 0.7 |
| CESC	| 0.7 |
| LIHC	| 0.7 |
| SARC	| 0.7 |
| HNSC	| 0.7 |
| KIRP	| 0.7 |
| GBM	| 0.7 |
| PCPG	| 0.7 |
| LUAD	| 0.7 |
| STAD	| 0.7 |
| LUSC	| 0.7 |
| ESCA	| 0.6 |
| TGCT	| 0.6 |
| ACC	| 0.7 |
| BLCA	| 0.7 |
| PAAD	| 0.33 |

 
 
## GDC Portal File Download
Prior to running the first script, the following files must be in the directory:
- gdc portal manifest file
- gdc portal samplesheet file

The cancer patients' RNAseq data must be obtained from GDC portal
1. Go to https://portal.gdc.cancer.gov/
2. Click Repository
3. Under 'Cases' select primary site and project name (ex: breast, TCGA-BRCA)
4. Under 'Files' select experimental strategy and workflow type (ex: RNA-seq, HTSeq - Counts)
5. Add all files to cart (make sure cart is initially empty)
6. Go to cart and download sample sheet and download manifest (under Download tab). Save in data collection directory.


Here are the descriptions for the scripts.

## 01_extract_purity_manifest.R

This file extracts patients of a certain purity cutoff from the list of patients from the manifest and sample sheet. The output is a purity-filtered manifest and sample sheet.


## 02_unzip.sh

This script uses the GDC portal command line interface to download the RNA-seq raw counts files of the purity patients. The counts files are saved in the raw_data folder



## 03_raw_counts_dataframe.R

This script combines all patients' raw counts files into one dataframe. The output is saved into the RNA-seq folder (up one directory)

