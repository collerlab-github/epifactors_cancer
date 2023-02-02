# Mutation and CNV Analysis for Clinically Prognostic Genes

Prerequisites: 05_single_gene_clinical_analysis and its prerequisites must be completed.

This folder contains scripts for gathering and analyzing mutation and CNV data from the CBioPortal. Here are the descriptions for the scripts.

## Collecting Cbioportal Mutation and CNA data

Copy patient ids from 'survival_patients.txt' from 05_gene_clinical_prediction_analysis
Go to cbioportal.org/study/summary?id=study_id (ex https://www.cbioportal.org/study/summary?id=prad_tcga_pub)
Under 'Custom Selection' select 'By Patient ID' and paste patient ids.

Download mutation and CNA data from respective tables.

## mutation_and_cna_gene.R

This script formats mutation and cnv data from CBioPortal.

## cnv_fisher_heatmap.R

This script calculates cnv enrichment in the high or low groups for each prognostic gene in each cancer type and visualizes results with a heatmap.

## mutation_fisher_heatmap.R

This script calculates mutation enrichment in the high or low groups for each prognostic gene in each cancer type and visualizes results with a heatmap.
