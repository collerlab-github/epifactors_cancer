# Mutation and CNV Analysis for Clinically Prognostic Genes

Prerequisites: 05_single_gene_clinical_analysis and its prerequisites must be completed.

This folder contains scripts for gathering and analyzing mutation and CNV data from the CBioPortal. Here are the descriptions for the scripts.

## Collecting Cbioportal Mutation and CNA data

Copy patient ids from 'survival_patients.txt' from 05_gene_clinical_prediction_analysis
Go to cbioportal.org/study/summary?id=study_id (ex https://www.cbioportal.org/study/summary?id=prad_tcga_pub)
Under 'Custom Selection' select 'By Patient ID' and paste patient ids.

Download mutation and CNA data from respective tables.

## 01_survminer.R

This script conducts uses survminer to stratify patients into high and low expression of each epigene and then performs survival analysis with cox regression to correct for patient age and gender.

## 02_survminer_one_gene.R

This script conducts uses survminer to stratify patients into high and low expression of a user-specified epigene and then performs survival analysis with cox regression to correct for patient age and gender.

## HR_forest_plot.R

This script plots the hazard ratios of the most significant endpoint-prognostic gene for each cancer.

## circos_plot.R

This script plots a circos plot for the most frequent endpoint-prognostic genes among all cancers.
