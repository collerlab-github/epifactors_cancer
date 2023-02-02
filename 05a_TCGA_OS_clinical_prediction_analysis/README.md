# Prognostic Gene Survival Analysis for Overall Survival in TCGA Cancers

Prerequisites: 03_nmf and its prerequisites must be completed.

This folder contains scripts for single-gene stratification and survival analysis of patients for the OS endpoints with different covariates for TCGA cancer types. Here are the descriptions for the scripts.


## 01_survminer.R

This script conducts uses survminer to stratify patients into high and low expression of each epigene and then performs survival analysis with cox regression to correct for patient age and gender.

