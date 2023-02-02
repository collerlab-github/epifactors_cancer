# Prognostic Gene Survival Analysis

Prerequisites: 03_nmf and its prerequisites must be completed.

This folder contains scripts for single-gene stratification and survival analysis of patients for the OS and PFI endpoints with different covariates.


Here are the descriptions for the scripts.


## 01_survminer.R

This script conducts uses survminer to stratify patients into high and low expression of each epigene and then performs survival analysis with cox regression to correct for patient age and gender.

## 01_survminer_one_gene.R

This script conducts uses survminer to stratify patients into high and low expression of a user-specified epigene and then performs survival analysis with cox regression to correct for patient age and gender.

## 01_survival_tp53.R

This script conducts survival analysis with cox regression to correct for patient age, gender, and p53 expression.
