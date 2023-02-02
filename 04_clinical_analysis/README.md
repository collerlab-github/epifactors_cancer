# Survival Analysis based on NMF clusters

Prerequisites: 03_nmf and its prerequisites must be completed.

This folder contains scripts for running survival analysis for the OS, DSS, and PFI endpoints with different covariates. Here are the descriptions for the scripts.


## 01_survival.R

This script conducts survival analysis with cox regression to correct for patient age and gender.

## 01_survival_purity.R

This script conducts survival analysis with cox regression to correct for patient age, gender, and purity information.

## 01_survival_tp53.R

This script conducts survival analysis with cox regression to correct for patient age, gender, and p53 expression.
