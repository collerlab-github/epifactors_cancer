# Prognostic Gene Survival Analysis

Prerequisites: 03_nmf and its prerequisites must be completed.

This folder contains scripts for single-gene stratification and survival analysis of patients for the OS and PFI endpoints with different covariates.


Here are the descriptions for the scripts.


## 01_survminer.R

This script conducts uses survminer to stratify patients into high and low expression of each epigene and then performs survival analysis with cox regression to correct for patient age and gender.

## 01_survminer_one_gene.R

This script conducts uses survminer to stratify patients into high and low expression of a user-specified epigene and then performs survival analysis with cox regression to correct for patient age and gender.

## HR_forestplot.R

This script plots the hazard ratios of the most significant endpoint-prognostic gene for each cancer.

## circos.R

This script plots a circos plot for the most frequent endpoint-prognostic genes among all cancers.
