# Prognostic Gene Survival Analysis with proliferation-dependent gene correction

Prerequisites: 05_nmf and its prerequisites must be completed.

This folder contains scripts for single-gene stratification and survival analysis of patients for the clinical endpoints. Here are the descriptions for the scripts.


## 01_survminer_metaPCNA.R

This script calculates the median log2 normalize expression of proliferation-dependent genes as a metaPCNA proliferation-dependent score.
Then it conducts survminer to stratify patients into high and low expression of each epigene and performs survival analysis with cox regression to correct for patient age and gender and metaPCNA score.

