# NMF

Prerequisites: 01_data_collection and 02_filtering must be completed.

This folder contains scripts for running NMF and generating initial visualizations. The overall process includes:
1. Running NMF, picking an optimal rank, and retrieving patient cluster memberships and top NMF genes.
2. Generating heatmaps for top NMF genes and patient groups.
3. Conduct PCA on all epigenes and compare NMF clusters in PC space.


Here are the descriptions for the scripts.


## 01_NMF

This script runs NMF on the patient x epigene matrix, picks an optimal rank, and retrieves patients' cluster memberships and top NMF genes.

## 02_heatmap.R

This script generates heatmaps for top NMF genes and patient groups.

## 03_pca.R

This script performs PCA on the log2 normalized epigenes counts and visualizes patients NMF cluster membership the principal component space.