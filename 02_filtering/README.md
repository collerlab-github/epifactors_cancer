# Filtering

Prerequisites: 01_data_collection must be completed.

This folder contains scripts for filtering the raw counts data. The overall process includes:
1. Raw counts normalized using DESeq2 and saved normalized data frame. (Normalized)
2. Log2 transformation is performed and saved as another data frame. (Log2 normalized)
3. Subset the epigenes for both normalized and log2 normalized data frames.
3. Only keep epigenes in which at least 25% of patients have a normalized count >=5 (75th percentile >=5). Subset these genes from both normalized and log-normalized data and save as separate data frames.
4. Only keep epigenes which have a standard deviation greater than or equal to a cutoff, as determined by the log2 normalized 75th percentile filtered matrix. The cutoff is set to include at least 450 epigenes. Subset these epigenes from both normalized-percentile filtered and log2-normalized-percentile filtered data and save as separate data frames.
The final counts matrix used for further analyses is the one that is log2-normalized, percentile filtered, and standard deviation filtered.

Here are the descriptions for the scripts.


## 01_normalize_counts.R

This script normalizes the cancer patient raw counts dataframe using DESeq2 and log2 transforms using the built in log2 function. Both normalized and log2-normalized data frames are outputted.



## 02_epigenes_filter.R

This script filters normalized and log2 normalized data frames by epigenes (data frames 1 and 2), then subsets them based on epigenes with at least 25% of patients >=5 normalized counts (data frames 3 and 4), then subsets them based on the epigenes with log2 normalized expression standard deviation >= cutoff, which is set to include at least 450 epigenes (data frames 5 and 6). All data frames are outputted, as well as a standard deviation cutoff table.

