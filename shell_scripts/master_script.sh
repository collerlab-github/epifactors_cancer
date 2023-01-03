#!/bin/bash

if [ "$#" -ne 2 ];
then
	echo "Usage: bash master_script.sh [cancer] [purity]"
	exit 1
fi

cancer=$1
purity=$2

# data filtering
echo 'Data Filtering'
bash data_filter.sh $cancer $purity

# geneset_analysis
echo 'Geneset analysis'
bash geneset_analysis.sh $cancer

# single gene analysis
echo 'Single gene analysis'
bash single_gene_analysis.sh $cancer