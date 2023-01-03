#!/bin/bash

. /u/local/Modules/default/init/modules.sh

#$ -cwd
#$ -V
#$ -t 1:1:1
#$ -N get_methyl
#$ -l h_rt=3:00:00,h_data=64G
#$ -j y
#$ -o joblog/get_methyl.$JOB_ID.$TASK_ID

module load R/4.1.0-BIO

counter=1
while read cancer
do
	if [ $SGE_TASK_ID = $counter ]
	then
		echo "$cancer"

		Rscript get_methylation_data.R $cancer
		break
	fi
	counter=$((${counter}+1))
done < cancers.txt

echo "sleeping"
sleep 5m
