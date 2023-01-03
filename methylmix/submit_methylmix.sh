#!/bin/bash

. /u/local/Modules/default/init/modules.sh

#$ -cwd
#$ -V
#$ -t 1:1:1
#$ -N methylmix
#$ -l h_rt=6:00:00,h_data=16G
#$ -pe shared 4
#$ -j y
#$ -o joblog/methylmix.$JOB_ID.$TASK_ID

module load R/4.1.0-BIO

counter=1
while read cancer
do
	if [ $SGE_TASK_ID = $counter ]
	then
		echo "$cancer"

		Rscript methylmix.R $cancer
		break
	fi
	counter=$((${counter}+1))
done < cancers.txt

echo "sleeping"
sleep 5m
