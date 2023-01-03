if [ "$#" -ne 2 ];
then
	echo "Usage: bash master_script.sh [project] [cancer]"
	exit 1
fi

project=$1

if [[ $2 == "all" ]];
then
	echo 'all cancers'
	# all cancers
	cancer=('ACC' 'THCA' 'OV' 'LGG' 'PRAD' 'SKCM' 'UCEC'
		'KIRC' 'CRC' 'CESC' 'LIHC' 'SARC' 'HNSC' 'KIRP' 'GBM' 'PCPG'
		'LUAD' 'STAD' 'LUSC' 'ESCA' 'TGCT' 'ACC' 'BLCA' 'PAAD'
		)
else
	cancer=$2
fi


for c in ${cancer[@]}
do
	echo $c
	# 05 survminer
	echo 'Single Gene Survminer Analysis'
	cd 05_gene_clinical_prediction_analysis
	Rscript 01_survminer.R $project $c 
	cd ../
	echo 'finished'
done


