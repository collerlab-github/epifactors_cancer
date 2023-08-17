# if [ "$#" -ne 2 ];
# then
# 	echo "Usage: bash master_script.sh [project] [cancer]"
# 	exit 1
# fi

project="TCGA"

cancer=('BRCA' 'THCA' 'OV' 'LGG' 'PRAD' 'SKCM' 'UCEC'
'KIRC' 'CRC' 'CESC' 'LIHC' 'SARC' 'HNSC' 'KIRP' 'GBM' 'PCPG'
'LUAD' 'STAD' 'LUSC' 'ESCA' 'TGCT' 'ACC' 'BLCA' 'PAAD')

for c in ${cancer[@]}
do
	echo $c
	# 05 survminer
	echo 'Single Gene Survminer Analysis'
	/Library/Frameworks/R.framework/Resources/bin/Rscript survminer_random.R $project $c 
	echo 'finished'
done