if [ "$#" -ne 1 ];
then
	echo "Usage: bash master_script.sh [cancer]"
	exit 1
fi

if [[ $1 == "all" ]];
then
	echo 'all cancers'
	# all cancers
	cancer=('ACC' 'THCA' 'OV' 'LGG' 'PRAD' 'SKCM' 'UCEC'
		'KIRC' 'CRC' 'CESC' 'LIHC' 'SARC' 'HNSC' 'KIRP' 'GBM' 'PCPG'
		'LUAD' 'STAD' 'LUSC' 'ESCA' 'TGCT' 'ACC' 'BLCA' 'PAAD'
		)
else
	cancer=$1
fi

for c in ${cancer[@]}
do
	echo $c
	# 03 NMF
	echo 'NMF'
	cd 03_nmf
	Rscript 01_nmf.R $c
	echo 'Heatmap of Ranks 2 NMF top epigene expression'
	Rscript 02_top_genes_heatmap.R $c 2
	echo 'Heatmap of Ranks 3 NMF top epigene expression'
	Rscript 02_top_genes_heatmap.R $c 3
	echo 'PCA of patients epigene expression with NMF cluster projection'
	Rscript 03_top_genes_pca.R $c 2
	Rscript 03_top_genes_pca.R $c 3

	cd ../

	# 04 survival analysis
	echo 'Survival Analysis'
	cd 04_clinical_analysis
	Rscript 01_survival.R $c 2
	Rscript 01_survival.R $c 3
	cd ../
	echo 'finished'
done
