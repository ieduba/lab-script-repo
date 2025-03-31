#! /bin/bash
#SBATCH -N 1
#SBATCH -n 4

source activate microc

# annotate gene list with compartment score, collect summary stats
echo 'average_diff_cscore gene_list N_genes' > gene-cscore-summary.csv
for genelist in *genebody.bed; do
	name=`echo $genelist | sed 's/noheader-//' | sed 's/.bed//'`
	bedtools intersect -a $genelist -b scr-corr-cscore.bedgraph -wao | cut -f 1,2,3,4,5,6,7,11 > $name-scr-cscore.temp
	bedtools intersect -a $genelist -b low-corr-cscore.bedgraph -wao | cut -f 11 > $name-low-cscore.temp
	echo 'chr	pos1	pos2	gene_ID	gene_type	base_mean	log2fold_change	scr_cscore	low_cscore' > $name-cscore.bedgraph
	paste $name-scr-cscore.temp $name-low-cscore.temp >> $name-cscore.bedgraph
	rm $name*temp

	Ngene=`wc -l $genelist`
	avgdc=`awk '$9~/0./ {n++;diff+=$9-$8} END {print diff/n}' $name-cscore.bedgraph`
	echo $avgdc $Ngene >> gene-cscore-summary.csv
done
