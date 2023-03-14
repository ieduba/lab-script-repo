celltype1=$1
celltype2=$2

awk 'NR==FNR{A[$1]=$1;next} A[$1]{print}' $celltype1\_st_genes.csv $celltype2\_st_genes.csv > consensuslist.txt

for celltype in $celltype1 $celltype2; do
	for q in q1 q2 q3 q4; do
		newlist=$celltype\_st_genes_$q\.filt.csv
		awk 'NR==FNR{A[$1]=$1;next} A[$1]{print}' $celltype\_st_genes_$q\.csv consensuslist.txt > $newlist
		awk -v q=$q '{print $0, "'$q'"}' $newlist > $celltype\_st_gene_$q\.lab.csv
	done
	
	biglist=$celltype\_st_gene.lab.csv
	stbiglist=$celltype\_st_gene.lab.st.csv	
	cat $celltype\_st_gene_q*.lab.csv > $biglist
	sort $biglist > $stbiglist
done




	
