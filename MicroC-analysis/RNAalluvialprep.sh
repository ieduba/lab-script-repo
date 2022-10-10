celltype1=$1
celltype2=$2

awk 'NR==FNR{A[$1]=$1;next} A[$1]{print}' $celltype1\_gene_tpm.csv $celltype2\_gene_tpm.csv > consensuslist.txt

for celltype in $celltype1 $celltype2; do
	for q in q1 q2 q3 q4; do
		newlist=$celltype\_gene_tpm_$q\.filt.csv
		awk 'NR==FNR{A[$1]=$1;next} A[$1]{print}' $celltype\_gene_tpm_$q\.csv consensuslist.txt > $newlist
		awk -v q=$q '{print $0, "'$q'"}' $newlist > $celltype\_gene_tpm_$q\.lab.csv
	done
	
	biglist=$celltype\_gene_tpm.lab.csv
	stbiglist=$celltype\_gene_tpm.st.csv	
	cat $celltype\_gene_tpm_q*.lab.csv > $biglist
	sort $biglist > $stbiglist
done




	
