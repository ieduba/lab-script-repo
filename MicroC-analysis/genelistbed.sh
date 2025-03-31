genelist=$1 #csv from DESeq2
transcriptome=$2

name=`basename $genelist .csv`
cat $genelist | sed 's/"//g' | sed 's/,/\t/g' > $name\.tsv

awk 'NR==FNR{A[$1]=$1;next} ($7 in A){print}' $name\.tsv $transcriptome | cut -f 1-3,7,8 | sort -k 4 > temp.st.bed
tail -n +2 $name\.tsv | sort -k 1 | join temp.st.bed - -1 4 -2 1 | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$5"\t"$6"\t"$7}' | uniq | bedtools sort -i - > $name-genebody.bed
rm temp.st.bed
