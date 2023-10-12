genelist=$1
transcriptome=$2
out=`echo $genelist | sed 's/.csv/-genebody.bed/'`
awk 'NR==FNR{A[$1]=$1;next} ($7 in A){print}' $genelist $transcriptome | cut -f 1-3,7 > $out

