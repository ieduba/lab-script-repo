for res in 50000 100000; do
	bed=hg38-$res\.bed
	for binfile in $res*txt; do
		out=`echo $binfile | sed 's/.txt/.bed/'`
		awk 'NR==FNR{A[$1]=$1;next} ($1 in A){print}' $binfile $bed | cut -f 2-4 > $out
	done
done
