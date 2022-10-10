for file in */aligned/merged_nodups.txt; do
	outfile=`echo $file | sed 's/merged_nodups.txt/hic.summary/'`
	awk '{print $15,"\t", $2,"\t", $3,"\t",$1,"\t",$6,"\t",$7,"\t",$5}' $file > $outfile
done
