for bam in *.bam; do
	out=`echo $bam | sed 's/bam/len.txt/'`
	samtools view -f 35 $bam | cut -f 9 > $out
done
