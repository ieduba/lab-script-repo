#! /bin/bash
# SBATCH -N 1
# SBATCH -n 2

dir=$1
genelist=$2
genelistname=$3

source activate microc

cd $dir
for bedpe in differential_loops1.bedpe differential_loops2.bedpe; do
	cat $bedpe | tail +3 | cut -f 1-3 > temp-anchors-$bedpe
	cat $bedpe | tail +3 | cut -f 4-6 >> temp-anchors-$bedpe

	bedtools sort -i temp-anchors-$bedpe > temp-anchors-st-$bedpe
	awk '{print "chr"$0}' temp-anchors-st-$bedpe > temp-anchors-chr-st-$bedpe
	bedtools slop -i temp-anchors-chr-st-$bedpe -g /lustre/fs4/risc_lab/store/risc_data/downloaded/hg38/genome/chrom.sizes -b 5000 > temp-anchors-chr-st-5kb-$bedpe 
	bedtools intersect -a temp-anchors-chr-st-5kb-$bedpe -b $genelist -wb -bed > $genelistname-5kb-$bedpe
	rm temp-anchors*bedpe
done
cd -

