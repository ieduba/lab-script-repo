#! /bin/bash
# SBATCH -N 1
# SBATCH -n 2

dir=$1
genelist=$2
genelistname=$3

source activate microc

cd $dir
for bedpe in merged_loops.bedpe noCTCF-loops.bedpe CTCF-loops.bedpe; do
	cat $bedpe | tail +3 | cut -f 1-3 > temp-anchors-$bedpe
	cat $bedpe | tail +3 | cut -f 4-6 >> temp-anchors-$bedpe

	bedtools sort -i temp-anchors-$bedpe > temp-anchors-st-$bedpe
	awk '{print "chr"$0}' temp-anchors-st-$bedpe > temp-anchors-chr-st-$bedpe
	#bedtools slop -i anchors-chr.st.bed -g /lustre/fs4/risc_lab/store/risc_data/downloaded/hg38/genome/chrom.sizes -b 5000 > anchors-chr.st.5kb.bed 
	bedtools intersect -a temp-anchors-chr-st-$bedpe -b $genelist -wb -bed > $genelistname-in-$bedpe
	rm temp-anchors*bedpe
done
cd -

