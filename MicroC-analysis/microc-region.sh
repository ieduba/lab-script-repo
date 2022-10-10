
# make file names. bedin is ChIP regions
bedin=$1
sample=$2
bedname=`basename $bedin | sed 's/.bed//'`

# making more file names
bam=$sample-mapped.PT.bam
bamout=$bedname-$bam
pair=$sample-mapped.pairs
pairout=$bedname-$pair

# make new bams with only reads that fall in heterochromatin or euchromatin regions 
bedtools intersect -a $bam -b $bedin -wa > $bamout

# turn bams to sams for awk manipluation
samtools view $bamout > temp.sam

# filter .pairs files to include only lines from intersected bams (matched by names) to make
# 	region .pairs file for contact probability analysis
awk 'NR==FNR{A[$1]=$2;next} A[$1]{print}' temp.sam $pair > $pairout

rm temp.sam

