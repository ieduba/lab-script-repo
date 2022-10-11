genome=mm10
kb=100

binsize=$(($kb*1000))

#divide genome into 100kb windows and remove blacklist
bedtools makewindows -g /rugpfs/fs0/risc_lab/store/risc_data/downloaded/$genome/genome/chrom.sizes -w $binsize > $genome_$kb\kbBin_$genome.bed
bedtools subtract -a $genome_$kb\kbBin.bed -b /rugpfs/fs0/risc_lab/store/risc_data/downloaded/$genome/blacklist/$genome-blacklist.v2.bed > $genome_$kb\kbBin_$genome.blrm.bed
#get GC content across the binned genome (remove other columns and row added by nuc)
#bedtools nuc -fi /rugpfs/fs0/risc_lab/store/risc_data/downloaded/$genome/genome/Sequence/WholeGenomeFasta/genome.fa -bed $genome_$kb\kbBin.bed > $genome_$kb\kbBin.nuc.bed
#cut -f 1-3,5 $genome_$kb\kbBin.blrm.nuc.bed > $genome_$kb\kbBin_GC.blrm.bed
#awk '!/usercol/' $genome_$kb\kbBin_GC.blrm.bed > $genome_$kb\kbBin_GC2.blrm.bed
#add column with unique ID for each bin
#awk '3{$5=NR}1' $genome_$kb\kbBin_GC2.blrm.bed > $genome_$kb\kbBin_GC.blrm.ID.bed
