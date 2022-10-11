import pybedtools as pbt
import subprocess

#evetually make these input arguments
genome = "hg19"
kb = 1000

binsize = kb*1000

#create pybedtools objects
chromsizes = pbt.BedTool("/rugpfs/fs0/risc_lab/store/risc_data/downloaded/{0}/genome/chrom.sizes".format(genome))
blacklist = pbt.BedTool("/rugpfs/fs0/risc_lab/store/risc_data/downloaded/{0}/blacklist/{0}-blacklist.v2.bed".format(genome))
#fasta = pbt.BedTool("/rugpfs/fs0/risc_lab/store/risc_data/downloaded/{0}/genome/Sequence/WholeGenomeFasta/genome.fa".format(genome))

##make bins of given size across genome, removing blacklisted bins
blrm = chromsizes.window_maker(g=chromsizes.fn, w=binsize).subtract(blacklist)
##troubleshoot this - getting python error that fasta is corrupted
#gc = fasta.nucleotide_content(blrm)

#system command workaround for getting GC content. then take just GC column, remove top line added by nuc, and add unique ID to each row as fifth column, sort by GC content, and remove wonky chromosomes (have underscores in name)
blrm.saveas("blrm.bed")
subprocess.run("bedtools nuc -fi /rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/genome/Sequence/WholeGenomeFasta/genome.fa -bed blrm.bed | cut -f 1-3,5 | awk '!/usercol/' | awk '3{$5=NR}1' | sort -k4 | awk '$1 !~ /_/' > blrm.gc.st.bed", shell=True)

