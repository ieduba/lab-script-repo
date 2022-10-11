cd ATACseq
find . -name '*.rmdup.bam' >> libs.txt	
samtools merge LSA1merge.bam `cat libs.txt`
macs2 callpeak --nomodel -t LSA1merge.bam -n peaksfrommerge --nolambda --keep-dup all --call-summits --slocal 10000
