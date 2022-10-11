#Cp /rugpfs/fs0/risc_lab/scratch/ascortea/LSA5/*/intermediates/*.qft.bam .
#cp /rugpfs/fs0/risc_lab/scratch/npagane/pipes/andrew_LSA3/*/intermediates/*.qft.bam .

echo -e "FILE NAME \t SEQUENCING DEPTH \t UNIQUE READS" >> thresholds.txt

for file in *qft.bam; do
	sbatch -p risc,hpc preseq.sh $file
done
