cd ~/linker-histone/Bcell
files=`ls BJ*.bam`
cd ~/scripts

for file in $files; do
	sbatch -p risc,hpc RICCannotations.sh $file
done
