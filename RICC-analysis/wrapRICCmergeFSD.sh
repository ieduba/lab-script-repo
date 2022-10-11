cd ~/linker-histone/MESR3/corrected
directories=`ls -d */`
cd ~/scripts

for directory in $directories; do
        sbatch -p risc,hpc RICCmergeFSD.sh $directory
done
