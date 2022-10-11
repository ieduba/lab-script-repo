cd ~/linker-histone/Bcell
directories=`ls -d Res*/`
cd ~/scripts

for directory in $directories; do
        sbatch -p risc,hpc RICCannoFSD.sh $directory
done
