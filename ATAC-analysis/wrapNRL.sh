cd LSA2/mergedLanes
files=`ls -d */`
cd ~/scripts

for file in $files; do
	sbatch NRLfinder.sh $file
done
