for dir in */aligned; do
	sbatch pseudorep-batch.sh $dir
done
