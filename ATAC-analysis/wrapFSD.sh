cd LSA2/mergedLanes
files=`ls -d */`

cd ~/labscratch/npagane/pipes/irene_senescence
for file in $files; do
        fileshort=`echo $file | cut -d'-' -f1`
        cd $fileshort
        cp *.rmdup.bam ~/scripts/LSA2/pipelineRerun/chromState/${file}genomewide-frag_E025.bam
        cd ..
done

cd ~/scripts
for file in $files; do
	sbatch runFSD.sh $file
done
