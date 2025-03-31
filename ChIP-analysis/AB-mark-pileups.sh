#! /bin/bash
#SBATCH -N 1
#SBATCH -n 16

source activate cut_tag

# pileups of all the marks around A and B compartment regions (merged chunks so bed regions correspond to compartment boundaries)
for bw in H1low_H1p2_r1_S38_001.trim.st.all.blft.qft.rmdup.bw \
	Scrm_H1p2_r3_S15_001.trim.st.all.blft.qft.rmdup.bw \
	H1low_H3k27me3_r3_S27_001.trim.st.all.blft.qft.rmdup.bw \
	Scrm_H3k27me3_r3_S3_001.trim.st.all.blft.qft.rmdup.bw \
	Scrm_H3k36me2_r2_S5_001.trim.st.all.blft.qft.rmdup.bw \
	H1low_H3k9me3_r3_S33_001.trim.st.all.blft.qft.rmdup.bw \
	Scrm_H3k9me3_r1_S7_001.trim.st.all.blft.qft.rmdup.bw \
	H1low_H3_r3_S43_001.trim.st.all.blft.qft.rmdup.bw; do	
	name=`echo $bw | sed 's/_001.trim.st.all.blft.qft.rmdup.bw//'`

	computeMatrix scale-regions -S $bw -R KH3-A-corrected-merge.bed KH3-B-corrected-merge.bed \
	--beforeRegionStartLength 25000 --afterRegionStartLength 25000 --binSize 1000 --regionBodyLength 100000 \
	-p 16 -o K562_$name\_CT_AB.mat.gz

	plotHeatmap -m K562_$name\_CT_AB.mat.gz -out K562_$name\_CT_AB-pileup.pdf --colorMap "viridis" --startLabel "start" --endLabel "end" --heatmapHeight 10

done

# pileups of all marks from previous C&T round around A and B compartments (this for loop and the above one can be combined)
for bw in H1low_H3k27ac_r1_S16_001.trim.st.all.blft.qft.rmdup.sorted_extendreads.bw \
        H1low_H3k4me1_r1_S13_001.trim.st.all.blft.qft.rmdup.sorted_extendreads.bw \
	H1low_H3k4me3_r1_S15_001.trim.st.all.blft.qft.rmdup.sorted_extendreads.bw \
	Scrm_H3k27ac_r2_S12_001.trim.st.all.blft.qft.rmdup.sorted_extendreads.bw \
	Scrm_H3k4me1_r1_S5_001.trim.st.all.blft.qft.rmdup.sorted_extendreads.bw \
	Scrm_H3k4me3_r1_S7_001.trim.st.all.blft.qft.rmdup.sorted_extendreads.bw; do
        name=`echo $bw | sed 's/_001.trim.st.all.blft.qft.rmdup.sorted_extendedreads.bw//'`

        computeMatrix scale-regions -S $bw -R KH3-A-corrected-merge.bed KH3-B-corrected-merge.bed \
        --beforeRegionStartLength 25000 --afterRegionStartLength 25000 --binSize 1000 --regionBodyLength 100000 \
        -p 16 -o K562_$name\_CT_AB.mat.gz

        plotHeatmap -m K562_$name\_CT_AB.mat.gz -out K562_$name\_CT_AB-pileup.pdf --colorMap "viridis" --startLabel "start" --endLabel "end" --heatmapHeight 10

done

# piileups of H1.2 in heterochromatic peaks. should I change this to merged bed files of peaks?
for bed in Scrm_H3k9me3_r1_S7/Scrm_H3k9me3_r1_S7_001.trim.st.all.blft.qft.rmdup-W200-G2000-FDR0.1-island.st.bed \
	Scrm_H3k27me3_r3_S3/Scrm_H3k27me3_r3_S3_001.trim.st.all.blft.qft.rmdup-W200-G2000-FDR0.1-island.st.bed \
        Scrm_H3k27ac_r1_S19/Scrm_H3k27ac_r1_S19_001.trim.st.all.blft.qft.rmdup-W200-G2000-FDR0.1-island.st.bed; do

	for bw in Scrm_H1p2_r3_S15_001.trim.st.all.blft.qft.rmdup.bw \
		Scrm_H1p2_r3_S15_001.trim.st.all.blft.qft.rmdup.bw \
		H1low_H1p2_r1_S38_001.trim.st.all.blft.qft.rmdup.bw \
		H1low_H1p2_r2_S39_001.trim.st.all.blft.qft.rmdup.bw \
		H1low_H1p2_r3_S40_001.trim.st.all.blft.qft.rmdup.bw; do
		
		name=`basename $bw _001.trim.st.all.blft.qft.rmdup.bw`-over-`basename $bed _001.trim.st.all.blft.qft.rmdup-W200-G2000-FDR0.1-island.st.bed`
		computeMatrix scale-regions -S $bw -R $bed --beforeRegionStartLength 2500 --afterRegionStartLength 2500 \
		--binSize 100 -p 16 -o $name\.mat.gz
 
		plotHeatmap -m $name\.mat.gz -out $name-pileup.pdf --colorMap "viridis" --startLabel "start" --endLabel "end" --heatmapHeight 10
	done
done

#multiBamSummary bins --bamfiles *rmdup.st.bam --outFileName bam-summary-results.npz
 
#plotCorrelation -in bam-summary-results.npz --corMethod pearson --whatToPlot heatmap --plotTitle "Pearson Correlation Heatmap" --colorMap 'viridis' --plotNumbers --plotFileFormat pdf -o bams-correlation_heatmap.pdf
