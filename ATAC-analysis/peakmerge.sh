##takes first three columns of narrowpeak file (chr#, peak start, peak end) and merges peaks that overlap. Ready for creating count matrix.

mkdir bed3s
for filename in */peakCalls/*.narrowPeak; do
	bedfile=`basename $filename | sed 's/.narrowPeak/.bed/'`
	cat $filename | cut -f 1,2,3 > $bedfile
	mv $bedfile bed3s
done
cd bed3s
cat *.bed > masterpeaklist.bed
bedtools sort -i masterpeaklist.bed > masterpeaklist.st.bed
bedtools merge -i masterpeaklist.st.bed -c 1 -o count > mergedpeaks.bed
wc -l *_peaks.bed
wc -l masterpeaklist.bed
wc -l mergedpeaks.bed
cd ..
