bed1=$1
bed2=$2

name1=`echo $bed1 | sed 's/.bed//'`
name2=`echo $bed2 | sed 's/.bed//'`

bedtools sort -i $bed1 > $name1\.st.bed
bedtools sort -i $bed2 > $name2\.st.bed
bedtools intersect -a $name1\.st.bed -b $name2\.st.bed -wa -wb > $name1-$name2-one-overlap.bed
awk '{OFS="\t"}{print $4,$5,$6,$1,$2,$3}' $name1-$name2-one-overlap.bed > $name1-one-overlap.bed
awk '{OFS="\t"}{print $10,$11,$12,$7,$8,$9}' $name1-$name2-one-overlap.bed > $name2-one-overlap.bed
bedtools sort -i $name1-one-overlap.bed > $name1-one-overlap.st.bed
bedtools sort -i $name2-one-overlap.bed > $name2-one-overlap.st.bed
bedtools intersect -a $name1-one-overlap.st.bed -b $name2-one-overlap.st.bed -wa -wb > $name1-$name2-two-overlap.bed
uniq $name1-$name2-two-overlap.bed > $name1-$name2-two-overlap.uq.bed
awk '{OFS="\t"} ($2 == $8 && $5 == $11){print $4,$5,$6,$1,$2,$3}' $name1-$name2-two-overlap.uq.bed > $name1-$name2-two-real-overlap.bed
sort $name1-$name2-two-real-overlap.bed | uniq > $name1-$name2-two-real-overlap.uq.bed
bedtools sort -i $name1-$name2-two-real-overlap.uq.bed > $name1-$name2-common.bed

#housekeeping
rm $name1*overlap*.bed
rm $name2*overlap*.bed
