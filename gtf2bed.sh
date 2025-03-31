gtf=$1
filename=`basename $gtf .gtf`

cut -f 1 $gtf > chr
cut -f 3 $gtf > feat
cut -f 4 $gtf > pos1
cut -f 5 $gtf > pos2
cut -f 6 $gtf > score
cut -f 7 $gtf > strand
cut -f 9 $gtf | cut -f 1 -d ';' | cut -f 2 -d ' ' | sed 's/"//g' > name
cut -f 9 $gtf | cut -f 3 -d ';' | cut -f 3 -d ' ' | sed 's/"//g' > transtype

paste chr pos1 pos2 name score strand transtype > $filename.bed
rm chr feat pos1 pos2 name score strand transtype

awk '/protein_coding/' $filename.bed | uniq > $filename-proteincoding.bed
