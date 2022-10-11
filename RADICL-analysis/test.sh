SAMPLENAME=$1

read1fastq=($SAMPLENAME/$SAMPLENAME\_S*_R1_001.fastq.gz)
echo "fastq path: $read1fastq"
read1=`basename $read1fastq ".fastq.gz"`
echo "fastq name: $read1"
