#! /bin/bash
# SBATCH -N 1
# SBATCH -n 20

echo This is $SLURM_JOB_ID
echo Test permission

source activate "microc"

rawpairs=$1
#echo $pairs
juicerpath=/rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/juicer
juicertools=$juicerpath/scripts/common/juicer_tools

pairs=`echo $rawpairs | sed 's/mapped/mapped-filt/'`
#awk '$2 !~ /chrM|chrY|chrU.+|chr.+_.+/ && $4 !~ /chrM|chrY|chrU.+|chr.+_.+/ {print $0}' $rawpairs > $pairs

noinpairs=`echo $pairs | sed 's/.pairs/-noin.pairs/'`
awk '$6 != "+" || $7 != "-" { print $0 }' $pairs > $noinpairs

hictxt=`echo $noinpairs | sed 's/.pairs/.hic.txt/'`
echo $hictxt
hic=`echo $hictxt | sed 's/.hic.txt/-mq30.hic/'`
echo $hic

#echo making txt file
awk '{print $6,"\t", $2,"\t", $3,"\t0\t",$7,"\t",$4,"\t",$5,"\t1"}' $pairs > $hictxt

echo making hic file
bash $juicertools pre -q 30 -r 1000000,100000,25000,10000,5000,1000,500 -t /tmp $pairs $hic /rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/chrom.sizes
#bash $juicertools pre -q 30 -t /tmp $pairs $hic /rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm10/genome/chrom.sizes
