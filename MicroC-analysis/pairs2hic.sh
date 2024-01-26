#! /bin/bash
# SBATCH -N 1
# SBATCH -n 16

echo This is $SLURM_JOB_ID
echo Test permission

source activate "microc"

pairs=$1
echo $pairs
juicerpath=/rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/juicer
juicertools=$juicerpath/scripts/common/juicer_tools

#pairs=`echo $rawpair | sed 's/mapped/mapped-filt/'`
#awk '$2 !~ /chrM|chrY|chrU.+|chr.+_.+/ && $4 !~ /chrM|chrY|chrU.+|chr.+_.+/ {print $0}' $rawpair > $pairs
hictxt=`echo $pairs | sed 's/.pairs/.hic.txt/'`
echo $hictxt
hic=`echo $hictxt | sed 's/.hic.txt/-mq30.hic/'`
echo $hic

echo making txt file
awk '{print $6,"\t", $2,"\t", $3,"\t0\t",$7,"\t",$4,"\t",$5,"\t1"}' $pairs > $hictxt
echo making hic file
bash $juicertools pre -q 30 $pairs $hic /rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/chrom.sizes

