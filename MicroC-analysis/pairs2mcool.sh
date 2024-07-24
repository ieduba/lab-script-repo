#! /bin/bash

#SBATCH -N 1
#SBATCH -n 8

source activate microc

rawpairs=$1 #path to pairs file

echo $rawpairs

pairs=`echo $rawpairs | sed 's/mapped/mapped-filt/'`
#awk '$2 !~ /chrM|chrY|chrU.+|chr.+_.+/ && $4 !~ /chrM|chrY|chrU.+|chr.+_.+/ {print $0}' $rawpairs > $pairs

noinpairs=`echo $pairs | sed 's/.pairs/-noin.pairs/'`
#awk '$6 != "+" || $7 != "-" { print $0 }' $pairs > $noinpairs

cool=`echo $noinpairs | sed 's/-noin.pairs/_500.cool/'`
mcool=`echo $cool | sed 's/.cool/.mcool/'`

echo making cool file $cool
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly hg38 /lustre/fs4/risc_lab/store/risc_data/downloaded/hg38/genome/chrom.sizes:500 $pairs $cool

echo making mcool file $mcool
cooler zoomify --nproc 8 --balance --out $mcool --resolutions 1000000,100000,10000,5000,1000,500 $cool
