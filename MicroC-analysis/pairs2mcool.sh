#! /bin/bash

#SBATCH -N 1
#SBATCH -n 8

source activate microc

pairs=$1 #path to pairs file

echo $pairs

cool=`echo $pairs | sed 's/-mapped.pairs/_5000.cool/'`
mcool=`echo $cool | sed 's/.cool/.mcool/'`

echo making cool file $cool
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly hg38 /lustre/fs4/risc_lab/store/risc_data/downloaded/hg38/genome/chrom.sizes:5000 $pairs $cool

echo making mcool file $mcool
cooler zoomify --nproc 8 --balance --out $mcool --resolutions 10000,5000 $cool
