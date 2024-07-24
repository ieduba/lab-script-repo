#! /bin/bash
#SBATCH -N 1
#SBATCH -n 8

source activate mustache

python mustache.py -f scr-combo/scr-combo-mapped-filt_500.mcool -r 1000 -o scr-combo-mustache-1kb.tsv -st 0.7 -pt 0.1

python mustache.py -f dH1-combo/dH1-combo-mapped-filt_500.mcool -r 1000 -o dH1-combo-mustache-1kb.tsv -st 0.7 -pt 0.1
