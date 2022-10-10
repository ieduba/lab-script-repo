#! /bin/bash
# SBATCH -N 1
# SBATCH -n 12

echo This is job $SLURM_JOB_ID
echo Test permission

python -u binned-contactprob.py -b 750kb-mm10-filt.bed -p mESC_WT_CKDL210025175-1a-2_HVWV3DSX2_L1-mapped-filt.pairs -c 12
python -u binned-contactprob.py -b 750kb-mm10-filt.bed -p mESC_low_CKDL210025175-1a-4_HVWV3DSX2_L1-mapped-filt.pairs -c 12

