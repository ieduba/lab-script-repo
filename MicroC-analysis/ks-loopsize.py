import numpy as np
import csv
from scipy import stats

## HARDCODED VARIABLES ##
maxdist = 1000000
scronly='scr-dH1-all-1-only-loops-subsample.bed-dist'
lowonly='scr-dH1-all-2-only-loops-subsample.bed-dist'
common='scr-dH1-all-common-loops-subsample.bed-dist'

## MAIN CODE ##
with open('subsampled-looplength-KSresults.csv', 'w') as csvfile:
	dists1 = np.loadtxt(scronly, dtype = int)
	dists2 = np.loadtxt(lowonly, dtype = int)
	dists3 = np.loadtxt(common, dtype = int)

	ks1 = stats.ks_2samp(dists1[dists1 < maxdist], dists2[dists2 < maxdist])
	ks2 = stats.ks_2samp(dists1[dists1 < maxdist], dists3[dists3 < maxdist])
	ks3 = stats.ks_2samp(dists2[dists2 < maxdist], dists3[dists3 < maxdist])

	csvfile.write(f'scr only vs low only: {ks1}\nscr only vs common: {ks2}\nlow only vs common: {ks3}\n')
