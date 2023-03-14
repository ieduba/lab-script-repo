import numpy as np
import csv
from scipy import stats
from optparse import OptionParser

## SAMPLE NAMES AS INPUT ARGS ##
opts = OptionParser()
opts.add_option("--s1", help = "<sample 1> Name of first sample, as in list of contact distances up until -mapped...")
opts.add_option("--s2", help = "<sample 2> Name of second sample, as in list of contact distances up until -mapped...")
options, arguments = opts.parse_args()

## HARDCODED VARIABLES ##
orients = ['outward', 'tandemminus']
maxdist = 1500

## MAIN CODE ##
with open(f'{options.s1}-{options.s2}-KSresults.csv', 'w') as csvfile:
	for orient in orients:
		file1 = f'{options.s1}-mapped-filt-{orient}-dist.txt'
		file2 = f'{options.s2}-mapped-filt-{orient}-dist.txt'

		dists1 = np.loadtxt(file1, dtype = int)
		dists2 = np.loadtxt(file2, dtype = int)

		ks = stats.ks_2samp(dists1[dists1 < maxdist], dists2[dists2 < maxdist])
		csvfile.write(f'{orient}: {ks}\n')
