from itertools import combinations
import re
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import pandas as pd
import cooler
import cooltools
import bioframe
from optparse import OptionParser

plt.rcParams.update({'font.size': 18})

opts = OptionParser()
opts.add_option('-m', help = '<.mcool> paths to mcool files, separated by comma and no space')
opts.add_option('-r', help = '<resolution> resolution in bp of mcool to use for plotting P(s) (usually 10 kb or 1 kb)')
options, arguments = opts.parse_args()

def makecvd(clr, hg38_arms):
	## select only those chromosomes available in cooler
	hg38_arms = hg38_arms[hg38_arms.chrom.isin(clr.chromnames)].reset_index(drop=True)

	## cvd == contacts-vs-distance
	cvd_smooth_agg = cooltools.expected_cis(
		clr=clr,
		view_df=hg38_arms,
		smooth=True,
		aggregate_smoothed=True,
		nproc=1
	)

	cvd_smooth_agg['s_bp'] = cvd_smooth_agg['dist']* resolution
	cvd_smooth_agg['balanced.avg.smoothed.agg'].loc[cvd_smooth_agg['dist'] < 2] = np.nan

	## Just take a single value for each genomic separation
	cvd_merged = cvd_smooth_agg.drop_duplicates(subset=['dist'])[['s_bp', 'balanced.avg.smoothed.agg']]
	## Calculate derivative in log-log space
	der = np.gradient(np.log(cvd_merged['balanced.avg.smoothed.agg']),
			  np.log(cvd_merged['s_bp']))
	
	return (cvd_merged, der)


#### MAIN CODE ####

mc_files = options.m.split(',')
mc_names = re.sub(",","_",options.m)
resolution = int(options.r)

## Use bioframe to fetch the genomic features from the UCSC.
hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_cens = bioframe.fetch_centromeres('hg38')
## create a view with chromosome arms using chromosome sizes and definition of centromeres
hg38_arms = bioframe.make_chromarms(hg38_chromsizes,  hg38_cens)

## Use makecvd funciton to make contact vs distance data for each mcool, then plot

f, axs = plt.subplots(
	figsize=(8,11),
	nrows=2,
	gridspec_kw={'height_ratios':[6,2]},
	sharex=True)

for mc in mc_files:
	print(mc)
	name = f'{mc.split("-")[0]}-{mc.split("-")[1]}'	
	clr = cooler.Cooler(f'{mc}::/resolutions/{resolution}')
	(cvd_merged, der) = makecvd(clr, hg38_arms)

	ax = axs[0]
	ax.loglog(cvd_merged['s_bp'], cvd_merged['balanced.avg.smoothed.agg'], '-', markersize=5, label=name)
	ax.set(ylabel='IC contact frequency', xlim=(1e3,1e8))
	ax.set_aspect(1.0)
	ax.grid(lw=1)

	ax = axs[1]
	ax.semilogx(cvd_merged['s_bp'], der)
	ax.set(xlabel='separation, bp', ylabel='slope')
	ax.grid(lw=1)

axs[0].legend()
plt.savefig(f'{mc_names}_Ps_deriv.pdf', format='pdf')
