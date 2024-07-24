import numpy as np
import csv
from optparse import OptionParser

opts = OptionParser()
opts.add_option("-c", help = "<control> 2 column csv  of RNAseq gene names and base means for control genes -- sorted by basemean")
opts.add_option("-u", help = "<upreg> 2 column csv of RNAseq gene names and base means for upregulated genes -- sorted by basemean")
options, arguments = opts.parse_args()

ctrl = np.loadtxt(options.c, dtype=str)
upreg = np.loadtxt(options.u, dtype=str)

ctrlgenes = list(ctrl[:,0])
upgenes = list(upreg[:,0])
ctrlmeans = list(ctrl[:,1].astype(float))
upmeans = list(upreg[:,1].astype(float))

print(f'original control mean: {np.mean(ctrlmeans)}')
print(f'upregulated mean: {np.mean(upmeans)}')

if np.mean(ctrlmeans) > np.mean(upmeans):
	while np.mean(ctrlmeans) > np.mean(upmeans):
		ctrlmeans.remove(ctrlmeans[-1])
		ctrlgenes.remove(ctrlgenes[-1])
else:
	while np.mean(ctrlmeans) < np.mean(upmeans):
		ctrlmeans.remove(ctrlmeans[0])
		ctrlgenes.remove(ctrlgenes[0])

mean_corrected_ctrl = np.column_stack((ctrlgenes,ctrlmeans))

with open(f'up-mean-{options.c}', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(mean_corrected_ctrl)

print(f'new control mean: {np.mean(ctrlmeans)}')
