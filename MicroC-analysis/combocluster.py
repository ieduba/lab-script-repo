import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from optparse import OptionParser
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

opts = OptionParser()
usage = 'usage: %prog [options] [inputs]'
opts = OptionParser(usage=usage)
opts.add_option("-c", help = "<celltype> Name of cell type as in csv names") 
options, arguments = opts.parse_args()

celltype = options.c

allscaled = np.zeros((0,152))
orients = ['tandemminus', 'tandemplus', 'in', 'out']
for i in range(0,3):
	datafile = f'{celltype}-{orients[i]}dist.csv'
	data = np.loadtxt(datafile, dtype = int, delimiter = ',')
	scaled_data = np.zeros((data.shape))
	scaled_rown = np.zeros((data.shape[0], data.shape[1]+1))
	scaled_or = np.zeros((data.shape[0], data.shape[1]+2))
	for r in range(0,len(data)):
		rowr = data[r]
		if max(rowr) == 0:
			scaled_rown[r] = np.insert(scaled_data[r], 0, r)
			continue
		scaled_data[r] = (rowr - min(rowr)) / (max(rowr) - min(rowr))
		scaled_rown[r] = np.insert(scaled_data[r], 0, r)
	scaled_or = np.insert(scaled_rown, 0, i, axis = 1)
	print(scaled_or.shape)
	print(allscaled.shape)
	allscaled = np.append(allscaled, scaled_or,0)
print(allscaled.shape)
np.savetxt('{celltype}-allscaled.csv',allscaled, delimiter =',')

kmeans = KMeans(init = 'kmeans++', n_clusters = 16, n_init = 50, max_iter = 500
kmeans.fit(allscaled)
