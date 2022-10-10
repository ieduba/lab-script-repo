import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from optparse import OptionParser
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from matplotlib import cm

########## DEFINE FUNCTIONS ##########
## function to iterate through k values to determine optimum number of clusters using SSE and silhouette score
def optimization():
	print('running clustering optimization')
	sse = []
	silhouette = []
	for k in range(2,16):
		kmeans = KMeans(init = "k-means++", n_clusters = k, n_init = 50, max_iter = 500)
		kmeans.fit(scaled_data)
		sse.append(kmeans.inertia_)
		score = silhouette_score(scaled_data, kmeans.labels_)
		silhouette.append(score)

	np.savetxt(f'{celltype}-{orient}-SSE.csv', sse, delimiter=',')	
	np.savetxt(f'{celltype}-{orient}-silhouette.csv', silhouette, delimiter=',')
	
	## plot SSE vs k to look for elbow
	plt.style.use('fivethirtyeight')
	plt.plot(range(2,16),sse)
	plt.xticks(range(2,16))
	plt.xlabel('Number of clusters')
	plt.ylabel('SSE')
	plt.savefig(f'{celltype}-{orient}-SSE.png', bbox_inches='tight')
	plt.clf()	

	## plot silhouette score vs k to look for high score
	plt.style.use('fivethirtyeight')
	plt.plot(range(2,16), silhouette)
	plt.xticks(range(2,16))
	plt.xlabel('Number of clusters')
	plt.ylabel('Silhouette coefficient')
	plt.savefig(f'{celltype}-{orient}-silhouette.png', bbox_inches='tight')
	plt.clf()

## function to run kmeans clustering and PCA, saving relevant results and plots
def clustering(scaled, wid):
	## run kmeans with 50 kmeans++ initializations and 500 max iterations to convergence
	kmeans = KMeans(init = 'k-means++', n_clusters = nclust, n_init = 50, max_iter = 500)
	kmeans.fit(scaled)

	## add cluster labels to array w IDs, aggregate by cluster to count # bins in each
	ids = kmeans.labels_.reshape((kmeans.labels_.shape[0],1))
	datawids = np.append(wid, ids, 1) #so now cluster is nhbin+2nd column
	countsdf = pd.DataFrame(datawids).groupby(wid.shape[1]).agg('count')

	## run PCA (to get PCs, explained variance, and means) and flattened PCA (for plotting)
	pca = PCA()
	pcadata = pca.fit(scaled)
	transformed = pca.fit_transform(scaled)

	## save cluster and PCA data 
	countsdf.to_csv(f'{celltype}-{orient}-{nclust}-clustersizes.csv', columns = (0,1))
	np.savetxt(f'{celltype}-{orient}-{nclust}-clustercenters.csv', kmeans.cluster_centers_, delimiter=',')
	np.savetxt(f'{celltype}-{orient}-{nclust}-binclusters.csv', datawids, delimiter=',')
	qc = [f'min SSE: {kmeans.inertia_}', f'number of iterations: {kmeans.n_iter_}', f'cluster means shape: {kmeans.cluster_centers_.shape}']
	with open(f'{celltype}-{orient}-{nclust}-clusterQC.txt', 'w') as f:
		for line in qc:
			f.write(line)
			f.write('\n')
#	np.savetxt(f'{celltype}-{orient}-components.csv', pcadata.components_, delimiter = ',')
#	np.savetxt(f'{celltype}-{orient}-explainedvariance.csv', pcadata.explained_variance_, delimiter = ',')
#	np.savetxt(f'{celltype}-{orient}-PCAtransform.csv', transformed, delimiter = ',')
#	np.savetxt(f'{celltype}-{orient}-PCAmeans.csv', pcadata.mean_, delimiter = ',')	

	## return labeled data, cluster centers, cluster counts, and flattened PCA for plotting
	return(datawids, kmeans.cluster_centers_, countsdf, transformed)

## function to plot cluster centers (means)
def clusterplot(clusters, cols, counts):
	## set up number of rows/columns for subplots based on number of clusters
	if nclust <= 4: 
		nr = 2
	elif nclust <= 9:
		nr = 3
	elif nclust <= 16:
		nr = 4
	
	## initialize nr x nr subplots and fill in with that row of cluster center array (using specified columns)
	figure, axis = plt.subplots(nr,nr)
	r = 0
	while r<nr:
		c = 0
		while c<nr:
			clust = (nr*r)+c
			if clust >= nclust:
				break
			axis[r,c].plot(range(0,nhbin), clusters[clust, cols])
			axis[r,c].set_title(f'cluster {clust}')
			axis[r,c].text(nhbin, max(clusters[clust, cols]), counts.loc[clust,1], ha='right', va='top')
			c+=1
		r+=1
	plt.suptitle('scaled cluster means')
	plt.tight_layout()
	figure.add_subplot(111, frame_on = False)
	plt.tick_params(labelcolor="none", bottom = False, left = False)
	plt.xlabel('contact distance (bp)') 
	plt.ylabel('min/max scaled counts')
	plt.savefig(f'{celltype}-{orient}-{nclust}-clustercenters.png', bbox_inches='tight')
	plt.close()

## function to plot PCA
def pcaplot(transformed, colorby):
	figure, axis = plt.subplots()
	plot = axis.scatter(transformed[:,0], transformed[:,1], c = colorby[0], s = 5)
	plt.xlabel('PC 1')
	plt.ylabel('PC 2')
	axis.legend(handles = plot.legend_elements()[0], labels = colorby[1])
	plt.savefig(f'{celltype}-{orient}-{nclust}-PCA.png', bbox_inches='tight')
	plt.close()

########## SET INPUT ARGUMENTS ##########
opts = OptionParser()
usage = "usage: %prog [options] [inputs]"
opts = OptionParser(usage=usage)
opts.add_option("-c", help = "<celltype> Name of cell type as in csv names") #input files should be [celltype]-[orient]dist.csv
opts.add_option("-k", help = "<nclusters> Number of clusters to use")
opts.add_option("--nhbin", default = "150", help = "number of histogram bins. default = 150")
opts.add_option("--ngbin", default = "4053", help = "number of genomic bins, default = 4053")
options, arguments = opts.parse_args()

########## MAIN CODE ##########
## read in arguments
celltype = options.c
print(f'cell type: {celltype}')
nclust = int(options.k)
nhbin = int(options.nhbin)
ngbin = int(options.ngbin)

## initialize combined data arrays
alldatav = np.zeros((0,nhbin))
alldatavwid = np.zeros((0,nhbin+2))
alldatah = np.zeros((ngbin, 0))

## loop through orientations
orients = ['tandemminus', 'tandemplus', 'in', 'out']
for i in range(0,4):
	## read in data and initialize data arrays
	orient = orients[i]
	datafile = f'{celltype}-{orient}dist.csv'
	data = np.loadtxt(datafile, dtype = int, delimiter=',')
	scaled_data = np.zeros((data.shape))
	datawid = np.zeros((data.shape[0], data.shape[1]+2))
	## min-max normalize each row (genomic bin), add bin number as extra column to array w ID 
	for r in range(0,len(data)):
		rowr = data[r]
		if max(rowr) == 0:
			continue
		scaled_data[r] = (rowr - min(rowr)) / (max(rowr) - min(rowr))
		datawid[r] = np.append(scaled_data[r], [i, r]) #add orientation as nhbin-th column, row number as nhbin+1st column
	## add scaled data to combined data sets
	alldatav = np.append(alldatav, scaled_data, 0)
	alldatavwid = np.append(alldatavwid, datawid, 0)
	alldatah = np.append(alldatah, scaled_data, 1)

##### OPTIMIZING NUMBER OF CLUSTERS - comment out if already optimized #####
#	optimization()

## orientation-specific clustering and PCA
	print(f'running clustering on {orient}')
	(datawids, clusters, counts, transformed) = clustering(scaled_data, datawid)
	cols = slice(0, clusters.shape[1],1)
	clusterplot(clusters, cols, counts)
	colorby = (datawids[:,-1], range(0,nclust,1)) #color by last annotation, which is cluster
	pcaplot(transformed, colorby)

## vertically combined clustering and PCA
print('running clustering on vert combined')
orient = "allvert"
(alldatavwids, clusters, counts, transformed) = clustering(alldatav, alldatavwid)
cols = slice(0, clusters.shape[1],1)
clusterplot(clusters, cols, counts)

## plot PCA to check how orientations cluster
colorby = (alldatavwids[:,nhbin], ['tandem (entry-entry)', 'tandem (exit-exit)', 'inward', 'outward']) #color by first annotation, which is ligation orientation where 0=tm, 1=tp, 2=i, 3=0
pcaplot(transformed, colorby)

## horizontally combined clustering and PCA
print('running clustering on ho combined')
orient = "allho"
(alldatahwids, clusters, counts, transformed) = clustering(alldatah, alldatah)

for i in range(0,4):
	orient = f"all-{orients[i]}"
	start = nhbin*i
	end = nhbin + (nhbin*i)
	cols = slice(start, end, 1)
	clusterplot(clusters, cols, counts)

colorby = (alldatahwids[:,-1], range(0,nclust,1)) #color by last annotation, which is cluster
pcaplot(transformed, colorby)
