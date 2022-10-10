import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from optparse import OptionParser
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from matplotlib import cm
from scipy import stats
from scipy import interpolate
from scipy.signal import savgol_filter, argrelextrema, find_peaks
from matplotlib.axes._axes import _log as matplotlib_axes_logger

plt.rcParams.update({'font.size': 22})

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
	plt.savefig(f'{celltype}-{orient}-SSE.pdf', format = 'pdf', bbox_inches='tight')
	plt.clf()	

	## plot silhouette score vs k to look for high score
	plt.style.use('fivethirtyeight')
	plt.plot(range(2,16), silhouette)
	plt.xticks(range(2,16))
	plt.xlabel('Number of clusters')
	plt.ylabel('Silhouette coefficient')
	plt.savefig(f'{celltype}-{orient}-silhouette.pdf', format = 'pdf', bbox_inches='tight')
	plt.clf()

## function to run kmeans clustering and PCA, saving relevant results and plots
def clustering(scaled, anno):
	## run kmeans with 50 kmeans++ initializations and 500 max iterations to convergence
	kmeans = KMeans(init = 'k-means++', n_clusters = nclust, n_init = 50, max_iter = 500)
	kmeans.fit(scaled)

	## add cluster labels to array w IDs, aggregate by cluster to count # bins in each
	ids = kmeans.labels_.reshape((kmeans.labels_.shape[0],1))
	anno = np.append(anno, ids, 1) #also add cluster labels to annotation matrix
	datawids = np.append(scaled, ids, 1) #so now cluster is last column
	df = pd.DataFrame(datawids)
	countsdf = df.groupby(scaled.shape[1]).agg('count') 

	## make new df for each cluster, save as list indexed by cluster number
	clustdfs = []
	for clust in range(nclust):
		clustdf = df[df[scaled.shape[1]]==clust]
		noid = clustdf.drop(scaled.shape[1], 'columns')
		clustdfs.append(noid)
	
	## filter data by removing noisy cluster
	noisei = countsdf.index[countsdf[1]==min(countsdf[1])].tolist()
	nonoisedf = df[df[df.shape[1]-1]!=noisei[0]]
	nonoise = nonoisedf.to_numpy()[:,0:scaled.shape[1]]
	
	annodf = pd.DataFrame(anno)
	nonoiseannodf = annodf[annodf[anno.shape[1]-1]!=noisei[0]]
	nonoiseanno = nonoiseannodf.to_numpy()

	## save cluster data 
	countsdf.to_csv(f'{celltype}-{orient}-{nclust}-clustersizes.csv', columns = (0,1))
	np.savetxt(f'{celltype}-{orient}-{nclust}-clustercenters.csv', kmeans.cluster_centers_, delimiter=',')
	np.savetxt(f'{celltype}-{orient}-{nclust}-binclusters.csv', datawids, delimiter=',')
	qc = [f'min SSE: {kmeans.inertia_}', f'number of iterations: {kmeans.n_iter_}', f'cluster means shape: {kmeans.cluster_centers_.shape}']
	with open(f'{celltype}-{orient}-{nclust}-clusterQC.txt', 'w') as f:
		for line in qc:
			f.write(line)
			f.write('\n')
        ## return annotations, cluster centers, cluster counts, and cluster dfs for plotting
	return(anno, kmeans.cluster_centers_, countsdf, clustdfs, nonoise, nonoiseanno)

def runpca(scaled):
        ## run PCA (to get PCs, explained variance, and means) and flattened PCA (for plotting)
	pca = PCA()
	pcadata = pca.fit(scaled)
	transformed = pca.fit_transform(scaled)

#	np.savetxt(f'{celltype}-{orient}-components.csv', pcadata.components_, delimiter = ',')
#	np.savetxt(f'{celltype}-{orient}-explainedvariance.csv', pcadata.explained_variance_, delimiter = ',')
#	np.savetxt(f'{celltype}-{orient}-PCAtransform.csv', transformed, delimiter = ',')
#	np.savetxt(f'{celltype}-{orient}-PCAmeans.csv', pcadata.mean_, delimiter = ',')	

	## return flattened PCA for plotting
	return(transformed)

## function to smooth contact prob curves and get info about peaks
def peaks(data):
	hbin = 1500/nhbin
	smoothed = savgol_filter(data, 7, 3) #7 determined by trial and error - consider trying different levels of smoothing
	interp = interpolate.interp1d(hbins, smoothed, 'cubic')
	newbins = np.linspace(0,1500-hbin,(nhbin*5))
	smoothint = interp(newbins)
#	maxi = argrelextrema(smoothed, np.greater)
	maxi, _ = find_peaks(smoothint, distance = 150/hbin)
#	maxi = np.arange(0,smoothed.shape[0])[np.r_[True, smoothed[1:] > smoothed[:-1]] & np.r_[smoothed[:-1] > smoothed[1:], True]]
	maxx = [newbins[i] for i in maxi]
	maxes = [smoothint[i] for i in maxi]
	return(smoothed, maxx, maxes)

## function to plot cluster centers (means)
def clusterplot(clusters, cols, counts, maxx, maxy, smoothed):
	hbin = 1500/nhbin
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
			axis[r,c].plot(hbins, clusters[clust, cols], zorder = 0)
			axis[r,c].plot(hbins, smoothed[clust], zorder = 1)
			axis[r,c].scatter(maxx[clust],maxy[clust], s = 10, c = 'black', zorder = 2)
			axis[r,c].set_title(f'cluster {clust}')
			axis[r,c].text(1500, max(clusters[clust, cols]), counts.loc[clust,1], ha='right', va='top')
			c+=1
		r+=1
	plt.suptitle('scaled cluster means')
	plt.tight_layout()
	figure.add_subplot(111, frame_on = False)
	plt.tick_params(labelcolor="none", bottom = False, left = False)
	plt.xlabel('contact distance') 
	plt.ylabel('scaled contact probability')
	plt.savefig(f'{celltype}-{orient}-{nclust}-clustercenters.pdf', format = 'pdf', bbox_inches='tight')
	plt.close()

	## plot all cluster averages on one plot
	matplotlib_axes_logger.setLevel('ERROR') #to stop printing dumb error about using RGB for color
	cmap = plt.get_cmap('viridis')
	colors = cmap(np.linspace(0, 1, nclust))
	for clust in range(nclust):
		if counts.loc[clust,1] == min(counts.loc[:,1]):
			continue
		plt.plot(hbins, clusters[clust,cols], c = colors[clust])
	plt.xlabel('contact distance')
	plt.ylabel('scaled contact probability')
	plt.savefig(f'{celltype}-{orient}-{nclust}-allclustercenters.pdf', format = 'pdf', bbox_inches = 'tight')
	plt.close()

	## plot all maxes on one plot
	for clust in range(nclust):
		if counts.loc[clust,1] == min(counts.loc[:,1]):
			continue
		plt.plot(maxx[clust], maxy[clust], c = colors[clust])
		plt.scatter(maxx[clust], maxy[clust], c = colors[clust])
	plt.xlabel('contact distance')
	plt.ylabel('scaled peak height')
	plt.savefig(f'{celltype}-{orient}-{nclust}-allpeaks.pdf', format = 'pdf', bbox_inches = 'tight')
	plt.close()

## function to plot PCA
def pcaplot(transformed, colorby):
	figure, axis = plt.subplots()
	plot = axis.scatter(transformed[:,0], transformed[:,1], c = colorby[0], s = 4, alpha = 0.5)
	plt.xlabel('PC 1')
	plt.ylabel('PC 2')
#	axis.legend(handles = plot.legend_elements()[0], labels = colorby[1])
	plt.savefig(f'{celltype}-{orient}-{colorby[2]}-PCA.pdf', format = 'pdf', bbox_inches='tight')
	plt.close()

## function to plot heatmap of all rows in each cluster
def clusterheatmap(dfs, cols):
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
			df = dfs[clust]
			axis[r,c].imshow(df.loc[:,cols], aspect='auto',interpolation='nearest')
			axis[r,c].set_title(f'cluster {clust}')
			c+=1
		r+=1
	plt.suptitle('all contact prob')
	plt.tight_layout()
	figure.add_subplot(111, frame_on = False)
	plt.tick_params(labelcolor="none", bottom = False, left = False)
	plt.xlabel('contact distance (bp)')
	plt.ylabel('genomic bin')
	plt.savefig(f'{celltype}-{orient}-{nclust}-heatmaps.pdf', format = 'pdf', bbox_inches='tight')
	plt.close()

## function to plot cluster means, heatmaps, and annotations
def comboplot(clusters, cols, counts, dfs, anno):
	grid = (3, nclust*5)
	plt.figure(figsize = (nclust*3,4))
	for clust in range(nclust):
		andf = pd.DataFrame(anno)
		clustanno = andf[andf[6]==clust]
		noid = clustanno.drop(6, 'columns')	
	#	annoz = np.zeros((noid.shape))
	#	for c in range(noid.shape[1]):
	#		annoz[:,c] = stats.zscore(noid.loc[:,c])
		
		annost = pd.DataFrame(noid).sort_values(0,0) #sort annos by GC	
		df = dfs[clust]
		dfst = df.reindex(index=noid.index)

		noid.to_csv(f'{celltype}-{nclust}-cluster{clust}-annotations.csv')		
		
		p1 = plt.subplot2grid(grid,(0,clust*5),1,3)
		p1.axes.xaxis.set_visible(False)
		p1.text(1500, max(clusters[clust, cols]), counts.loc[clust,1], ha='right', va='top')
		p1.set_title(f'cluster {clust}')
		p1.set_ylabel('contact prob')
		p2 = plt.subplot2grid(grid,(1,clust*5),2,3)
		p2.axes.yaxis.set_visible(False)
		p2.set_xlabel(f'contact distance x{hbin}')
		p3 = plt.subplot2grid(grid,(1,(clust*5)+3),2,1)
		p3.axes.yaxis.set_visible(False)
		p3.axes.xaxis.set_visible(False)
		
		p1.plot(hbins, clusters[clust,cols])
		p2.imshow(dfst.loc[:,cols], aspect = 'auto', interpolation = 'nearest')
		p3.imshow(annost, aspect = 'auto', interpolation = 'nearest')
#	plt.tight_layout()
	plt.savefig(f'{celltype}-{orient}-{nclust}-comboplot.pdf', format = 'pdf', bbox_inches = 'tight')
	plt.close()


########## SET INPUT ARGUMENTS ##########
usage = "usage: %prog [options] [inputs]"
opts = OptionParser(usage=usage)
opts.add_option("-c", help = "<celltype> Name of cell type as in csv names") #input files should be [celltype]-[orient]dist.csv
opts.add_option("-k", help = "<nclusters> Number of clusters to use")
opts.add_option("-g", help = "<gcbins> GC content in bins as contact probability")
opts.add_option("--nhbin", default = "150", help = "number of histogram bins. default = 150")
opts.add_option("--ngbin", default = "4053", help = "number of genomic bins, default = 4053")
options, arguments = opts.parse_args()

########## MAIN CODE ##########
## read in arguments
celltype = options.c
print(f'cell type: {celltype}')
nclust = int(options.k)
print(f'n clusters: {nclust}')
nhbin = int(options.nhbin)
print(f'n histogram bins: {nhbin}')
ngbin = int(options.ngbin)
print(f'n genomic bins: {ngbin}')
gc = np.loadtxt(options.g)
hbin = 1500/nhbin #size of histogram bins
hbins = np.linspace(0,1500-hbin,nhbin) #starts of histogram bins

## initialize combined data arrays
#alldatav = np.zeros((0,nhbin))
#alldatavwid = np.zeros((0,nhbin+1))
alldatah = np.zeros((ngbin, 0))
allcounts = np.zeros((ngbin))
annotations = np.zeros((ngbin,6))

## loop through orientations
orients = ['tandemminus', 'tandemplus', 'in', 'out']
for i in range(0,4):
	## read in data and initialize data arrays
	orient = orients[i]
	datafile = f'{celltype}-{orient}dist.csv'
	data = np.loadtxt(datafile, dtype = int, delimiter=',')
	scaled_data = np.zeros((data.shape))
	datawid = np.zeros((data.shape[0], data.shape[1]+1))
	orcounts = np.zeros((ngbin))
	## min-max normalize each row (genomic bin), add bin number as extra column to array w ID 
	for r in range(0,len(data)):
		rowr = data[r]
		if orient == 'in':
			rowr[0:int(250/hbin)] = 0 #trim inward <250 bp to get rid of self ligation
		if max(rowr) == 0:
			continue
		scaled_data[r] = (rowr - min(rowr)) / (max(rowr) - min(rowr))
		smoothed, maxi, maxes = peaks(scaled_data[r]) #get peak maxes for this bin and orientation
		annotations[r,2+i] = maxes[0]-maxes[1] #save y distance from first to second peak as annotation
		orcounts[r] = sum(rowr) 
#		datawid[r] = np.append(scaled_data[r], i) #add orientation as nhbinth column
	## add scaled data to combined data sets
#	alldatav = np.append(alldatav, scaled_data, 0)
#	alldatavwid = np.append(alldatavwid, datawid, 0)
	alldatah = np.append(alldatah, scaled_data, 1)
	allcounts = np.add(allcounts,orcounts)
	
annotations[:,0] = gc
annotations[:,1] = allcounts

##### OPTIMIZING NUMBER OF CLUSTERS - comment out if already optimized #####
#	optimization()

##### RUNNING CLUSTERING AND PCA #####
## vertically combined PCA to see how orientations cluster
#transformed = runpca(alldatav)
#colorby = (alldatavwids[:,nhbin], ['tandem (entry-entry)', 'tandem (exit-exit)', 'inward', 'outward'], 'orientation') #color by first annotation, which is ligation orientation where 0=tm, 1=tp, 2=i, 3=0
#pcaplot(transformed, colorby)

## horizontally combined clustering and PCA
print('running clustering on ho combined')
orient = "allho"

(anno, clusters, counts, clustdfs, nonoise, nonoiseanno) = clustering(alldatah, annotations)
transformed = runpca(nonoise)
transformedall = runpca(alldatah)

#annotations: col 0 = gc, col 1 = fragment count, col 2-5 = slopes between 1st and 2nd peaks (2 = tm, 3 = tp, 4 = in, 5 = out), col 6 = cluster assignment 
colorby = (anno[:,6], np.unique(anno[:,6]).tolist(), 'cluster-all')
pcaplot(transformedall, colorby)
colorby = (nonoiseanno[:,6], np.unique(nonoiseanno[:,6]).tolist(), 'cluster')
pcaplot(transformed, colorby)

stallcounts = np.sort(allcounts)
countlabels = []
for i in np.linspace(0, len(allcounts)-1, 8, dtype = int):
        countlabels.append(str(stallcounts[i]))
colorby = (nonoiseanno[:,1], countlabels, 'counts')
pcaplot(transformed, colorby)

gcst = np.sort(gc)
gclabels = []
for i in np.linspace(0, len(gc)-1, 8, dtype=int):
        gclabels.append(str(gcst[i]))
colorby = (nonoiseanno[:,0], gclabels, 'gc')
pcaplot(transformed,colorby)

stslope = np.sort(nonoiseanno[:,5])
slopelabels = []
for i in np.linspace(0, nonoiseanno.shape[0]-1, 8, dtype = int):
	slopelabels.append(str(stslope[i]))
colorby = (nonoiseanno[:,5], slopelabels, 'slope_1-2')
pcaplot(transformed,colorby)

for i in range(4):
	orient = f"all-{orients[i]}"
	start = nhbin*i
	end = nhbin + (nhbin*i)
	transformed = runpca(nonoise[:,start:end])
	#transformed = runpca(alldatah[:,start:end])

	maxx = [[] for j in range(nclust)]
	maxy = [[] for j in range(nclust)]
	smoothed = [[] for j in range(nclust)]
	for j in range(nclust):
		line = clusters[j, start:end]
		smoothed[j], maxx[j], maxy[j] = peaks(line)
		with open(f'{celltype}-{orient}-clust{j}-peakinfo.txt', 'w') as f:
			f.write('X \t Y \n')
			for k in range(len(maxy[j])):
				f.write(f'{str(maxx[j][k])}\t{str(maxy[j][k])}\n')
			f.write(f'\n NRL: {np.average(np.diff(maxx[j]))}')

	## plot cluster contact prob curves
	cols = slice(start, end, 1)	
	clusterplot(clusters, cols, counts, maxx, maxy, smoothed)
	clusterheatmap(clustdfs,cols)
	comboplot(clusters, cols, counts, clustdfs, anno)

	## plot individual PCAs with various color schemes
#	colorby = (nonoiseanno[:,6], np.unique(anno[:,6]).tolist(),'cluster')
#	pcaplot(transformed,colorby)
#	colorby = (nonoiseanno[:,1], countlabels, 'counts')
#	pcaplot(transformed, colorby)
#	colorby = (nonoiseanno[:,0], gclabels, 'gc')
#	pcaplot(transformed, colorby)


