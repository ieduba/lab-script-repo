import numpy as np
from matplotlib import pyplot as plt
from optparse import OptionParser
from scipy.signal import savgol_filter, find_peaks

plt.rcParams.update({'font.size': 22})

usage = "usage: %prog [options] [inputs]"
opts = OptionParser(usage=usage)
opts.add_option("-f", help = "<file> txt file with list of contact distances")
opts.add_option("--binsize", default = 30, help = "<binsize> size of histogram bins in bp, default = 30")
opts.add_option("--f2", default = None, help = "<file2> second txt file to compare, optional")
opts.add_option("--f3", default = None, help = "<file3> third txt file to compare, optional")
opts.add_option("--f4", default = None, help = "<file4> fourth txt file to compare, optional")
options, arguments = opts.parse_args()

def peaks(histdata, sample):
	print('calculating peaks')
	smoothed = savgol_filter(histdata, 7, 3)
	maxi, _ = find_peaks(smoothed, distance = 150/hbin)
	maxes = [smoothed[i] for i in maxi]
	nrl = np.average(np.diff(maxi))*hbin
	print(f'saving NRL for {sample}')
	with open(f'{sample}-peakinfo.txt', 'w') as f:
		f.write('X \t Y \n')
		for i in range(len(maxes)):
			f.write(f'{str(maxi[i]*hbin)}\t{str(maxes[i])}\n')
		f.write(f'\n NRL: {nrl}')
	return(smoothed, maxi, maxes)

def cpplot(histdata, maxi, maxes, smoothed):
	print('plotting curves')
	plt.plot(hbins, histdata)
	i = 0
	histdata=[]
	while extras[i] != None:
		distse = np.loadtxt(extras[i])
		histdatae, _ = np.histogram(distse, bins = nhbin, range = (0,1500), density = True)
		histdata.append(histdatae)
		plt.plot(hbins, histdata[i])
		i+=1
	plt.xlabel('contact distance')
	plt.ylabel('density')
	plt.xlim(0,1500)
	plt.ylim(bottom = 0)
	plt.savefig(f'{name}-cpp.pdf', format = 'pdf', bbox_inches = 'tight')
	plt.close()

	print('plotting maxes')
	plt.plot(maxi*hbin, maxes)
	plt.scatter(maxi*hbin, maxes)
	i = 0
	while extras[i] != None:
		sample = "-".join(extras[i].split("-")[0:7]).split("/")[1]
		smoothede, maxie, maxese, = peaks(histdata[i], sample)
		plt.plot(maxie*hbin, maxese)
		plt.scatter(maxie*hbin, maxese)
		i+=1
	plt.xlabel('contact distance')
	plt.ylabel('peak height')
	plt.xlim(0,1500)
	plt.ylim(0.0004,0.0014)
	plt.savefig(f'{name}-peaks.pdf', format = 'pdf',  bbox_inches = 'tight')
	plt.close()	

print(f'loading file: {options.f}')
dists = np.loadtxt(options.f)
hbin = options.binsize
extras = [options.f2, options.f3, options.f4, None]
annos = [options.f.split('-')[0].split('/')[1]]
#annos = [options.f.split('_')[1]]
i = 0
while extras[i] != None:
	annos.append(extras[i].split('-')[0].split('/')[1])
	#annos.append(extras[i].split('_')[1])
	i += 1
name = f'{"-".join(annos)}-{"-".join(options.f.split(".")[0].split("-")[1:7])}'
#name = f'{"-".join(annos)}-{"-".join(options.f.split(".")[0].split("/")[1].split("-")[0:6])}'
nhbin = int(1500 / hbin)
hbins = np.linspace(0, 1500, nhbin)
histdata, _ = np.histogram(dists, bins = nhbin, range = (0,1500), density = True)
sample = "-".join(options.f.split("-")[0:7]).split("/")[1]
smoothed, maxi, maxes = peaks(histdata, sample)
cpplot(histdata, maxi, maxes, smoothed)
