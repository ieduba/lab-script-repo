import numpy as np
from matplotlib import pyplot as plt
from optparse import OptionParser
from scipy.signal import savgol_filter, find_peaks

plt.rcParams.update({'font.size': 12})

usage = "usage: %prog [options] [inputs]"
opts = OptionParser(usage=usage)
opts.add_option("-l", help = "<list> comma separated list of txt files of fragment lengths")
opts.add_option("-b", default = 5, help = "<binsize> size of histogram bins in bp, default = 5")
opts.add_option("-o", help = "<output> output file prefix")
options, arguments = opts.parse_args()

def cpplot(histdata, names):
	print('plotting curves')
	for i in range(len(histdata)):
		plt.plot(hbins, histdata[i], label=names[i])
	plt.xlabel('Fragment length (bp)')
	plt.ylabel('Normalized fragment counts')
	plt.xlim(0,1000)
	plt.yscale('log')
	plt.legend()
	plt.savefig(f'{options.o}.pdf', format = 'pdf', bbox_inches = 'tight')
	plt.close()

lenfiles = options.l.split(",")
print(f'loading files: {lenfiles}')
lens = [np.loadtxt(lenfile) for lenfile in lenfiles]

hbin = options.b
nhbin = int(1000 / hbin)
hbins = np.linspace(0, 1000, nhbin)

histdata = [np.histogram(samplelen, bins = nhbin, range = (0,1000), density = True)[0] for samplelen in lens]
names = [lenfile.split(".")[0] for lenfile in lenfiles]
cpplot(histdata, names)
