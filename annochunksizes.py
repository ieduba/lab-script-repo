import numpy as np
from matplotlib import pyplot as plt
from optparse import OptionParser

usage = "usage: %prog [options] [inputs]"
opts = OptionParser(usage=usage)
opts.add_option("-b", help = "<bedfilesizes> 1D list of sizes of bed file chunks")
options, arguments = opts.parse_args()

bed = np.loadtxt(options.b)
print('total: ', sum(bed))
print('mean: ', np.average(bed))
print('max: ',max(bed))
plt.hist(bed, bins = np.linspace(0,1000000,50))
plt.xlabel('annotation chunk size')
plt.ylabel('counts')
plt.suptitle(options.b)
plt.savefig(f'{options.b}.pdf', format='pdf')
