import numpy as np
import random
from optparse import OptionParser
from matplotlib import pyplot as plt

plt.rcParams.update({'font.size': 18})

opts = OptionParser()
opts.add_option("-f", help = "<filt> alignment sizes from filtered reads")
opts.add_option("-u", help = "<unfilt> alignment sizes from unfiltered reads")
opts.add_option("-p", help = "<plot> plot output file name")
options,arguments = opts.parse_args()

seqlenfilt = []
filtsum = 0
with open(options.f, 'r') as filt:
	line = filt.readline()
	while line:
		seq = int(line.rstrip())
		seqlenfilt.append(seq)
		filtsum += seq
		line = filt.readline()

seqlenunfilt = []
unfiltsum = 0
with open(options.u, 'r') as unfilt:
	line = unfilt.readline()
	while line:
		seq = int(line.rstrip())
		seqlenunfilt.append(seq)
		unfiltsum += seq
		line = unfilt.readline()

#plt.boxplot([seqlenfilt,seqlenunfilt])
#for i in random.sample(range(len(seqlenfilt)), 1000):
#	y = seqlenfilt[i]
#	x = np.random.normal(1, 0.04, 1)
#	plt.scatter(x,y,s=0.5,c='Blue')
#for j in random.sample(range(len(seqlenunfilt)), 1000):
#	y = seqlenunfilt[j]
#	x = np.random.normal(2, 0.04, 1)
#	plt.scatter(x,y,s=0.5,c='Blue')
#plt.savefig(options.p, format = 'pdf', bbox_inches = 'tight') 

f_w = np.empty(np.asarray(seqlenfilt).shape)
u_w = np.empty(np.asarray(seqlenunfilt).shape)
f_w.fill(1/np.asarray(seqlenfilt).shape[0])
u_w.fill(1/np.asarray(seqlenunfilt).shape[0])

plt.hist([seqlenfilt, seqlenunfilt], np.linspace(0,150,30), weights = [f_w, u_w], color=['coral', 'darkturquoise'], label = ['euchromatin','heterochromatin'])
plt.legend(loc='upper left')
plt.xlabel('alignment length (bp)')
plt.ylabel('density')
plt.savefig(options.p, format = 'pdf', bbox_inches = 'tight') 
