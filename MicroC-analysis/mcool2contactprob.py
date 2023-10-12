import numpy as np
import cooler
import cooltools
from matplotlib import pyplot as plt
from optparse import OptionParser

opts = OptionParser()
opts.add_option("-m", help = "<mcool> mcool file with multiple resolutions")
opts.add_option("-r", help = "<resolution> bp resolution to make contact prob")
opts.add_option("-b", help = "<balance> True or False, perform ICE balancing on contact matrix")
opts.add_option("--chroms", help = "<chrom.sizes> chrom sizes file for your genome, default: hg38", default = '/lustre/fs4/risc_lab/store/risc_data/downloaded/hg38/genome/chrom.sizes')
opts.add_option("--chunksize", help = "size of submatrices in which to chunk contact matrix to calculate P(s), in bp. Default = 100,000", default = 100000)
opts.add_option("--plotmax", help = "max bp contact for P(s) curve. Must be multiple of resolution and smaller than chunk size. Default = 1500", default = 1500)
options, arguments = opts.parse_args()

infile = options.m
resolution = int(options.r) 
bal = eval(options.b)
csraw = np.loadtxt(options.chroms, dtype='str')
chunksize = int(options.chunksize)
plotmax = int(options.plotmax)

## load in mcool file at chosen resolution
cooldata = cooler.Cooler(f'{infile}::resolutions/{resolution}')
csdict = {csraw[n][0]:int(csraw[n][1]) for n in range(len(csraw))}

## set name for this contact prob
sp = options.m.split('.')
if bal == True:
	btxt = 'bal'
else:
	btxt = 'nobal'
name = f'{sp[0]}_{options.r}_{options.chunksize}_{btxt}contactprob'

## iterate through chromosomes and extract balanced contact matrix for that chrom from cooler data
chlist = ['chr'+str(n) for n in range(1,22)]
chlist.append('chrX')
allprobs = {}
for ch in chlist:
	print(ch)

	cs = csdict[ch]
	chunks = list(range(0,cs,chunksize)) #load in limited bp at a time to not overload memory
	chunks.append(cs)
	
	chprobs = {}
	for c in range(len(chunks)-1):
  
		chmat = cooldata.matrix(balance=bal).fetch(f'{ch}:{chunks[c]}-{chunks[c+1]}')
#		print(f'chunk: {chunks[c]}-{chunks[c+1]}')	
	
		## iterate through contact matrix starting at central diagonal (ie j=i) and moving
		## outward (ie j=i+1, j=i+2 ... j=i+d) to build dictionary of contact probability as 
		## a function of inter-locus distance d 
		probs = {}
		for i in range(len(chmat)):
			d = 0
			j = i + d
			while j < len(chmat[0]):
				if np.isnan(chmat[i][j]):
					val = 0
				else: 
					val = chmat[i][j]
				if d in probs.keys():
					probs[d].append(val)
				else:
					probs[d] = [val]
				d += 1
				j = i + d

		## make dict of mean contact probs for this chunk and add to this chromosome's dict
		## of contact probs
		meanprobs = {d:sum(probs[d])/len(probs[d]) for d in probs}
		for d in range(len(meanprobs)):
			if d in chprobs.keys():
				chprobs[d].append(meanprobs[d])
			else:
				chprobs[d] = [meanprobs[d]]

	## make dict of mean contact probs for this chrom and add to the dict of all contact probs 
	meanchprobs = {d:sum(chprobs[d])/len(chprobs[d]) for d in chprobs}
	
	for d in range(len(meanchprobs)):
		if d in allprobs.keys():
			allprobs[d].append(meanchprobs[d])
		else:
			allprobs[d] = [meanchprobs[d]]
	
allmeanprobs = [sum(allprobs[d])/len(allprobs[d]) for d in allprobs]

## save as csv with name based on cooler file name

np.savetxt(f'{name}.csv', allmeanprobs, delimiter=',', fmt='%.0f')

## plot contact probability curve
plt.plot(range(2*resolution,plotmax,resolution),allmeanprobs[2:int(plotmax/resolution)])
plt.savefig(f'{name}.pdf', format = 'pdf', bbox_inches = 'tight')
plt.clf()
