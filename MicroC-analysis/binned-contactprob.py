import numpy as np
import time
import math
import subprocess
from optparse import OptionParser
from multiprocessing import Pool

## set command line arguments
opts = OptionParser()
usage = "usage: %prog [options] [inputs]"
opts = OptionParser(usage=usage)
opts.add_option("-b", help = "<bed> Bed file of genomic bins")
opts.add_option("-p", help = "<pairs> .pairs file of contact pairs")
opts.add_option("-c", help = "<cores> Number of cores to run on")
opts.add_option("--maxdist", default = "1500", help = "max contact distance to count (bp), default=1500")
opts.add_option("--histbin", default = "10", help = "size of histogram bin (bp), default=10")
options, arguments = opts.parse_args()

## define function to assign pairs to bins and make array of contact distance counts in each orientation
def pairsbin(start):
	## open pairs file to read in line by line
	pairs = open(options.p, "r")
	
	## initialize empty counts arrays
	end = min(start+chunksize, len(bins))
	thischunksize = end - start
#	tandemdist = np.zeros((thischunksize,nhistbins))
	tandempdist = np.zeros((thischunksize,nhistbins))
	tandemmdist = np.zeros((thischunksize,nhistbins))
	indist = np.zeros((thischunksize,nhistbins))
	outdist = np.zeros((thischunksize,nhistbins))

	## loop through lines of pairs file 
	line = pairs.readline()
	while line:
		## check if line is header, if so, move to next line
		if line[0] == '#':
			line = pairs.readline()
			continue

                ## split line into columns and calculate contact distance, move to next line if greater than max distance or interchrom contact
		linecol = line.split('\t')
		pos1= int(linecol[2])
		pos2= int(linecol[4])
		dist = pos2 - pos1
		if (dist >= maxcontact or linecol[1] != linecol[3]):
			line = pairs.readline()
			continue

		## loop through lines of bed file to test if pairs line falls in that bin. break loop when match is found. move to the next pairs line if no match is found
		gbin = "NA"
		for b in range(start, end):
			if (linecol[1]==bins[b][0] and bins[b][1] <= pos1 < bins[b][2]):
				gbin = b - start
				break
		if gbin == "NA":
			line = pairs.readline()
			continue

		## assign the contact to a genomic bin and distance bin in the appropriate array and increase the count in that cell
		distbin = math.floor(dist/histbinsize)
		if (linecol[5]=='+' and linecol[6]=='-'):
			indist[gbin,distbin] += 1
		if (linecol[5]=='-' and linecol[6]=='+'):
			outdist[gbin,distbin] += 1
		if (linecol[5] == '+' and linecol[6] == '+'):
			tandempdist[gbin,distbin] += 1
		if (linecol[5] == '-' and linecol[6] == '-'):
			tandemmdist[gbin,distbin] += 1
#		else:
#			tandemdist[gbin,distbin] += 1

		## iterate to next line in pairs file
		line = pairs.readline()
	
	## close pairs file	
	pairs.close()
#	np.savetxt(f'tandemdist-{start}.csv', tandemdist, delimiter=',', fmt='%.0f')
#	np.savetxt(f'indist-{start}.csv', indist, delimiter=',', fmt='%.0f')
#	np.savetxt(f'outdist-{start}.csv', outdist, delimiter=',', fmt='%.0f')

	## return tuple containing all three counts arrays
	return (tandemmdist, tandempdist, indist, outdist)


#### MAIN CODE ####

## start timer
tic = time.perf_counter()

## print input file names
print(f"bed file: {options.b}")
print(f"pairs file: {options.p}")

## read in command line arguments
bins = np.genfromtxt(options.b, dtype=None, encoding=None)
pairs = open(options.p, "r")
ncores = int(options.c)
maxcontact = int(options.maxdist)
histbinsize = int(options.histbin)

## calculate histogram bin size for counting and chunk size/start coordinates for parallelization 
nhistbins = int(maxcontact/histbinsize)
print(f"nhistbins: {nhistbins}")
chunksize = int(len(bins)/(ncores-1))
print(f"chunksize: {chunksize}")
starts = range(0, len(bins), chunksize)
print(f"nchunks: {len(starts)}")

## split bed file and run parallelized pairsbin function to assign pairs to bins and count contact distances
print("making counts submatrices")
pool = Pool(processes = ncores)
output = pool.map(pairsbin, starts)
print("counts submatrices made, combining")

## combine counts matrices
alltandemm = np.empty((len(bins), nhistbins))
alltandemp = np.empty((len(bins), nhistbins))
allin = np.empty((len(bins), nhistbins))
allout = np.empty((len(bins), nhistbins))
for s in range(len(starts)):
	chunkstart = starts[s]
	chunkend = min(chunkstart+chunksize,len(bins))
	alltandemm[chunkstart:chunkend,:] = output[s][0]
	alltandemp[chunkstart:chunkend,:] = output[s][1]
	allin[chunkstart:chunkend,:] = output[s][2]
	allout[chunkstart:chunkend,:] = output[s][3]

## parse the .pairs name to get name for output
pairname = options.p
sp = pairname.split('-')
tu = ((sp[0], sp[1]))
csvname = '-'.join(tu)

## save contact distance counts matrices as csv
np.savetxt(f'{csvname}-tandemminusdist.csv', alltandemm, delimiter=',', fmt='%.0f')
np.savetxt(f'{csvname}-tandemplusdist.csv', alltandemp, delimiter=',', fmt='%.0f')
np.savetxt(f'{csvname}-indist.csv', allin, delimiter=',', fmt='%.0f')
np.savetxt(f'{csvname}-outdist.csv', allout, delimiter=',', fmt='%.0f')

## stop timer and print run time
toc = time.perf_counter()
print(f"Run time {toc - tic:0.4f} seconds")
