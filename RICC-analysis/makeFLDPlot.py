#!/usr/bin/env python

# makeFLDplot.py

# Author: Viviana Risca, Stanford University 
# modified from script by Jason Buenrostro and Alicia Schep
# May 2015
#
# Will make an insert size distribution plot from bed regions 
# (not just symmetric about centers) without adding Tn5 offsets
# Note: To avoid double-counting reads, this code counts only reads mapping
#       to the plus strand, and within a region, a read is counted only if 
#       its left coordinate falls within the region. 
# June 22, 2015 - Changed page size for figures
# July 12, 2015 - Version 2.0 - added support for user-defined number of columns
# September 11, 2015 - Version 3.0 Also added to V1.0 and V2.0: changed read mapping so left fragment left position is explicitly in the interval. 
# September 12, 2015 - Version 3.0 Added column of window names at 0 position

##### IMPORT MODULES #####
# import necessary for python
import os
import sys
import numpy as np
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
from optparse import OptionParser
import random

#### OPTIONS ####
# read options from command line
opts = OptionParser()
usage = "usage: %prog [options] [inputs]"
opts = OptionParser(usage=usage)
opts.add_option("-a", help="<Reads> Accepts sorted BAM file")
opts.add_option("-b", help="<Bed> Accepts bed file")
opts.add_option("-o", help="OutputFile")
opts.add_option("-c", default="20", help="number of threads to use, default=20")
opts.add_option("-l", default="700", help="maximum fragment length (in bp) to count, default=700")
opts.add_option("-v", action="store_true", default=False, help="Print insert size across all regions")
opts.add_option("-u", action="store_true", default=False, help="Print uncompressed output file")
opts.add_option("-n", action="store_true", default=False, help="With -u, print region names as first column")
opts.add_option("--window",default='20',help="window size for ploting")
options, arguments = opts.parse_args()

# return usage information if no argvs given
if len(sys.argv)==1:
    os.system(sys.argv[0]+" --help")
    sys.exit()

##### DEFINE FUNCTIONS #####
# assign mat
def asn_mat(l_pos,r_pos,mat,s_int,e_int,t,i,starti,weight):
    if (float(l_pos)>=s_int and float(l_pos)<e_int) and t<(cols-1):
        #print('i in asn_mat is: '+str(i)) # debug
        mat[i-starti][t] += weight
    return mat

#compute insert size dist for a particular chunk of BED entries
def sub_Mat(start):
    # get chunk end coordinate (note, it's exclusive end, so that coord is actually coord of next start)
    end = min(start+chunksize,len(p1_ints))
    thischunksize = end-start
    # print('thischunksize: '+str(thischunksize))
    # print('start: '+str(start))
    # print('end: '+str(end))

    # initialize data matrix for the chunk
    mat = np.zeros([thischunksize,cols])
    # loop through the intervals and get relevent info
    bamfile = pysam.Samfile(options.a, "rb")
    
    for i in range(start,end): # (actually runs from start to end-1)
        # get interval as num
        s_int=int(p1_ints[i][1])
        e_int=int(p1_ints[i][2])
        # loop through rds
        for p2_rds in bamfile.fetch(str(p1_ints[i][0]), s_int, e_int):
            #check mapping quality
            if p2_rds.mapq<30:# or p2_rds.is_proper_pair==False:
                continue
            # get read positions
            if p2_rds.is_reverse:
		continue
            else:
		l_pos = p2_rds.pos
		ilen = abs(p2_rds.tlen)
		r_pos=l_pos+ilen-1
                mat = asn_mat(l_pos,r_pos,mat,s_int,e_int,ilen,i,start,1)
    return mat

##### INPUTS AND OUTPUTS #####
# get intervals
p1_ints = np.loadtxt(options.b,'str')

##### SCRIPT #####
# open and read BAM file
#bamfile = pysam.Samfile(options.a, "rb")  # messes with file handle, add into func instead

# determine number of columns for matrix
cols = int(options.l)

# split bedfile into chunks
maxi=len(p1_ints)
chunksize=maxi/int(options.c)
# chunks=maxi/chunksize # not actually used
starts=range(0,maxi,chunksize)
rows = maxi
# print('rows: '+str(rows)) # debug

# parallel processed computation of matrix for each chunk
if __name__ == "__main__":
    pool = Pool(processes=int(options.c))
    sub_mats=pool.map(sub_Mat, starts,1)

# print('Submat assignment finished')

# place matrix of each chunk into matrix for all
mat = np.zeros([rows,cols])
for i in range(len(starts)):
    chunkstart = starts[i]
    chunkend = min(chunkstart+chunksize,maxi)
    # print('submat rows: '+str(np.size(sub_mats[i],0))) # debug
    # print('submat columns: '+str(np.size(sub_mats[i],1))) # debug
    mat[chunkstart:chunkend,:] = sub_mats[i]

# get row vector
if options.v == True:
    mat = np.sum(mat,0)

# save matrix
if not options.o:
    n1=os.path.basename(options.a)
    n2=os.path.basename(options.b)
    if options.v == True: options.o=n1+'.'+n2+'.vect'
    else: options.o=n1+'.'+n2+'.vplot'
if options.u == True:
    if options.n == True:
        intnames=np.reshape(p1_ints[:,3],(rows,1))
        # print('mat rows: '+str(np.size(mat,0))) # debug 
        # print('p1_ints rows: '+str(np.size(intnames,0))) # debug
        # print('p1_ints columns: '+str(np.size(intnames,1))) # debug
        np.savetxt(options.o,np.hstack((intnames,mat)),delimiter='\t',fmt='%s')
    else:
        if np.size(p1_ints,1) >= 4:
            intnames=np.reshape(p1_ints[:,3],(rows,1))
            np.savetxt(options.o+'.intnames',intnames,delimiter='\t',fmt='%s')
        np.savetxt(options.o,mat,delimiter='\t',fmt='%s')
else:
    np.save(options.o,mat)

# plot
fig=plt.figure(figsize=(5.0, 8.0))
#xran=min(500,int(options.e))
#yran=min(500,rows)
if options.v == True:
    plt.plot(mat,'k-')
    plt.plot(np.convolve(mat,np.ones(int(options.window)),'same')/int(options.window),'r')
    plt.xlabel('Fragment length (bp)')
    plt.ylabel('Counts over all regions')
else:
    plt.imshow(mat,origin='lower',aspect='equal')
    plt.xlabel('Fragment length (bp)')
    plt.ylabel('Region')

# save figure
fig.savefig(options.o+'.png')
plt.close(fig)
