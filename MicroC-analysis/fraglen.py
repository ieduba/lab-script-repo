import re
from optparse import OptionParser
from matplotlib import pyplot as plt

opts = OptionParser()
opts.add_option("-s", help = "<sam> sam alignment file")
opts.add_option("-o", help = "<out> txt file output name")
opts.add_option("-p", help = "<plot> plot output file name")
options,arguments = opts.parse_args()

print('calculating fragment lengths')

seqlen = []
with open(options.s, "r") as sam, open(options.o, 'w') as outtxt:
	line = sam.readline()
	while line:
		nbases = []
		if line[0] == '@':
			line = sam.readline()
			continue
		linecol = line.split('\t')
		cigar = linecol[5]
		matches = re.findall(r'(\d*)M', cigar) #number of matches 
		dels = re.findall(r'(\d*)D', cigar) #number of deletions
		subs = re.findall(r'(\d*)X', cigar) #number of mismatches
		nbases = matches + dels + subs
		intbases = [int(i) for i in nbases]
		totalNbases = sum(intbases) 
		seqlen.append(totalNbases)
		outtxt.write(f'{totalNbases}\n')
		line = sam.readline()

print('plotting histogram')
plt.hist(seqlen,20)
plt.savefig(options.p, format = 'pdf', bbox_inches = 'tight') 
