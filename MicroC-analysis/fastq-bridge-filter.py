import re
from optparse import OptionParser
from itertools import islice

opts = OptionParser()
opts.add_option('--r1', help = '<read1> unzipped R1 fastq')
opts.add_option('--r2', help = '<read2> unzipped R2 fastq')
options, arguments = opts.parse_args()

r1filt = re.sub('.fastq', '-bridgefilt.fastq', options.r1)
r2filt = re.sub('.fastq', '-bridgefilt.fastq', options.r2)

bridgestart = 'AGGTTCGT'
bridgeend = 'ACGAACCT'

with open(options.r1,'r') as r1, open(options.r2,'r') as r2, open(r1filt,'w') as r1out, open(r2filt,'w') as r2out:
	lines1 = list(islice(r1,4))
	lines2 = list(islice(r2,4))
	while lines1:
		seq1 = lines1[1]
		seq2 =  lines2[1]
		if bridgestart in seq1 or bridgeend in seq1 or bridgestart in seq2 or bridgeend in seq2:
			for line1 in lines1:
				r1out.write(line1)
			for line2 in lines2:
				r2out.write(line2)
		lines1 = list(islice(r1,4))
		lines2 = list(islice(r2,4))
