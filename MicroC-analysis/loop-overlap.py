from optparse import OptionParser

opts = OptionParser()
opts.add_option("--b1", help = "<bed1> Path to first bed file") 
opts.add_option("--b2", help = "<bed2> Path to second bed file")
opts.add_option("--op", help = "<out-prefix> Name of comparison for output file prefixes")
options, arguments = opts.parse_args()

bed1 = open(options.b1)
bed1lines = open(options.b1).read()
bed2 = open(options.b2)
bed2lines = open(options.b2).read()
only1 = open(f'{options.op}-1-only-loops.bed', 'w')
only2 = open(f'{options.op}-2-only-loops.bed', 'w')
common = open(f'{options.op}-common-loops.bed', 'w')

for line in bed1:
	if line in bed2lines:
		common.write(line)
	else:
		only1.write(line)

for line in bed2:
	if line not in bed1lines:
		only2.write(line)

bed1.close()
bed2.close()
only1.close()
only2.close()
common.close()	
