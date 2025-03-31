from optparse import OptionParser

opts = OptionParser()
opts.add_option('-l', help = '<loop file> path to loop file')
options,arguments = opts.parse_args()

input_file = options.l
output = input_file.rstrip('.bed')
output_file = f'{output}.longrange'
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
	for line in infile:
		c1,s1,e1,c2,s2,e2,count,exp_dv,exp_vv,exp_hv,exp_llv,exp_dq,exp_vq,exp_hq,exp_llq,region,cs1,cs2,c_label,c_size,r1,r2 = line.strip().split(',')
		# make new format, with zero base coordinates and intensity variable
		line1 = f"{c1}\t{int(s1)-1}\t{int(e1)-1}\t{c2}:{int(s2)-1}-{int(e2)-1},{count}"
		line2 = f"{c2}\t{int(s2)-1}\t{int(e2)-1}\t{c1}:{int(s1)-1}-{int(e1)-1},{count}"
		outfile.write(line1 + "\n")
		outfile.write(line2 + "\n")
