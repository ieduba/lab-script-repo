import alluvial
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np

list1 = open(f'Scr-1_st_gene.lab.st.csv', 'r')
list2 = open(f'H1-1_st_gene.lab.st.csv', 'r')
line1 = list1.readline()
line2 = list2.readline()
quantiles = []
up = []
down = []
while line1:
	line1col = line1.split(' ')
	line2col = line2.split(' ')
	quantiles.append([f'WT-{line1col[-1][1]}', f'low-{line2col[-1][1]}'])
	q1 = int(line1col[-1][1].split()[0])
	q2 = int(line2col[-1][1].split()[0])
	if q1 > q2:
		down.append(line1col[0])
	if q1 < q2:
		up.append(line2col[0])	

	line1 = list1.readline()
	line2 = list2.readline()
ax = alluvial.plot(quantiles, disp_width = True, wdisp_sep = ' '*2, figsize=(7,5), a_sort=["WT-4","WT-3","WT-2","WT-1"], b_sort=["low-4","low-3","low-2","low-1"])
plt.savefig(f'18v31-RNAseq-alluvial.pdf', format='pdf')

with open(f'WTtolowup_gene.csv', 'w') as upfile:
	for item in up:
		upfile.write("%s\n" % item)
with open(f'WTtolowdown_gene.csv', 'w') as downfile:
	for item in down:
		downfile.write("%s\n" % item)

