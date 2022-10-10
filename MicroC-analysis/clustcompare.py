import alluvial
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np

celltype = 'K562_WT'
nclust = 4
tr1 = open(f'TR1-{celltype}-allho-{nclust}-binclusters.csv', 'r')
tr2 = open(f'TR2-{celltype}-allho-{nclust}-binclusters.csv', 'r')
line1 = tr1.readline()
line2 = tr2.readline()
clusts = []
while line1:
	line1col = line1.split(',')
	line2col = line2.split(',')
	clusts.append([f'TR1-{line1col[-1][0]}', f'TR2-{line2col[-1][0]}'])
	line1 = tr1.readline()
	line2 = tr2.readline()
ax = alluvial.plot(clusts, disp_width = True, wdisp_sep = ' '*2, figsize=(7,5))
plt.savefig(f'{celltype}-{nclust}-alluvial.pdf', format='pdf')
