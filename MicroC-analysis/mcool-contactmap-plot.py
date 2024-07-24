import numpy as np
import matplotlib.pyplot as plt
import cooler
import cooltools
import bioframe
import cooltools.lib.plotting
from matplotlib.ticker import EngFormatter
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from optparse import OptionParser

opts = OptionParser()
opts.add_option('-m', help = '<.mcool> path to mcool file')
opts.add_option('-r', help = '<resolution> resolution in bp of mcool to use for plotting contact map')
opts.add_option('-c', default = 'all', help = '<chr#> chromosome number for contact map region. default = full genome')
opts.add_option('-b', default = 'True', help = '<balanced?> plot ICE balanced contact map? default = True')
opts.add_option('-s', default = '0', help = '<start> optional. start position of contact map region, requires chromosome # passed to -c. default = 0')
opts.add_option('-e', default ='end', help = '<end> optional. end position of contact map region, requires chromosome # passed to -c. default = end of chromosome.')
options, arguments = opts.parse_args()

mc_file = options.m
res = options.r
kbres = int(int(res)/1000)

print(f'loading {mc_file} at resolution {kbres}kb')

bp_formatter = EngFormatter('b')
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)

clr = cooler.Cooler(f'{mc_file}::/resolutions/{res}')

bal = eval(options.b)
chrom = options.c
start = options.s
end = options.e

if bal == True:
	norm = LogNorm(vmax=0.1)
	balstr = 'balanced'
else:
	norm = LogNorm(vmax = 100)
	balstr = 'raw'

print(f'plotting {chrom} from {start} to {end}, {balstr}')

chromstarts = []
for i in clr.chromnames:
    chromstarts.append(clr.extent(i)[0])

f,ax = plt.subplots(figsize=(7,5))

ax = ax
ax.set_title(f'{chrom}: {start}-{end}, {kbres}kb resolution')
ax.set_ylabel(balstr, fontsize=16)

if chrom == 'all':
	im = ax.matshow(clr.matrix(balance=bal)[:], norm=norm, cmap='fall', aspect='equal')
	ax.set(yticks=[chromstarts[i] for i in range(0,23,2)],
		yticklabels=[clr.chromnames[i] for i in range(0,23,2)])
	ax.set(xticks=[chromstarts[i] for i in range(0,23,2)],
		xticklabels=[clr.chromnames[i] for i in range(0,23,2)])
	ax.xaxis.tick_bottom()
	ax.tick_params(axis='x',rotation=45)
else:
	if end == 'end':
		end = clr.chromsizes[chrom]
	region = (chrom, start, end)
	im = ax.matshow(clr.matrix(balance=bal).fetch(region), norm=norm, cmap='fall', 
		aspect='equal', extent=(int(start),int(end),int(end),int(start)))
	format_ticks(ax)

plt.colorbar(im, ax=ax, fraction = 0.046, pad = 0.04, label=f'{balstr} counts')

plt.savefig(f'{mc_file}_{chrom}-{start}-{end}_{balstr}map.pdf', format='pdf')
