import numpy as np
import matplotlib.pyplot as plt
import cooler
import cooltools
import bioframe
import cooltools.lib.plotting
from matplotlib.ticker import EngFormatter
from matplotlib.colors import LogNorm


def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)

clr = cooler.Cooler('K562-WT_1000.mcool::/resolutions/1000000')

chromstarts = []
for i in clr.chromnames:
    chromstarts.append(clr.extent(i)[0])

f, axs = plt.subplots(
    figsize=(14,4),
    ncols=3)
bp_formatter = EngFormatter('b')

norm = LogNorm(vmax=100000)
ax = axs[0]
im = ax.matshow(clr.matrix(balance=False)[:], 
    norm=norm,
    cmap = 'fall');
plt.colorbar(im, ax=ax ,fraction=0.046, pad=0.04, label='raw counts');
ax.set_xticks(chromstarts)
ax.set_xticklabels(clr.chromnames)
ax.set_yticks(chromstarts)
ax.set_yticklabels(clr.chromnames)
ax.xaxis.tick_bottom()
ax.set_title('All data (1 Mb resolution)')

norm = LogNorm(vmax=100000)
ax = axs[1]
im = ax.matshow(
    clr.matrix(balance=False).fetch('chr17'),
    norm=norm,
    extent=(0,clr.chromsizes['chr17'], clr.chromsizes['chr17'], 0),
    cmap = 'fall'
);
plt.colorbar(im, ax=ax ,fraction=0.046, pad=0.04, label='raw counts');
ax.set_title('chr17 (1 Mb resolution)', y=1.08)
ax.set_ylabel('position, Mb')
format_ticks(ax)

clr = cooler.Cooler('K562-WT_1000.mcool::/resolutions/100000')
norm = LogNorm(vmax=1000)
ax = axs[2]
start, end = 30_000_000, 60_000_000
region = ('chr17', start, end)
im = ax.matshow(
    clr.matrix(balance=False).fetch(region),
    norm=norm,
    extent=(start, end, end, start),
    cmap = 'fall'
);
ax.set_title('chr17: 30 Mb - 60 Mb (100 kb resolution)', y=1.08)
plt.colorbar(im, ax=ax ,fraction=0.046, pad=0.04, label='raw counts');
format_ticks(ax)
plt.tight_layout()

plt.savefig('K562-WT_100kb_rawmaps.pdf', format='pdf')
