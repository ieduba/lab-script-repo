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
opts.add_option('--r1', help = '<resolution 1> resolution in bp of mcool to use for plotting full contact map (100kb or larger)')
opts.add_option('--r2', help = '<resolution 2> resolution in bp of mcool to use for plotting zoomed contact map')
options, arguments = opts.parse_args()

chrom = 'chr17'
start = 30_000_000
end = 32_000_000
region = (chrom, start, end)

mc_file = options.m
res1 = options.r1
kbres1 = int(int(res1)/1000)
res2 = options.r2
kbres2 = int(int(res2)/1000)

bp_formatter = EngFormatter('b')
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)

clr1 = cooler.Cooler(f'{mc_file}::/resolutions/{res1}')
clr2 = cooler.Cooler(f'{mc_file}::/resolutions/{res2}')

chromstarts = []
for i in clr1.chromnames:
    chromstarts.append(clr1.extent(i)[0])
plt_width=4
f, axs = plt.subplots(
    figsize=( plt_width+plt_width+2, plt_width+plt_width+1),
    ncols=4,
    nrows=3,
    gridspec_kw={'height_ratios':[4,4,1],"wspace":0.01,'width_ratios':[1,.05,1,.05]},
    constrained_layout=True
)

norm = LogNorm(vmax=0.001)
norm_raw = LogNorm(vmin=1, vmax=100)

ax = axs[0,0]
im = ax.matshow(
    clr1.matrix(balance=False)[:],
    norm=norm_raw,
    cmap='fall',
    aspect='auto'
);
ax.xaxis.set_visible(False)
ax.set_title(f'full matrix, {kbres1}kb resolution')
ax.set_ylabel('raw', fontsize=16)

cax = axs[0,1]
plt.colorbar(im, cax=cax, label='raw counts')

ax = axs[1,0]
im = ax.matshow(
    clr1.matrix(balance=True)[:],
    norm=norm,
    cmap='fall',
);
ax.xaxis.set_visible(False)
ax.set_ylabel('balanced', fontsize=16)

cax = axs[1,1]
plt.colorbar(im, cax=cax, label='corrected freqs')

ax1 = axs[2,0]
weights = clr1.bins()[:]['weight'].values
ax1.plot(weights)
ax1.set_xlim([0, len(clr1.bins()[:])])
ax1.set_xlabel('position, bins')

ax1 = axs[2,1]
ax1.set_visible(False)

ax = axs[0,2]
im = ax.matshow(
        clr2.matrix(balance=False).fetch(region),
    norm=LogNorm(vmin=1,vmax=100),
    cmap='fall'
);
ax.set_title(f'{chrom}: {int(int(start)/1000)} kb - {int(int(end)/1000)} kb, {kbres2} kb resolution')
ax.xaxis.set_visible(False)

cax = axs[0,3]
plt.colorbar(im, cax=cax, label='raw counts');

ax = axs[1,2]
im = ax.matshow(
    clr2.matrix(balance=True).fetch(region),
    norm=norm,
    cmap='fall',
    extent=(start, end, end, start)
);
ax.xaxis.set_visible(False)

cax = axs[1,3]
plt.colorbar(im, cax=cax, label='corrected frequencies');

ax1 = axs[2,2]
weights = clr2.bins().fetch(region)['weight'].values
ax1.plot(
    np.linspace(start, end, len(weights)),
    weights
)
format_ticks(ax1, y=False, rotate=False)
ax1.set_xlim(start, end);
ax1.set_xlabel('chr17 position, bp')

ax1 = axs[2,3]
ax1.set_visible(False)

plt.savefig(f'{mc_file}_rawmaps.pdf', format='pdf')
