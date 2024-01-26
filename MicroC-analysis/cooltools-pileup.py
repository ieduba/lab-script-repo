import numpy as np
import matplotlib.pyplot as plt
import cooler
import cooltools
import bioframe
import cooltools.lib.plotting
from optparse import OptionParser

opts = OptionParser()
opts.add_option('-m', help = '<.mcool> path to mcool file')
opts.add_option('-r', help = '<resoltion> resolution of mcool to use for p
lotting contact map')
options, arguments = opts.parse_args()


sample = 'K562_WT_50_test'
resolution = 50
flank = 1000
tss_file = '/ru-auth/local/home/iduba/linker-histone/Micro-C/annotations/hg38_TSS_export_bedrearrange_gene.bed'
uptss_file = '/ru-auth/local/home/iduba/linker-histone/RNA-seq/human/DESeq2_results/WTvslowup-fixed-symbol-TSS.bed'
downtss_file = '/ru-auth/local/home/iduba/linker-histone/RNA-seq/human/DESeq2_results/WTvslowdown-fixed-symbol-TSS.bed'
ctcf_file = '/rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/annotations/chip/K562/ENCFF396BZQ.bed' 

def pileupmat(infile):
	allsites = bioframe.read_table(infile, schema='bed')
	filtsites = bioframe.cluster(allsites, min_dist=resolution).drop_duplicates('cluster').reset_index(drop=True)
	stack = cooltools.pileup(clr, filtsites, view_df = hg38_arms, flank = flank)
	mask = np.array(filtsites.strand == '-', dtype = bool)
	stack[:, :, mask] = stack[::-1, ::-1, mask]
	mtx = np.nanmean(stack, axis=2)
	return mtx

def plotpileup(mtx, anno):
	plt.imshow(np.log10(mtx),cmap='autumn',interpolation='none',vmin=-5,vmax=-2.5)
	plt.colorbar(label = 'log10 mean ICed Hi-C')
	ticks_pixels = np.linspace(0, flank*2//resolution,5)
	ticks_kbp = ((ticks_pixels-ticks_pixels[-1]/2)*resolution//1000).astype(int)
	plt.xticks(ticks_pixels, ticks_kbp)
	plt.yticks(ticks_pixels, ticks_kbp)
	plt.xlabel('relative position, kbp')
	plt.ylabel('relative position, kbp')
	plt.savefig(f'{sample}-{anno}.pdf', format='pdf', bbox_inches='tight')
	plt.close()


clr = cooler.Cooler(f'{sample}.mcool::/resolutions/{resolution}')
hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_cens = bioframe.fetch_centromeres('hg38')
hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)

tssmat = pileupmat(tss_file)
uptssmat = pileupmat(uptss_file)
downtssmat = pileupmat(downtss_file)
ctcfmat = pileupmat(ctcf_file)

plotpileup(tssmat, 'TSS-1kb')
plotpileup(uptssmat, 'upgenes-TSS-1kb')
plotpileup(downtssmat, 'downgenes-TSS-1kb')
plotpileup(ctcfmat, 'CTCF-1kb')
