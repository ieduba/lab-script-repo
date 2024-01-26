library(chromVAR)
library(rtracklayer)
library(SummarizedExperiment)
library(BiocParallel)
library(DESeq2)
register(SerialParam())
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(apeglm)

rawrundir <- "/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23-DE/DEseq"

setwd(paste0(rawrundir))
libs = list.files(rawrundir)
bamFilesVect = libs[grep(".bam", libs)]
bedfile = import("mergedpeaks.bed", format = "BED")

##creates matrix of fragment counts that fall in peaks
fragmentCounts = getCounts(bamFilesVect, bedfile, paired=TRUE, by_rg=FALSE, format="bam")
fragmentCounts = addGCBias(fragmentCounts, genome = BSgenome.Hsapiens.UCSC.hg38)

rownames(fragmentCounts) <- 1:length(bedfile)
colData(fragmentCounts)$Replicates = rep(c("rep1","rep2"), 2)
colData(fragmentCounts)$Treatment = c(rep("H1low",2), rep("scr",2))
colData(fragmentCounts)$sample_id = c("low1","low2","scr1","scr2")

rawCounts = assays(fragmentCounts)[[1]]

counts_filtered = filterSamples(fragmentCounts, min_depth=1500, min_in_peaks=0.15, shiny =FALSE)
counts_filtered = filterPeaks(counts_filtered)

##calls DEseq to find which peaks are differentially expressed between conditions
dds = DESeqDataSetFromMatrix(countData = as.matrix(rawCounts), colData = colData(fragmentCounts), 
                             design= ~Treatment)
dds$Treatment <- relevel(dds$Treatment, ref="scr")

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c('Treatment', 'H1low', 'scr'))

save(dds, file ="DEseqDataSet.RData")
save(res, file ="DEseqResults.RData")

pdf(file = "PlotDisp.pdf")
plotDispEsts(dds)
dev.off()

scrvslowLFC <- lfcShrink(dds, coef="Treatment_H1low_vs_scr", type="apeglm")
#scrvslowLFC <- lfcShrink(dds, coef="Treatment_H1low_vs_scr", type="normal")
scrvslowup <- subset(subset(scrvslowLFC, padj < 0.05), log2FoldChange > 1)
scrvslowdown <- subset(subset(scrvslowLFC, padj < 0.05), log2FoldChange < -1)
write.csv(scrvslowup, "scrvslow-up.csv")
write.csv(scrvslowdown, "scrvslow-down.csv")

pdf(file = "MA.pdf")
plotMA(scrvslowLFC,alpha=0.05,colSig='red')
dev.off()

beddf <- data.frame(peak=1:length(bedfile),chr=seqnames(bedfile), start=start(ranges(bedfile)), end=end(ranges(bedfile)))
updf <- data.frame(peak=rownames(scrvslowup))
upbed <- merge(beddf, updf, by='peak', all=FALSE)
downdf <- data.frame(peak=rownames(scrvslowdown))
downbed <- merge(beddf, downdf, by='peak', all=FALSE)
write.table(upbed[,2:4], "scrvslow-up.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep='\t')
write.table(downbed[,2:4], "scrvslow-down.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep='\t')


widths <- data.frame(width(ranges(bedfile)))
l2fc <- data.frame(scrvslowLFC$log2FoldChange)
widthvsl2fc <- merge(widths, l2fc, by='row.names', all=FALSE)

pdf(file = "l2fc-vs-peaksize.pdf")
plot(widthvsl2fc$width.ranges.bedfile.., widthvsl2fc$scrvslowLFC.log2FoldChange, pch = '.', xlab='peak width (bp)', ylab='log2fold change')
dev.off()
