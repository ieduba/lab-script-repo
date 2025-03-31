library(chromVAR)
library(rtracklayer)
library(SummarizedExperiment)
library(BiocParallel)
library(DESeq2)
register(SerialParam())
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(apeglm)

rawrundir <- "/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/DEseq"

setwd(paste0(rawrundir))
#libs = list.files(rawrundir)
#bamFilesVect = libs[grep("-no5kbTSS.bam", libs)]
bamFilesVect=c("scr-KA2.bam", "scr-KA3.bam", "dH1-KA2.bam", "dH1-KA3.bam")
bedfile = import("reps-mergedpeaks.bed", format = "BED")

##creates matrix of fragment counts that fall in peaks
fragmentCounts = getCounts(bamFilesVect, bedfile, paired=TRUE, by_rg=FALSE, format="bam")
fragmentCounts = addGCBias(fragmentCounts, genome = BSgenome.Hsapiens.UCSC.hg38)

rownames(fragmentCounts) <- 1:length(bedfile)
colData(fragmentCounts)$Replicates = rep(c("rep1","rep2"), 2)
colData(fragmentCounts)$Treatment = c(rep("scr",2), rep("H1low",2))
colData(fragmentCounts)$sample_id = c("scr1","scr2","low1","low2")

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

save(dds, file ="rep-peaks-DEseqDataSet.RData")
save(res, file ="rep-peaks-DEseqResults.RData")

pdf(file = "rep-peaks-PlotDisp.pdf")
plotDispEsts(dds)
dev.off()

scrvslowLFC <- lfcShrink(dds, coef="Treatment_H1low_vs_scr", type="apeglm")
#scrvslowLFC <- lfcShrink(dds, coef="Treatment_H1low_vs_scr", type="normal")
scrvslowup <- subset(subset(scrvslowLFC, padj < 0.05), log2FoldChange > 1)
scrvslowdown <- subset(subset(scrvslowLFC, padj < 0.05), log2FoldChange < -1)
scrvslownochange <- subset(subset(scrvslowLFC, padj > 0.05), log2FoldChange > -1 & log2FoldChange < 1)
write.csv(scrvslowup, "rep-peaks-scrvslow-up.csv")
write.csv(scrvslowdown, "rep-peaks-scrvslow-down.csv")
write.csv(scrvslowLFC, "rep-peak-allLFC.csv")
write.csv(scrvslownochange, "rep-peak-scrvslow-nochange.csv")

pdf(file = "rep-peaks-MA.pdf")
plotMA(scrvslowLFC,alpha=0.05,colSig='red')
dev.off()

beddf <- data.frame(peak=1:length(bedfile),chr=seqnames(bedfile), start=start(ranges(bedfile)), end=end(ranges(bedfile)))
updf <- data.frame(peak=rownames(scrvslowup), baseMean=scrvslowup$baseMean, log2FoldChange=scrvslowup$log2FoldChange, padj=scrvslowup$padj)
upbed <- merge(beddf, updf, by='peak', all=FALSE)
downdf <- data.frame(peak=rownames(scrvslowdown), baseMean=scrvslowdown$baseMean, log2FoldChange=scrvslowdown$log2FoldChange, padj=scrvslowdown$padj)
downbed <- merge(beddf, downdf, by='peak', all=FALSE)
nochangedf <- data.frame(peak=rownames(scrvslownochange), baseMean=scrvslownochange$baseMean, log2FoldChange=scrvslownochange$log2FoldChange, padj=scrvslownochange$padj)
nochangebed <- merge(beddf, nochangedf, by='peak', all=FALSE)
alldf <- data.frame(peak=rownames(scrvslowLFC), baseMean=scrvslowLFC$baseMean, log2FoldChange=scrvslowLFC$log2FoldChange, padj=scrvslowLFC$padj)
allbed <- merge(beddf, alldf, by='peak', all=FALSE)
write.table(upbed[,2:7], "rep-peaks-scrvslow-up.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep='\t')
write.table(downbed[,2:7], "rep-peaks-scrvslow-down.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep='\t')
write.table(nochangebed[,2:7], "rep-peaks-scrvslow-nochange.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep='\t')
write.table(allbed[,2:7], "rep-peaks-scrvslow-LFC.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep='\t')

widths <- data.frame(width(ranges(bedfile)))
l2fc <- data.frame(scrvslowLFC$log2FoldChange)
widthvsl2fc <- merge(widths, l2fc, by='row.names', all=FALSE)

pdf(file = "rep-peaks-l2fc-vs-peaksize.pdf")
plot(widthvsl2fc$width.ranges.bedfile.., widthvsl2fc$scrvslowLFC.log2FoldChange, pch = '.', xlab='peak width (bp)', ylab='log2fold change')
dev.off()
