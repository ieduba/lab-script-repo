#!/usr/bin/env Rscript

library(chromVAR)
library(rtracklayer)
library(SummarizedExperiment)
library(BiocParallel)
library(DESeq2)
register(SerialParam())
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)

rawrundir <- "/Users/irene/Documents/Rockefeller/Rotations/Risca Lab/TestRun"

# Turn this into function
setwd(paste0(rawrundir))
#libs = list.files(rawrundir)
#libs = libs[grep("LSATAC", libs)]
libs = basename(list.dirs(rawrundir, recursive = FALSE))

#listPeakFiles = list.files(rawrundir)
bedfile = import("LSATACmergedpeaks.bed", format = "BED")
#bedfile = listPeakFiles[grep("narrowPeak.bed", listPeakFiles)]
#narrowPeak = readNarrowpeaks(bedfile, width=500, non_overlapping=TRUE)	

bamFilesVect = c()

##creates vector that is list of paths of all completely filtered bam files
for (i in 1:length(libs)) {
  #setwd(paste0(rawrundir, "/", libs[i], "/00_source"))
  setwd(paste0(rawrundir, "/", libs[i]))
  bamFile = list.files(getwd(), pattern ="rmdup.bam")
  bamPath= paste0(getwd(), "/", bamFile)
  bamFilesVect = c(bamFilesVect, bamPath)
}

setwd(paste0(rawrundir))
fragmentCounts = getCounts(bamFilesVect, bedfile, paired=TRUE, by_rg=FALSE, format="bam")
fragmentCounts = addGCBias(fragmentCounts, genome = BSgenome.Hsapiens.UCSC.hg19)

colData(fragmentCounts)$Replicates = c(rep("rep1", 3), rep("rep2", 3), rep("rep3", 3))
colData(fragmentCounts)$Conditions = rep(c("Cycling", "Quiescent", "Senescent"), 3)
colData(fragmentCounts)$sample_id = libs

rawCounts = assays(fragmentCounts)[[1]]

counts_filtered = filterSamples(fragmentCounts, min_depth=1500, min_in_peaks=0.15, shiny =FALSE)
counts_filtered = filterPeaks(counts_filtered)

dds = DESeqDataSetFromMatrix(countData = as.matrix(rawCounts), colData = colData(fragmentCounts), 
                             design= ~Conditions, rowRanges= rowRanges(fragmentCounts))
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)

pdf(file = "PlotDisp.pdf")
plotDispEsts(dds)
dev.off()

CyclingvsQuiescent = results(dds, c("Conditions", "Cycling", "Quiescent"), format="GRanges", alpha=0.01)
CyclingvsQuiescent = CyclingvsQuiescent[order(CyclingvsQuiescent$padj),]

UpCyclingvsQuiescent = CyclingvsQuiescent[intersect(which(CyclingvsQuiescent$log2FoldChange > 2),which(CyclingvsQuiescent$padj < 0.01)),]
DownCyclingvsQuiescent = CyclingvsQuiescent[intersect(which(CyclingvsQuiescent$log2FoldChange < -2), which(CyclingvsQuiescent$padj < 0.01)),]
print(length(UpCyclingvsQuiescent))
print(length(DownCyclingvsQuiescent))
averageNumCQ = (length(UpCyclingvsQuiescent) + length(DownCyclingvsQuiescent))/2
NormCyclingvsQuiescent = CyclingvsQuiescent[intersect(which(CyclingvsQuiescent$log2FoldChange > -1 & CyclingvsQuiescent$log2FoldChange < 1 ), which(CyclingvsQuiescent$padj >0.1))] 
NormCyclingvsQuiescent = sample(NormCyclingvsQuiescent, averageNumCS)

pdf(file = "MA_CyclingvsQuiescent.pdf")
CyclingvsQuiescentx = results(dds, c("Conditions", "Cycling", "Quiescent"), lfcThreshold=1,  alpha =0.01)
plotMA(CyclingvsQuiescentx, alpha=0.01, ylim = c(-4,4))
dev.off()

export(NormCyclingvsQuiescent, "NormCyclingvsQuiescent.bed", format = "BED")
export(UpCyclingvsQuiescent, "UpCyclingvsQuiescent.bed", format="BED")
export(DownCyclingvsQuiescent, "DownCyclingvsQuiescent.bed", format="BED")

pdf(file = "MA_CyclingvsSenescent.pdf")
CyclingvsSenescent = results(dds, c("Conditions", "Cycling", "Senescent"), format="GRanges", alpha =0.01)
CyclingvsSenescent = CyclingvsSenescent[order(CyclingvsSenescent$padj)]
sumSigPeaksCS = sum(CyclingvsSenescent$padj < 0.01, na.rm=TRUE)
UpCyclingvsSenescent = CyclingvsSenescent[intersect(which(CyclingvsSenescent$log2FoldChange > 2), which(CyclingvsSenescent$padj <0.01)),]
DownCyclingvsSenescent = CyclingvsSenescent[intersect(which(CyclingvsSenescent$log2FoldChange < -2), which(CyclingvsSenescent$padj < 0.01)),]
print(length(UpCyclingvsSenescent))
print(length(DownCyclingvsSenescent))
averageNumCS = (length(UpCyclingvsSenescent) + length(DownCyclingvsSenescent))/2
NormCyclingvsSenescent = CyclingvsSenescent[intersect(which(CyclingvsSenescent$log2FoldChange > -1 & CyclingvsSenescent$log2FoldChange < 1 ), which(CyclingvsSenescent$padj >0.1))] 
NormCyclingvsSenescent = sample(NormCyclingvsSenescent, averageNumCS)

CyclingvsSenescentx = results(dds, c("Conditions", "Cycling", "Senescent"), lfcThreshold= 1, alpha =0.01)
plotMA(CyclingvsSenescentx, alpha=0.01, ylim = c(-4,4))
dev.off()

export(NormCyclingvsSenescent, "NormCyclingvsSenescent.bed", format = "BED")
export(UpCyclingvsSenescent, "UpCyclingvsSenescent.bed", format="BED")
export(DownCyclingvsSenescent, "DownCyclingvsSenescent.bed", format="BED")

pdf(file = "MA_QuiescentvsSenescent.pdf")
QuiescentvsSenescent = results(dds, c("Conditions", "Quiescent", "Senescent"), format="GRanges", alpha =0.01)
QuiescentvsSenescent = QuiescentvsSenescent[order(QuiescentvsSenescent$padj),]
sumSigPeaksQS = sum(QuiescentvsSenescent$padj < 0.01, na.rm=TRUE)
UpQuiescentvsSenescent = QuiescentvsSenescent[intersect(which(QuiescentvsSenescent$log2FoldChange > 2), which(QuiescentvsSenescent$padj < 0.01)),]
DownQuiescentvsSenescent = QuiescentvsSenescent[intersect(which(QuiescentvsSenescent$log2FoldChange < -2), which(QuiescentvsSenescent$padj < 0.01)),]

print(length(UpQuiescentvsSenescent))
print(length(DownQuiescentvsSenescent))
averageNumQS = (length(UpQuiescentvsSenescent) + length(DownQuiescentvsSenescent))/2
NormQuiescentvsSenescent = QuiescentvsSenescent[intersect(which(QuiescentvsSenescent$log2FoldChange > -1 & QuiescentvsSenescent$log2FoldChange < 1 ), which(QuiescentvsSenescent$padj >0.1))] 
NormQuiescentvsSenescent = sample(NormQuiescentvsSenescent, averageNumQS)


QuiescentvsSenescentx = results(dds, c("Conditions", "Quiescent", "Senescent"), lfcThreshold = 1,  alpha =0.01)
plotMA(QuiescentvsSenescentx, alpha = 0.01, ylim= c(-4,4))
dev.off()

export(NormQuiescentvsSenescent, "NormQuiescentvsSenescent.bed", format = "BED")
export(UpQuiescentvsSenescent, "UpQuiescentvsSenescent.bed", format="BED")
export(DownQuiescentvsSenescent, "DownQuiescentvsSenescent.bed", format="BED")
