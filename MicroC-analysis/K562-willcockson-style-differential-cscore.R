library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

datadir <- "/lustre/fs4/risc_lab/scratch/iduba/linker-histone/HiC/KH3/deepseq/combinedreps/chroms"
#datadir <- "/ru-auth/local/home/iduba/linker-histone/Micro-C/in-situ/old-new-combo/combo-bio-reps/compartments/chroms"
setwd(datadir)

correctsign <- function(df){
  scrrho <- cor.test(df$K36me3score, df$scrcscore, method = "spearman")
  lowrho <- cor.test(df$K36me3score, df$lowcscore, method = "spearman")
  if (scrrho$estimate < 0){df$scrcscore <- -df$scrcscore}
  if (lowrho$estimate < 0){df$lowcscore <- -df$lowcscore}
  df$diff <- df$lowcscore - df$scrcscore
  if (scrrho$p.value > 0.001 | lowrho$p.value > 0.01){df$diff <- NA}
  df
}

alldf <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("chr", "start", "end", "row", "K36me3score", "scrcscore", "lowcscore", "diff"))
scrplots <- list()
lowplots <- list()
for (chrom in c(1:22)){
  filename <- paste0("chr", chrom, "-scr-K36me3.bed")
  df <- data.frame(read.table(filename, col.names = c("chr", "start", "end", "row", "K36me3score", "scrcscore", "lowcscore")))

  corrdf <- correctsign(df)
  alldf <- rbind(alldf, corrdf)
  
  scrplots[[chrom]] <- ggplot(data = df, aes(x=scrcscore, y=K36me3score)) + geom_point(size=0.05) + ggtitle(paste0("chr",chrom))
  lowplots[[chrom]] <- ggplot(data = df, aes(x=lowcscore, y=K36me3score)) + geom_point(size=0.05) + ggtitle(paste0("chr",chrom))
}

write.table(alldf, 'all-cis-5kb_corrected_cscore.bedgraph', sep="\t", row.names=FALSE, quote=FALSE)

pdf("KIM-100kb-scr1-cscore-chip-corr.pdf", width = 10, height = 10)
grid.arrange(grobs=scrplots,ncol=4)
dev.off()

pdf("KIM-100kb-scr2-cscore-chip-corr.pdf", width = 10, height = 10)
grid.arrange(grobs=lowplots,ncol=4)
dev.off()

alldf$shift <- "None"
alldf$shift[is.na(alldf$diff)] <- "undefined"
alldf$shift[alldf$diff > 0.25 & alldf$scrcscore < 0 & alldf$lowcscore < 0] <- "decompaction"
alldf$shift[alldf$diff > 0.25 & alldf$scrcscore < 0 & alldf$lowcscore > 0] <- "B to A"
alldf$shift[alldf$diff > 0.25 & alldf$scrcscore > 0 & alldf$lowcscore > 0] <- "decompaction"
alldf$shift[alldf$diff < -0.25 & alldf$scrcscore > 0 & alldf$lowcscore > 0] <- "compaction"
alldf$shift[alldf$diff < -0.25 & alldf$scrcscore > 0 & alldf$lowcscore < 0] <- "A to B"
alldf$shift[alldf$diff < -0.25 & alldf$scrcscore < 0 & alldf$lowcscore < 0] <- "compaction"
write.table(alldf, file="KIM-100kb_diffcscore.csv", quote = FALSE, row.names = FALSE, sep="\t")

shiftcounts <- alldf %>% count(shift)
write.table(shiftcounts, file="scr-KIM-100kb_scatter-shiftcounts.txt", quote=FALSE, row.names = FALSE, sep = "\t")

plotcolors <- c("black", "maroon", "pink", "darkblue", "lightblue", "transparent")
names(plotcolors) <- c("None", "A to B", "compaction", "B to A", "decompaction", "undefined")

pdf("scr-KIM-100kb-cscore-scatter.pdf", width = 7, height = 5)
print(ggplot(data = alldf, aes(x=scrcscore, y=lowcscore, color = shift)) +
        geom_point(size=0.01) +
        scale_color_manual(values = plotcolors) +
        geom_vline(xintercept=0, size=1) + 
        geom_hline(yintercept=0, size=1)+
        geom_abline(intercept = 0, slope = 1, colour="black", linetype="longdash")+
        geom_abline(intercept = .25, slope = 1, colour="#787878", linetype="dashed")+
        geom_abline(intercept = -.25, slope = 1, colour="#787878", linetype="dashed")+
        xlab("scr cscore") +
        ylab("H1-low cscore") +
        theme_classic(base_size = 14) +
        guides(colour = guide_legend(override.aes = list(size=2))))
dev.off()

df_reshape <- data.frame("sample"=c(rep("scr1",length(alldf$scrcscore)),rep("scr2",length(alldf$lowcscore))), "cscore"=c(alldf$scrcscore,alldf$lowcscore))

pdf("scr-KIM-100kb-all-cscores.pdf", width = 5, height = 3)
print(ggplot(data = df_reshape, aes(x=cscore, color=sample)) +
  geom_density() +
  xlab("cscore") +
  ylab("frequency") +
  theme_classic(base_size = 14))
dev.off()


setwd("/ru-auth/local/home/iduba/linker-histone/HiC/KH3/deepseq/combinedreps/chroms/annotations")

genesdf <- data.frame(read.table("all-cis-5kb_corrected_cscore-upgenecounts.bedgraph", col.names = c("chr", "start", "end", "row", "K36me3score", "scrcscore", "lowcscore", "diff", "upgenes")))
genesdf <- genesdf[!is.na(genesdf$diff),]
genesdf <- genesdf[order(genesdf$upgenes),]
pdf("KH3-100kb-cscore-upgenes.pdf", width = 7, height = 5)
print(ggplot(data = genesdf, aes(x=scrcscore, y=lowcscore, color = upgenes)) +
        geom_point(size=0.1) +
        scale_color_gradient(low="navy", high="magenta", trans="log") +
        geom_vline(xintercept=0, size=0.6) + 
        geom_hline(yintercept=0, size=0.6)+
        xlab("scr cscore") +
        ylab("H1-low cscore") +
        ggtitle("number upregulated genes per bin (log scale)") +
        theme_classic(base_size = 14))
dev.off()

plotscore <- function(df, mark){
  df <- df[order(df$scrcnt),]
  pdf(paste0('KH3-100kb-cscore-scr',mark,'cnt.pdf'), width = 7, height = 5)
  print(ggplot(data = df, aes(x=scrcscore, y=lowcscore, color = scrcnt)) +
          geom_point(size=0.1) +
          scale_color_gradient(low="navy", high="magenta") +
          geom_vline(xintercept=0, size=0.6) + 
          geom_hline(yintercept=0, size=0.6)+
          xlab("scr cscore") +
          ylab("H1-low cscore") +
          ggtitle(paste0('average scr ', mark, ' Cut&Tag score per bin')) +
          theme_classic(base_size = 14))
  dev.off()
  df <- df[order(df$lowcnt),]
  pdf(paste0('KH3-100kb-cscore-low',mark,'cnt.pdf'), width = 7, height = 5)
  print(ggplot(data = df, aes(x=scrcscore, y=lowcscore, color = lowcnt)) +
          geom_point(size=0.1) +
          scale_color_gradient(low="navy", high="magenta") +
          geom_vline(xintercept=0, size=0.6) + 
          geom_hline(yintercept=0, size=0.6)+
          xlab("scr cscore") +
          ylab("H1-low cscore") +
          ggtitle(paste0('average low ', mark, ' Cut&Tag score per bin')) +
          theme_classic(base_size = 14))
  dev.off()
}

k9df <- data.frame(read.table("all-corrected-cscore-H3k9me3-cutntag.bedgraph", header=TRUE))
plotscore(k9df[!is.na(k9df$diff),], 'K9me3')
k27df <- data.frame(read.table("all-corrected-cscore-H3k27me3-cutntag.bedgraph", header=TRUE))
plotscore(k27df[!is.na(k27df$diff),], 'K27me3')
k36df <- data.frame(read.table("all-corrected-cscore-H3k36me2-cutntag.bedgraph", header=TRUE))
plotscore(k36df[!is.na(k36df$diff),], 'K36me2')
h1p2df <- data.frame(read.table("all-corrected-cscore-H1p2-cutntag.bedgraph", header=TRUE))
plotscore(h1p2df[!is.na(h1p2df$diff),], 'H1p2')
panh1df <- data.frame(read.table("all-corrected-cscore-panH1-cutntag.bedgraph", header=TRUE))
plotscore(panh1df[!is.na(panh1df$diff),], 'panH1')

l2fcdf <- data.frame(read.table("all-corrected-cscore-avgl2fc.bedgraph", col.names = c("chr", "start", "end", "row", "K36me3score", "scrcscore", "lowcscore", "diff", "l2fc")))
basemeandf <- data.frame(read.table("all-corrected-cscore-avgbasemean.bedgraph", col.names = c("chr", "start", "end", "row", "K36me3score", "scrcscore", "lowcscore", "diff", "basemean")))

#df <- l2fcdf[order(l2fcdf$l2fc),]
df <- l2fcdf
pdf('KH3-cscore-100kb-RNAl2fc.pdf', width = 7, height = 5)
print(ggplot(data = df, aes(x=scrcscore, y=lowcscore, color = l2fc)) +
        geom_point(size=0.2) +
        scale_color_gradient(low="blue", high="red") +
        geom_vline(xintercept=0, size=0.6) + 
        geom_hline(yintercept=0, size=0.6)+
        xlab("scr cscore") +
        ylab("H1-low cscore") +
        ggtitle("average expression log2 fold-change per bin") +
        theme_classic(base_size = 14))
dev.off()
#df <- basemeandf[order(basemeandf$basemean),]
df <- basemeandf
pdf(paste0('KH3-cscore-100kb-RNAbasemean.pdf'), width = 7, height = 5)
print(ggplot(data = df, aes(x=scrcscore, y=lowcscore, color = basemean)) +
        geom_point(size=0.2) +
        scale_color_gradient(low="blue", high="red", trans="log") +
        geom_vline(xintercept=0, size=0.6) + 
        geom_hline(yintercept=0, size=0.6)+
        xlab("scr cscore") +
        ylab("H1-low cscore") +
        ggtitle("average expression base mean per bin (log scale)") +
        theme_classic(base_size = 14))
dev.off()
