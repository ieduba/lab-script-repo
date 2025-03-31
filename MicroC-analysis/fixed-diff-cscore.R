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
  filename <- paste0("chr", chrom, "-K36me3.bed")
  df <- data.frame(read.table(filename, col.names = c("chr", "start", "end", "row", "K36me3score", "scrcscore", "lowcscore")))

  corrdf <- correctsign(df)
  alldf <- rbind(alldf, corrdf)
  
  scrplots[[chrom]] <- ggplot(data = df, aes(x=scrcscore, y=K36me3score)) + geom_point(size=0.05) + ggtitle(paste0("chr",chrom))
  lowplots[[chrom]] <- ggplot(data = df, aes(x=lowcscore, y=K36me3score)) + geom_point(size=0.05) + ggtitle(paste0("chr",chrom))
}

write.table(alldf, 'all-cis-5kb_corrected_cscore.bedgraph', sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

pdf("KIM-100kb-scr-cscore-chip-corr.pdf", width = 10, height = 10)
grid.arrange(grobs=scrplots,ncol=4)
dev.off()

pdf("KIM-100kb-low-cscore-chip-corr.pdf", width = 10, height = 10)
grid.arrange(grobs=lowplots,ncol=4)
dev.off()

alldf$shift <- "None"
alldf$shift[is.na(alldf$diff)] <- "undefined"
alldf$shift[alldf$diff > 0.25 & alldf$scrcscore < 0 & alldf$lowcscore < 0] <- "A-shift"
alldf$shift[alldf$diff > 0.25 & alldf$scrcscore < 0 & alldf$lowcscore > 0] <- "B to A"
alldf$shift[alldf$diff > 0.25 & alldf$scrcscore > 0 & alldf$lowcscore > 0] <- "A-shift"
alldf$shift[alldf$diff < -0.25 & alldf$scrcscore > 0 & alldf$lowcscore > 0] <- "B-shift"
alldf$shift[alldf$diff < -0.25 & alldf$scrcscore > 0 & alldf$lowcscore < 0] <- "A to B"
alldf$shift[alldf$diff < -0.25 & alldf$scrcscore < 0 & alldf$lowcscore < 0] <- "B-shift"
alldf$shift[alldf$diff > 0.4 & alldf$lowcscore > 0.7] <- "big-A-shift"
write.table(alldf, file="KH3-100kb_diffcscore.csv", quote = FALSE, row.names = FALSE, sep="\t")

shiftcounts <- alldf %>% count(shift)
write.table(shiftcounts, file="KH3-100kb_scatter-shiftcounts.txt", quote=FALSE, row.names = FALSE, sep = "\t")

plotcolors <- c("black", "maroon", "pink", "darkblue", "lightblue", "transparent", "purple")
names(plotcolors) <- c("None", "A to B", "B-shift", "B to A", "A-shift", "undefined", "big-A-shift")

pdf("KIM-100kb-cscore-scatter-weirddotthresh.pdf", width = 7, height = 5)
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

df_reshape <- data.frame("sample"=c(rep("scr",length(alldf$scrcscore)),rep("low",length(alldf$lowcscore))), "cscore"=c(alldf$scrcscore,alldf$lowcscore))

pdf("KIM-100kb-all-cscores.pdf", width = 5, height = 3)
print(ggplot(data = df_reshape, aes(x=cscore, color=sample)) +
  geom_density() +
  xlab("cscore") +
  ylab("frequency") +
  theme_classic(base_size = 14))
dev.off()
