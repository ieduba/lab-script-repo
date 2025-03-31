library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

datadir <- "/ru-auth/local/home/iduba/labWTatch/iduba/arnold-cscore/chroms"
setwd(datadir)

correctsign <- function(df){
  WTrho <- cor.test(df$K36me3score, df$WTcscore, method = "spearman")
  KOrho <- cor.test(df$K36me3score, df$KOcscore, method = "spearman")
  if (WTrho$estimate < 0){df$WTcscore <- -df$WTcscore}
  if (KOrho$estimate < 0){df$KOcscore <- -df$KOcscore}
  df$diff <- df$KOcscore - df$WTcscore
  if (WTrho$p.value > 0.001 | KOrho$p.value > 0.01){df$diff <- NA}
  df
}

alldf <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("chr", "start", "end", "row", "K36me3score", "WTcscore", "KOcscore", "diff"))
WTplots <- list()
KOplots <- list()
for (chrom in c(1:22)){
  filename <- paste0("combo-chr", chrom, "-K36me3.bed")
  df <- data.frame(read.table(filename, col.names = c("chr", "start", "end", "row", "K36me3score", "WTcscore", "KOcscore")))

  corrdf <- correctsign(df)
  alldf <- rbind(alldf, corrdf)
  
  WTplots[[chrom]] <- ggplot(data = df, aes(x=WTcscore, y=K36me3score)) + geom_point(size=0.05) + ggtitle(paste0("chr",chrom))
  KOplots[[chrom]] <- ggplot(data = df, aes(x=KOcscore, y=K36me3score)) + geom_point(size=0.05) + ggtitle(paste0("chr",chrom))
}

write.table(alldf, 'all-combined-cis-5kb_corrected_cscore.bedgraph', sep="\t", row.names=FALSE, quote=FALSE)

pdf("100kb-WT1-cscore-chip-corr.pdf", width = 10, height = 10)
grid.arrange(grobs=WTplots,ncol=4)
dev.off()

pdf("100kb-KO-cscore-chip-corr.pdf", width = 10, height = 10)
grid.arrange(grobs=KOplots,ncol=4)
dev.off()

alldf$shift <- "None"
alldf$shift[is.na(alldf$diff)] <- "undefined"
alldf$shift[alldf$diff > 0.25 & alldf$WTcscore < 0 & alldf$KOcscore < 0] <- "decompaction"
alldf$shift[alldf$diff > 0.25 & alldf$WTcscore < 0 & alldf$KOcscore > 0] <- "B to A"
alldf$shift[alldf$diff > 0.25 & alldf$WTcscore > 0 & alldf$KOcscore > 0] <- "decompaction"
alldf$shift[alldf$diff < -0.25 & alldf$WTcscore > 0 & alldf$KOcscore > 0] <- "compaction"
alldf$shift[alldf$diff < -0.25 & alldf$WTcscore > 0 & alldf$KOcscore < 0] <- "A to B"
alldf$shift[alldf$diff < -0.25 & alldf$WTcscore < 0 & alldf$KOcscore < 0] <- "compaction"
write.table(alldf, file="combo-corr_diffcscore.csv", quote = FALSE, row.names = FALSE, sep="\t")

shiftcounts <- alldf %>% count(shift)
write.table(shiftcounts, file="combo-scatter-cscore-shiftcounts.txt", quote=FALSE, row.names = FALSE, sep = "\t")

plotcolors <- c("black", "maroon", "pink", "darkblue", "lightblue", "transparent")
names(plotcolors) <- c("None", "A to B", "compaction", "B to A", "decompaction", "undefined")

pdf("100kb-combo-cscore-scatter.pdf", width = 7, height = 5)
print(ggplot(data = alldf, aes(x=WTcscore, y=KOcscore, color = shift)) +
        geom_point(size=0.01) +
        scale_color_manual(values = plotcolors) +
        geom_vline(xintercept=0, size=1) + 
        geom_hline(yintercept=0, size=1)+
        geom_abline(intercept = 0, slope = 1, colour="black", linetype="longdash")+
        geom_abline(intercept = .25, slope = 1, colour="#787878", linetype="dashed")+
        geom_abline(intercept = -.25, slope = 1, colour="#787878", linetype="dashed")+
        xlab("WT cscore") +
        ylab("H1-KO cscore") +
        theme_classic(base_size = 14) +
        guides(colour = guide_legend(override.aes = list(size=2))))
dev.off()

df_reshape <- data.frame("sample"=c(rep("WT",length(alldf$WTcscore)),rep("KO",length(alldf$KOcscore))), "cscore"=c(alldf$WTcscore,alldf$KOcscore))

pdf("100kb-combo-all-cscores.pdf", width = 5, height = 3)
print(ggplot(data = df_reshape, aes(x=cscore, color=sample)) +
  geom_density() +
  xlab("cscore") +
  ylab("frequency") +
  theme_classic(base_size = 14))
dev.off()

genesdf <- data.frame(read.table("all-cis-5kb_corrected_cscore-upgenecounts.bedgraph", col.names = c("chr", "start", "end", "row", "K36me3score", "WTcscore", "KOcscore", "diff", "upgenes")))
genesdf <- genesdf[order(genesdf$upgenes),]
pdf("KH3-100kb-cscore-upgenes.pdf", width = 7, height = 5)
print(ggplot(data = genesdf, aes(x=WTcscore, y=KOcscore, color = upgenes)) +
        geom_point(size=0.1) +
        scale_color_gradient(KO="blue", high="magenta", trans="log") +
        geom_vline(xintercept=0, size=0.6) + 
        geom_hline(yintercept=0, size=0.6)+
        xlab("WT cscore") +
        ylab("H1-KO cscore") +
        theme_classic(base_size = 14))
dev.off()
