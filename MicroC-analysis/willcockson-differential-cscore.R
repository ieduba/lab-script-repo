library(ggplot2)
library(dplyr)

datadir <- "/lustre/fs4/risc_lab/scratch/iduba/linker-histone/HiC/willcockson21/compartments/combined-reps/chroms"
setwd(datadir)

correctsign <- function(df){
  WTrho <- cor.test(df$K36me2score, df$WTcscore, method = "spearman")
  TKOrho <- cor.test(df$K36me2score, df$TKOcscore, method = "spearman")
  if (WTrho$estimate < 0){df$WTcscore <- -df$WTcscore}
  if (TKOrho$estimate < 0){df$TKOcscore <- -df$TKOcscore}
  df
}

alldf <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("chr", "start", "end", "row", "K36me2score", "WTcscore", "TKOcscore", "diff"))
for (chrom in c(1:19,'X')){
  filename <- paste0("chr", chrom, "-K36me2.bed")
  df <- data.frame(read.table(filename, col.names = c("chr", "start", "end", "row", "K36me2score", "WTcscore", "TKOcscore")))
  corrdf <- correctsign(df)
  corrdf$diff <- corrdf$TKOcscore - corrdf$WTcscore
  alldf <- rbind(alldf, corrdf)
}

alldf$shift <- "None"
alldf$shift[alldf$diff > 0.25 & alldf$WTcscore < 0 & alldf$TKOcscore < 0] <- "decompaction"
alldf$shift[alldf$diff > 0.25 & alldf$WTcscore < 0 & alldf$TKOcscore > 0] <- "B to A"
alldf$shift[alldf$diff > 0.25 & alldf$WTcscore > 0 & alldf$TKOcscore > 0] <- "decompaction"
alldf$shift[alldf$diff < -0.25 & alldf$WTcscore > 0 & alldf$TKOcscore > 0] <- "compaction"
alldf$shift[alldf$diff < -0.25 & alldf$WTcscore > 0 & alldf$TKOcscore < 0] <- "A to B"
alldf$shift[alldf$diff < -0.25 & alldf$WTcscore < 0 & alldf$TKOcscore < 0] <- "compaction"
write.table(alldf, file="willcockson100kb_diffcscore.csv", quote = FALSE, row.names = FALSE, sep="\t")

shiftcounts <- alldf %>% count(shift)
write.table(shiftcounts, file="willcockson100kb_scatter-shiftcounts.txt", quote=FALSE, row.names = FALSE, sep = "\t")

plotcolors <- c("black", "maroon", "pink", "darkblue", "lightblue")
names(plotcolors) <- c("None", "A to B", "compaction", "B to A", "decompaction")

figname <- "willcockson100kb-cscore-scatter.pdf"
pdf(figname, width = 7, height = 5)
print(ggplot(data = alldf, aes(x=WTcscore, y=TKOcscore, color = shift)) +
        geom_point(size=0.01) +
        scale_color_manual(values = plotcolors) +
        geom_vline(xintercept=0, size=1) + 
        geom_hline(yintercept=0, size=1)+
        geom_abline(intercept = 0, slope = 1, colour="black", linetype="longdash")+
        geom_abline(intercept = .25, slope = 1, colour="#787878", linetype="dashed")+
        geom_abline(intercept = -.25, slope = 1, colour="#787878", linetype="dashed")+
        xlab("scr cscore") +
        ylab("H1-TKO cscore") +
        theme_classic(base_size = 14) +
        guides(colour = guide_legend(override.aes = list(size=2))))
dev.off()

df_reshape <- data.frame("sample"=c(rep("WT",length(alldf$WTcscore)),rep("TKO",length(alldf$TKOcscore))), "cscore"=c(alldf$WTcscore,alldf$TKOcscore))

figname <- "willcockson100kb-all-cscores.pdf"
pdf(figname, width = 5, height = 3)
print(ggplot(data = df_reshape, aes(x=cscore, color=sample)) +
  geom_density() +
  xlab("cscore") +
  ylab("frequency") +
  theme_classic(base_size = 14))
dev.off()
