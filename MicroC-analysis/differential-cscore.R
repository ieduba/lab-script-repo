library(ggplot2)
library(matrixTests)
library(plyr)

correctsign <- function(df){
  scr1rho <- cor.test(df$K36me3score, df$scr1cscore, method = "spearman")
  low1rho <- cor.test(df$K36me3score, df$low1cscore, method = "spearman")
  scr2rho <- cor.test(df$K36me3score, df$scr2cscore, method = "spearman")
  low2rho <- cor.test(df$K36me3score, df$low2cscore, method = "spearman")
  if (scr1rho$estimate < 0){df$scr1cscore <- -df$scr1cscore}
  if (low1rho$estimate < 0){df$low1cscore <- -df$low1cscore}
  if (scr2rho$estimate < 0){df$scr2cscore <- -df$scr2cscore}
  if (low2rho$estimate < 0){df$low2cscore <- -df$low2cscore}
  if (scr1rho$p.value > 0.01 | low1rho$p.value > 0.01 | scr2rho$p.value > 0.01 | low2rho$p.value > 0.01){df$scr1cscore <- "NaN"}
  df
}

datadir <- "/lustre/fs4/risc_lab/scratch/iduba/linker-histone/HiC/KH3/deepseq/compartments/chroms"
setwd(datadir)

alldf <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), c("chr", "start", "end", "row", "K36me3score", "scr1cscore", "scr2cscore", "low1cscore", "low2cscore"))
for (chrom in c(1:22)){
  filename <- paste0("chr", chrom, "-all-K36me3.bed")
  df <- data.frame(read.table(filename, col.names = c("chr", "start", "end", "row", "K36me3score", "scr1cscore", "scr2cscore", "low1cscore", "low2cscore")))
  corrdf <- correctsign(df)
  alldf <- rbind(alldf, corrdf)
}

alldf_nonan <- alldf[alldf$scr1cscore != "NaN",]
dfW <- data.frame("rep1" = as.numeric(alldf_nonan$scr1cscore), "rep2" = alldf_nonan$scr2cscore)
dfL <- data.frame("rep1" = alldf_nonan$low1cscore, "rep2" = alldf_nonan$low2cscore)

tout <- row_t_welch(dfL, dfW)
main <- na.omit(data.frame("binnumber" = rownames(dfW_nonan), "meanWT" = tout$mean.y, "meanlow" = tout$mean.x, "deltac" = tout$mean.diff, "pvalue" = tout$pvalue))
main$logp <- -log10(main$pvalue)

main$shift <- "None"
main$shift[main$pvalue < 0.1 & main$deltac > 0.1 & main$meanWT < -0.1] <- "B to BwA"
main$shift[main$pvalue < 0.1 & main$deltac > 0.1 & main$meanWT > 0.1] <- "A to AwA"
main$shift[main$pvalue < 0.1 & main$deltac > 0.1 & main$meanWT < -0.1 & main$meanlow > 0] <- "B to A"
main$shift[main$pvalue < 0.1 & main$deltac < -0.1 & main$meanWT > 0.1] <- "A to AwB"
main$shift[main$pvalue < 0.1 & main$deltac < -0.1 & main$meanWT < -0.1] <- "B to BwB"
main$shift[main$pvalue < 0.1 & main$deltac < -0.1 & main$meanWT > 0.1 & main$meanlow < 0] <- "A to B"

write.table(main, file="corr-diffcscore.csv", quote = FALSE, row.names = FALSE, sep="\t")

shiftcounts <- count(main$shift)
write.table(shiftcounts, file="corr_shiftcounts.txt", quote=FALSE, row.names = FALSE, sep = "\t")

plotcolors <- c("black", "maroon", "red", "pink", "darkblue", "blue", "lightblue")
names(plotcolors) <- c("None", "B to A", "B to BwA", "A to AwA", "A to B", "A to AwB", "B to BwB")

pdf("corr-cscore-volcano.pdf", width = 5, height = 3)
print(ggplot(data = main, aes(x = deltac, y = logp, color = shift)) + 
  geom_point(size=0.1) +
  scale_color_manual(values = plotcolors) +
  xlab("delta cscore (H1-low - scr)") +
  ylab("-log10(p value)") +
  theme_classic(base_size = 14) +
  guides(colour = guide_legend(override.aes = list(size=2))))
dev.off()
