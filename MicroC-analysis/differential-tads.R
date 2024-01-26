library(ggplot2)
library(matrixTests)
library(plyr)

cell <- "K562"
datadir <- "~/linker-histone/Micro-C/in-situ/KIM1-2-combo/tad-compare"
setwd(datadir)


file1W <- "scr-12-cornerscores.txt"
file2W <- "scr-34-cornerscores.txt"
file1L <- "dH1-12-cornerscores.txt"
file2L <- "dH1-34-cornerscores.txt"
df1W <- read.table(file1W, col.names = c("cornerscore"))
df2W <- read.table(file2W, col.names = c("cornerscore"))
df1L <- read.table(file1L, col.names = c("cornerscore"))
df2L <- read.table(file2L, col.names = c("cornerscore"))

dfW <- data.frame("rep1" = df1W$cornerscore, "rep2" = df2W$cornerscore)
dfL <- data.frame("rep1" = df1L$cornerscore, "rep2" = df2L$cornerscore)
dfW_nonan <- dfW[dfW$rep1 != "NaN" & dfW$rep2 != "NaN" & dfL$rep1 != "NaN" & dfL$rep2 != "NaN",]
dfL_nonan <- dfL[dfL$rep1 != "NaN" & dfL$rep2 != "NaN" & dfW$rep1 != "NaN" & dfW$rep2 != "NaN",]

tout <- row_t_welch(dfW_nonan, dfL_nonan)
main <- na.omit(data.frame("binnumber" = rownames(dfW_nonan), "meanWT" = tout$mean.x, "meanlow" = tout$mean.y, "deltac" = tout$mean.diff, "pvalue" = tout$pvalue))
main$logp <- -log10(main$pvalue)

main$shift <- "None"
main$shift[main$pvalue < 0.1 & main$deltac > 0.05] <- "weakening"
main$shift[main$pvalue < 0.1 & main$deltac < -0.05] <- "strengthening"
write.table(main, file="diffcornerscore.csv", quote = FALSE, row.names = FALSE, sep="\t")

shiftcounts <- count(main$shift)
write.table(shiftcounts, file="shiftcounts.txt", quote=FALSE, row.names = FALSE, sep = "\t")

plotcolors <- c("black", "red", "blue")
names(plotcolors) <- c("None", "weakening", "strengthening" )

figname <- "TAD-volcano.pdf"
pdf(figname, width = 5, height = 3)
print(ggplot(data = main, aes(x = deltac, y = logp, color = shift)) + 
  geom_point(size=0.1) +
  scale_color_manual(values = plotcolors) +
  xlab("delta corner score (H1 low - scr)") +
  ylab("-log10(p value)") +
  theme_classic(base_size = 14) +
  guides(colour = guide_legend(override.aes = list(size=2))))
dev.off()

