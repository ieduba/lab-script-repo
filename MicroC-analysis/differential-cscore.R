library(ggplot2)
library(matrixTests)
library(plyr)

celllist <- c("K562", "mESC")
datadir <- "/Users/ireneduba/Dropbox (Dropbox @RU)/Risca Laboratory/Risca Laboratory/Users/Irene/_Work_In_Progress/1 - Linker histone project/mESC/Micro-C/dovetail/tech rep 1/deep sequencing/cscore"
setwd(paste0(datadir))

for (cell in celllist){
  file1W <- paste0(cell, "-WT-1-100kb_cscore.txt")
  file2W <- paste0(cell, "-WT-2-100kb_cscore.txt")
  file1L <- paste0(cell, "-low-1-100kb_cscore.txt")
  file2L <- paste0(cell, "-low-2-100kb_cscore.txt")
  df1W <- read.table(file1W, sep="\t", col.names = c("Bin number", "Cscore"))
  df2W <- read.table(file2W, sep="\t", col.names = c("Bin number", "Cscore"))
  df1L <- read.table(file1L, sep="\t", col.names = c("Bin number", "Cscore"))
  df2L <- read.table(file2L, sep="\t", col.names = c("Bin number", "Cscore"))
  
  dfW <- data.frame("rep1" = df1W$Cscore, "rep2" = df2W$Cscore)
  dfL <- data.frame("rep1" = df1L$Cscore, "rep2" = df2L$Cscore)

  tout <- row_t_welch(dfW, dfL)
  main <- na.omit(data.frame("meanWT" = tout$mean.x, "meanlow" = tout$mean.y, "deltac" = tout$mean.diff, "pvalue" = tout$pvalue))
  main$logp <- -log10(main$pvalue)
  
  main$shift <- "None"
  main$shift[main$pvalue < 0.1 & main$deltac > 0.1 & main$meanWT < -0.1 & main$meanlow > 0] <- "B to A"
  main$shift[main$pvalue < 0.1 & main$deltac > 0.1 & main$meanWT < -0.1 & main$meanlow < 0] <- "B to BwA"
  main$shift[main$pvalue < 0.1 & main$deltac > 0.1 & main$meanWT > 0.1] <- "A to AwA"
  main$shift[main$pvalue < 0.1 & main$deltac < -0.1 & main$meanWT > 0.1 & main$meanlow < 0] <- "A to B"
  main$shift[main$pvalue < 0.1 & main$deltac < -0.1 & main$meanWT > 0.1 & main$meanlow > 0] <- "A to AwB"
  main$shift[main$pvalue < 0.1 & main$deltac < -0.1 & main$meanWT < -0.1] <- "B to BwB"
  
  shiftcounts <- count(main$shift)
  write.table(shiftcounts, file=paste0(cell, "_shiftcounts.txt"), quote=FALSE, row.names = FALSE, sep = "\t")
  
  plotcolors <- c("black", "maroon", "red", "pink", "darkblue", "blue", "lightblue")
  names(plotcolors) <- c("None", "B to A", "B to BwA", "A to AwA", "A to B", "A to AwB", "B to BwB")
  
  figname <- paste0(cell, "_cscore-volcano.pdf")
  pdf(figname, width = 4, height = 3)
  print(ggplot(data = main, aes(x = deltac, y = logp, color = shift)) + 
    geom_point(size=0.5) +
    scale_color_manual(values = plotcolors) +
    xlab("delta cscore (WT - H1 low)") +
    ylab("log10(p value)") +
    theme_classic(base_size = 14) +
    guides(colour = guide_legend(override.aes = list(size=2))))
  dev.off()
}
