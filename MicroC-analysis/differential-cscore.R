library(ggplot2)
library(matrixTests)
library(plyr)

reslist <- c("100000")
datadir <- "/ru-auth/local/home/iduba/linker-histone/HiC/KH3/deepseq/compartments"
setwd(datadir)

for (res in reslist){
  file1W <- paste0("scr-1_cscore_", res, "-filt.txt")
  file2W <- paste0("scr-2_cscore_", res, "-filt.txt")
  file1L <- paste0("dH1-1_cscore_", res, "-filt.txt")
  file2L <- paste0("dH1-2_cscore_", res, "-filt.txt")
  df1W <- read.table(file1W, col.names = c("Cscore"))
  df2W <- read.table(file2W, col.names = c("Cscore"))
  df1L <- read.table(file1L, col.names = c("Cscore"))
  df2L <- read.table(file2L, col.names = c("Cscore"))
  
  dfW <- data.frame("rep1" = df1W$Cscore, "rep2" = df2W$Cscore)
  dfL <- data.frame("rep1" = df1L$Cscore, "rep2" = df2L$Cscore)
  dfW_nonan <- dfW[dfW$rep1 != "NaN" & dfW$rep2 != "NaN" & dfL$rep1 != "NaN" & dfL$rep2 != "NaN",]
  dfL_nonan <- dfL[dfL$rep1 != "NaN" & dfL$rep2 != "NaN" & dfW$rep1 != "NaN" & dfW$rep2 != "NaN",]
  
  tout <- row_t_welch(dfL_nonan, dfW_nonan)
  main <- na.omit(data.frame("binnumber" = rownames(dfW_nonan), "meanWT" = tout$mean.y, "meanlow" = tout$mean.x, "deltac" = tout$mean.diff, "pvalue" = tout$pvalue))
  main$logp <- -log10(main$pvalue)
  
  main$shift <- "None"
  main$shift[main$pvalue < 0.1 & main$deltac > 0.1 & main$meanWT < -0.1] <- "B to BwA"
  main$shift[main$pvalue < 0.1 & main$deltac > 0.1 & main$meanWT > 0.1] <- "A to AwA"
  main$shift[main$pvalue < 0.1 & main$deltac > 0.1 & main$meanWT < -0.1 & main$meanlow > 0] <- "B to A"
  main$shift[main$pvalue < 0.1 & main$deltac < -0.1 & main$meanWT > 0.1] <- "A to AwB"
  main$shift[main$pvalue < 0.1 & main$deltac < -0.1 & main$meanWT < -0.1] <- "B to BwB"
  main$shift[main$pvalue < 0.1 & main$deltac < -0.1 & main$meanWT > 0.1 & main$meanlow < 0] <- "A to B"
  
  write.table(main, file=paste0(res,"_diffcscore.csv"), quote = FALSE, row.names = FALSE, sep="\t")
  
  shiftcounts <- count(main$shift)
  write.table(shiftcounts, file=paste0(res, "_shiftcounts.txt"), quote=FALSE, row.names = FALSE, sep = "\t")
  
  plotcolors <- c("black", "maroon", "red", "pink", "darkblue", "blue", "lightblue")
  names(plotcolors) <- c("None", "B to A", "B to BwA", "A to AwA", "A to B", "A to AwB", "B to BwB")
 
  figname <- paste0(res, "_cscore-volcano.pdf")
  pdf(figname, width = 5, height = 3)
  print(ggplot(data = main, aes(x = deltac, y = logp, color = shift)) + 
    geom_point(size=0.1) +
    scale_color_manual(values = plotcolors) +
    xlab("delta cscore (H1-low - scr)") +
    ylab("-log10(p value)") +
    theme_classic(base_size = 14) +
    guides(colour = guide_legend(override.aes = list(size=2))))
  dev.off()
  
  main$shift <- "None"
  main$shift[main$deltac > 0.2 & main$meanWT < 0 & main$meanlow < 0] <- "decompaction"
  main$shift[main$deltac > 0.2 & main$meanWT < 0 & main$meanlow > 0] <- "B to A"
  main$shift[main$deltac > 0.2 & main$meanWT > 0 & main$meanlow > 0] <- "decompaction"
  main$shift[main$deltac < -0.2 & main$meanWT > 0 & main$meanlow > 0] <- "compaction"
  main$shift[main$deltac < -0.2 & main$meanWT > 0 & main$meanlow < 0] <- "A to B"
  main$shift[main$deltac < -0.2 & main$meanWT < 0 & main$meanlow < 0] <- "compaction"
  write.table(main, file=paste0(res,"_diffcscore.csv"), quote = FALSE, row.names = FALSE, sep="\t")
  
  shiftcounts <- count(main$shift)
  write.table(shiftcounts, file=paste0(res, "_scatter-shiftcounts.txt"), quote=FALSE, row.names = FALSE, sep = "\t")
  
  plotcolors <- c("black", "maroon", "pink", "darkblue", "lightblue")
  names(plotcolors) <- c("None", "A to B", "compaction", "B to A", "decompaction")
  
  figname <- paste0(res, "_cscore-scatter.pdf")
  pdf(figname, width = 5, height = 3)
  print(ggplot(data = main, aes(x=meanWT, y=meanlow, color = shift)) +
          geom_point(size=0.05) +
          scale_color_manual(values = plotcolors) +
          xlab("scr cscore") +
          ylab("H1-low cscore") +
          theme_classic(base_size = 14) +
          guides(colour = guide_legend(override.aes = list(size=2))))
  dev.off()
  
  df_allnonan <- data.frame("sample"=c(rep("scr-1",length(dfW_nonan$rep1)),rep("scr-2",length(dfW_nonan$rep2)),rep("dH1-1",length(dfL_nonan$rep1)),rep("dH1-2",length(dfL_nonan$rep2))), "cscore"=c(dfW_nonan$rep1,dfW_nonan$rep2,dfL_nonan$rep1,dfL_nonan$rep2))
  
  figname <- paste0(res, "_all-cscores.pdf")
  pdf(figname, width = 5, height = 3)
  print(ggplot(data = df_allnonan, aes(x=cscore, color=sample)) +
    geom_density() +
    xlab("cscore") +
    ylab("frequency") +
    theme_classic(base_size = 14))
  dev.off()
}
