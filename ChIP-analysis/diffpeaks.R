library('DESeq2')
library('dplyr')
library('ggplot2')

setwd("/ru-auth/local/home/iduba/linker-histone/CutnTag/hg38/diff-peaks")

### FUNCTIONS ###

ggplotMA <- function(data, Title){
  data$Color <- ifelse(data$padj > 0.05, "not significant", ifelse(abs(data$log2FoldChange) > 1, "padj < 0.05 & LFC > 1", "padj < 0.05"))
  pdf(paste0(Title,"-MA.pdf"), width = 8, height = 4)
  print(ggplot(as.data.frame(data), aes(x=baseMean, y=log2FoldChange, color=Color))+
   scale_colour_manual(values=c("grey", "pink", "red"))+
   geom_point(data=as.data.frame(data), aes(x=baseMean, y=log2FoldChange), size=0.2)+ 
   scale_x_continuous(trans='log10')+
   ggtitle(paste(Title, "MAPlot"))+
   xlab("log10 Base Mean")+
   theme_classic(base_size = 18)+
   ylim(c(min(-max(data$log2FoldChange),min(data$log2FoldChange)), max(max(data$log2FoldChange),-min(data$log2FoldChange))))+ #centers y axis around 0
   labs(color="Significance"))
  dev.off()
}

runPCA <- function(dds, Title){
  vsd <- vst(dds, blind=FALSE)
  pcaData <- plotPCA(vsd, intgroup=c("condition", "rep"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pdf(paste0(Title,"-PCA.pdf"), width = 5, height = 4)
  print(ggplot(pcaData, aes(PC1, PC2, color=condition, shape=rep)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed())
  dev.off()
}

markDE <- function(mark, reps){
  peak_count_table = read.table(paste0(mark,"_peak_counts.tsv")) #from diffpeaks-prep.sh

  # convert counts matrix to df based on number of replicates
  if (reps==3){
    colnames(peak_count_table) <- c("Chr","Start","End","low1","low2","low3","scr1","scr2","scr3")
    peak_count_table$key = paste(peak_count_table$Chr, peak_count_table$Start, peak_count_table$End, sep = ':') #this will be used to ID peaks
    input_df <- data.frame(low1=as.integer(peak_count_table$low1), low2=as.integer(peak_count_table$low2), low3=as.integer(peak_count_table$low3), scr1=as.integer(peak_count_table$scr1), scr2=as.integer(peak_count_table$scr2), scr3=as.integer(peak_count_table$scr3))
    rownames(input_df) <- peak_count_table$key
    colData = data.frame(row.names = colnames(input_df), conditions = c("low", "low", "low", "scr", "scr", "scr"), rep = c("1","2","3","1","2","3"))
  }
  if (reps==2.5){ #special case for H1.2 which has three low reps and two scr reps
    colnames(peak_count_table) <- c("Chr","Start","End","low1","low2","low3","scr1","scr2")
    peak_count_table$key = paste(peak_count_table$Chr, peak_count_table$Start, peak_count_table$End, sep = ':')
    input_df <- data.frame(low1=as.integer(peak_count_table$low1), low2=as.integer(peak_count_table$low2), low3=as.integer(peak_count_table$low3), scr1=as.integer(peak_count_table$scr1), scr2=as.integer(peak_count_table$scr2))
    rownames(input_df) <- peak_count_table$key
    colData = data.frame(row.names = colnames(input_df), conditions = c("low", "low", "low", "scr", "scr"), rep = c("1","2","3","1","3"))
  }
  if (reps==2){
    colnames(peak_count_table) <- c("Chr","Start","End","low1","low2","scr1","scr2")
    peak_count_table$key = paste(peak_count_table$Chr, peak_count_table$Start, peak_count_table$End, sep = ':')
    input_df <- data.frame(low1=as.integer(peak_count_table$low1), low2=as.integer(peak_count_table$low2), scr1=as.integer(peak_count_table$scr1), scr2=as.integer(peak_count_table$scr2))
    rownames(input_df) <- peak_count_table$key
    colData = data.frame(row.names = colnames(input_df), conditions = c("low", "low", "scr", "scr"), rep = c("1","2","1","2"))
  }
  colData$condition <- factor(colData$conditions, levels = c("scr", "low"))
  
  # run DESeq on counts matrix based on condition (scr or low)
  dds = DESeqDataSetFromMatrix(countData=input_df, colData= colData, design = ~condition)
  dds = DESeq(dds)
  runPCA(dds, paste0(mark,"-WTvslow-LFC"))
  res = results(dds)
  resLFC = lfcShrink(dds, coef="condition_low_vs_scr", type="normal")
  resLFC.noNA <- na.omit(resLFC)
  ggplotMA(resLFC.noNA, paste0(mark,"-WTvslow-LFC"))

  
  resLFC.df <- as.data.frame(resLFC.noNA)
  resLFC.df$type <- ifelse(resLFC.df$padj < 0.05 & resLFC.df$log2FoldChange > 1, "Up", ifelse(resLFC.df$padj < 0.05 & resLFC.df$log2FoldChange < -1, "Down", "Not significant"))
  write.csv(resLFC.df, paste0(mark,"-WTvslow-peaks.csv"), quote=FALSE, row.names=TRUE)
  
  WTvslowup <- subset(subset(resLFC.df, padj < 0.05), log2FoldChange > 1)
  WTvslowdown <- subset(subset(resLFC.df, padj < 0.05), log2FoldChange < -1)
  write.csv(WTvslowup, paste0(mark,"-WTvslowup.csv"), quote=FALSE, row.names=TRUE)
  write.csv(WTvslowdown, paste0(mark,"-WTvslowdown.csv"), quote=FALSE, row.names=TRUE)
}


### RUN DE ANALYSIS ON HISTONE MARKS ###

markDE("H1p2", 2.5)
markDE("H3k27me3", 3)
markDE("H3k36me2", 3)
markDE("H3k9me3", 3)
markDE("panH1", 3)
markDE("H3k27ac", 2)
markDE("H3k4me3", 2)
markDE("H3k4me1", 2)

