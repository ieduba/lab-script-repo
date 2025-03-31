library(ggplot2)
datadir <- "/ru-auth/local/home/iduba/linker-histone/RNA-seq/human/DESeq2_results"
setwd(datadir)

df <- read.table("WTvslow-fixed-mean-genebody.bed", sep="\t", col.names = c("chr","pos1","pos2","gene","mean","l2fc"))
df$kblength <- (df$pos2 - df$pos1)/1000

pdf("gene-length-vs-change.pdf", width = 4, height = 4)
print(ggplot(df, aes(x=kblength, y=l2fc)) +
        geom_point(size=0.7) +
        xlab("gene length (kb)") +
        ylab("log2fold change") +
        theme_classic())
dev.off()
