library(ggplot2)

setwd('/ru-auth/local/home/iduba/linker-histone/multi-results/HiCxRNA')

compact <- data.frame(read.table('compaction_RNA-lfc.txt', col.names = "L2FC"))
decompact <- data.frame(read.table('decompaction_RNA-lfc.txt', col.names = "L2FC"))
bigdecompact <- data.frame(read.table('big-decompaction_RNA-lfc.txt', col.names = "L2FC"))

compact$deltacompartment <- 'B-shift'
decompact$deltacompartment <- 'A-shift'
bigdecompact$deltacompartment <- 'A-shift-clust'
df <- rbind(compact,decompact,bigdecompact)

pdf("RNA-lfc_by_compartment-change.pdf", width = 4, height = 6)
print(ggplot(data = df, aes(x=deltacompartment, y=L2FC)) +
        geom_violin() +
        stat_summary(fun.y=median, geom="point", color="red") +
        ylab("RNA-seq log2fold change") +
        xlab("HiC compartment change") +
        theme_classic(base_size = 14))
dev.off()
