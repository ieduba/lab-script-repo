library(ggplot2)

setwd('/ru-auth/local/home/iduba/linker-histone/multi-results/HiCxATAC/l2fc')

compact <- data.frame(read.table('compaction_all-ATAC-lfc.txt', col.names = "L2FC"))
decompact <- data.frame(read.table('decompaction_all-ATAC-lfc.txt', col.names = "L2FC"))
bigdecompact <- data.frame(read.table('bigdecompaction_all-ATAC-lfc.txt', col.names = "L2FC"))

compact$deltacompartment <- 'B-shift'
decompact$deltacompartment <- 'A-shift'
bigdecompact$deltacompartment <- 'A-shift-clust'
df <- rbind(compact,decompact,bigdecompact)

pdf("allATAC-lfc_by_compartment-change.pdf", width = 4, height = 6)
print(ggplot(data = df, aes(x=deltacompartment, y=L2FC)) +
        geom_violin() +
        stat_summary(fun.y=median, geom="point", color="red") +
	ylab("ATAC-seq log2fold change") +
        xlab("HiC compartment change") +
        theme_classic(base_size = 14))
dev.off()
