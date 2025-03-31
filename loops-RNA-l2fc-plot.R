library(ggplot2)

setwd('/ru-auth/local/home/iduba/linker-histone/multi-results/loopxRNA')

scronly <- data.frame(read.table('l2fc-scr-vs-dH1-1-only-anchors.txt', col.names = "l2fc"))
lowonly <- data.frame(read.table('l2fc-scr-vs-dH1-2-only-anchors.txt', col.names = "l2fc"))
common <- data.frame(read.table('l2fc-scr-vs-dH1-common-anchors.txt', col.names = "l2fc"))

scronly$loopset <- 'scr-only'
lowonly$loopset <- 'low-ony'
common$loopset <- 'common'
df <- rbind(scronly,lowonly,common)

pdf("RNA-l2fc-by-loops.pdf", width = 4, height = 6)
print(ggplot(data = df, aes(x=loopset, y=l2fc)) +
        geom_violin() +
        stat_summary(fun.y=median, geom="point", color="red") +
        ylab("RNA-seq log2fold change") +
        xlab("loop set") +
        theme_classic(base_size = 14))
dev.off()
