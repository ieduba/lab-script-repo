library(ggplot2)

setwd('/ru-auth/local/home/iduba/linker-histone/multi-results/loopxRNA')

scronly <- data.frame(read.table('basemeans-scr-vs-dH1-1-only-anchors.txt', col.names = "basemean"))
lowonly <- data.frame(read.table('basemeans-scr-vs-dH1-2-only-anchors.txt', col.names = "basemean"))
common <- data.frame(read.table('basemeans-scr-vs-dH1-common-anchors.txt', col.names = "basemean"))

scronly$loopset <- 'scr-only'
lowonly$loopset <- 'low-ony'
common$loopset <- 'common'
df <- rbind(scronly,lowonly,common)

pdf("RNA-basemean-by-loops.pdf", width = 4, height = 6)
print(ggplot(data = df, aes(x=loopset, y=basemean)) +
        geom_violin() +
        stat_summary(fun.y=median, geom="point", color="red") +
	scale_y_continuous(trans='log10') +
        ylab("base mean transcription (log)") +
        xlab("loop set") +
        theme_classic(base_size = 14))
dev.off()
