library('ChIPseeker')
library('TxDb.Hsapiens.UCSC.hg38.knownGene')

ATACpeaks <- 

peakAnno <- annotatePeak(ATACpeaks, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)