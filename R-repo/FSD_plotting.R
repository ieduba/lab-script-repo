datadir <- "/Users/irene/Documents/Rockefeller/Risca Lab/fragment size distribution"
setwd(paste0(datadir))

timepoints <- c('CA-S1', 'CB-S10', 'CC-S19', '21A-S7', '21B-S16', '21C-S25')
states <- c(1:15, 'H3K9me3', 'H3K27me3')
statenames <- list('1'='Active TSS', '2'='Flanking active TSS', '3'='Transcr. at gene 5\' and 3\'', '4'='Strong transcription', 
                   '5'='Weak transcription', '6'='Genetic enhancers', '7'='Enhancers', '8'='ZNF genes and repeats', '9'='Heterochromatin', 
                   '10'='Bivalent/Poised TSS', '11'='Flanking bivalent TSS/Enh', '12'='Bivalent enhancer', '13'='Repressed polycomb', 
                   '14'='Weak repressed polycomb', '15'='Quiescent/Low', 'H3K9me3'='H3K9me3', 'H3K27me3'='H3K27me3')


## READING IN RAW DATA TABLES
rawdata <- list()
totalreads <- list()
normdata <- list()
for (tp in timepoints) {
  for (st in states) {
    namedir <- paste(tp, 'DIY-FSD', sep='-')
    namefile <- paste(st, 'DIY_FSD.txt', sep='-')
    datafile <- paste(namedir, namefile, sep='/')
    datatable <- read.table(datafile, col.names=c('Counts', 'Lengths'))
    rawdata[[tp]][[st]] <- datatable
  }
  ## COUNTING TOTAL READS FOR EACH SAMPLE AND NORMALIZING EACH STATE'S COUNTS BY THAT VALUE
  runningsum <- 0
  for (chmm in 1:15){
    runningsum <- runningsum + sum(rawdata[[tp]][[chmm]]$Counts)
  }
  totalreads[[tp]] <- runningsum
  for (st in states){
    normdata[[tp]][[st]] <- rawdata[[tp]][[st]]
    normdata[[tp]][[st]]$'Counts' <- (normdata[[tp]][[st]]$'Counts') / totalreads[[tp]]
  }
}

## LINEAR PLOTTING TRIPLICATES OF CYCLING AND SENESCENT NORMALIZED DATA FOR EACH CHROMATIN STATE
max1 <- list(); min1 <- list()
max2 <- list(); min2 <- list()
for (st in states){
  for (tp in timepoints){
    max1[[st]][[tp]] <- max(normdata[[tp]][[st]]$'Counts')
    min1[[st]][[tp]] <- min(normdata[[tp]][[st]]$'Counts')
  }
  max2[[st]] <- max(max1[[st]])
  min2[[st]] <- min(min1[[st]])

  figurename <- paste(st, 'normalizedFragSizeDist.pdf', sep='_')
  pdf(figurename)
  
  plot(normdata$'CA-S1'[[st]]$Lengths, normdata$'CA-S1'[[st]]$Counts, type="l", xlab='Fragment length', ylab='Normalized frequency', 
       main=statenames[[toString(st)]], ylim = c(0, max2[[st]]), xlim = c(0, 500))
  lines(normdata$'CB-S10'[[st]]$Lengths, normdata$'CB-S10'[[st]]$Counts) 
  lines(normdata$'CC-S19'[[st]]$Lengths, normdata$'CC-S19'[[st]]$Counts)
  lines(normdata$'21A-S7'[[st]]$Lengths, normdata$'21A-S7'[[st]]$Counts, col='red') 
  lines(normdata$'21B-S16'[[st]]$Lengths, normdata$'21B-S16'[[st]]$Counts, col='red')
  lines(normdata$'21C-S25'[[st]]$Lengths, normdata$'21C-S25'[[st]]$Counts, col='red')
  legend("topright", c("Cycling","21 days senescent"), fill=c("black","red"))
  
  dev.off() 
  
  ## LOG PLOTTING TRIPLICATES OF CYCLING AND SENESCENT NORMALIZED DATA FOR EACH CHROMATIN STATE
  figurename <- paste(st, 'normalizedFragSizeDist_log.pdf', sep='_')
  pdf(figurename)
  
  plot(normdata$'CA-S1'[[st]]$Lengths, normdata$'CA-S1'[[st]]$Counts, type="l", xlab='Fragment length', ylab='log normalized frequency', 
       main=statenames[[toString(st)]], log = 'y', ylim = c(min2[[st]], max2[[st]]))
  lines(normdata$'CB-S10'[[st]]$Lengths, normdata$'CB-S10'[[st]]$Counts) 
  lines(normdata$'CC-S19'[[st]]$Lengths, normdata$'CC-S19'[[st]]$Counts)
  lines(normdata$'21A-S7'[[st]]$Lengths, normdata$'21A-S7'[[st]]$Counts, col='red') 
  lines(normdata$'21B-S16'[[st]]$Lengths, normdata$'21B-S16'[[st]]$Counts, col='red')
  lines(normdata$'21C-S25'[[st]]$Lengths, normdata$'21C-S25'[[st]]$Counts, col='red')
  legend("topright", c("Cycling","21 days senescent"), fill=c("black","red"))
  
  dev.off()
}