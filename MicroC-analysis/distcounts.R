arg = commandArgs(trailingOnly=TRUE)

sample <- arg[1] #name of sample through L1 (before -mapped)
datadir <- arg[2]
nbins <- as.integer(arg[3])
setwd(paste0(datadir))

Infile <- paste0(sample, "-mapped-filt-inward-dist.txt")
Outfile <- paste0(sample, "-mapped-filt-outward-dist.txt")
Tandemfile <- paste0(sample, "-mapped-filt-tandem-dist.txt")
  
orientlist <- c("In", "Out", "Tandem")

for (orient in orientlist){
	orfile <- paste0(orient, "file")
	df <- read.table(get(orfile), sep="\t", col.names = "Distance")
	firstzero <- which.min(abs(df$Distance - 0))
	last1500 <- 'NA'
	last1500 <- which(df$Distance == 1500)[length(which(df$Distance == 1500))]
	if (last1500 == 'NA'){
		last1500 <- which.min(abs(df$Distance - 1500))
	}
	counts <- hist(df$Distance[firstzero:last1500], right=FALSE, breaks = nbins, plot=FALSE)
	countfile <- paste(sample, orient, "counts.txt", sep = '-')
	write.table(counts$counts, file=countfile, quote=FALSE, row.names = FALSE, col.names=FALSE)    
} 
