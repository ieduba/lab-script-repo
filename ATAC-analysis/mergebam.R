library(bedr)
library(Rsamtools)

setwd("/Users/ireneduba/Documents/Rockefeller/Risca Lab/Senescence project/LSA1-5combo")

for (expt in c("LSA2","LSA3","LSA")) {
  print(paste("experiment:", expt, sep=' '))
  for (timepoint in c("C",0,3,4,6,9,12,14,16,21,28)) {
    print paste("time point:" timepoint, sep=' ')
    name <- paste(timepoint, expt, sep='-')
    mergedname <- paste(name, "_merged.bam", sep = '')
    group <- grep(name, list.files(), value = TRUE)
    print(paste("   ", group))
    mergeBam(group, mergedname)
  }
}
