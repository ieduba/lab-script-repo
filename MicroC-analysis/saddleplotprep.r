#usage:Rscript saddleplotprep.r [merged_nodups.txt] [cscore.bedgraph] [resolution of cscore bed in bp]
args = commandArgs(trailingOnly=TRUE)
 
#function for rounding up positions in contact list to a given resolution
roundup <- function(pos, res){round(pos + (res/2), -log(res, 10))}

#read in contact list and cscore bedgraph
contacts <- read.table(args[1], sep = "\t", col.names = c("Str1", "Chr1", "Pos1", "Frag1", "Str2", "Chr2", "Pos2", "Frag2", "mapq1", "cigar1", "sequence1", "mapq2", "cigar2", "sequence2", "name1", "name2"))
cscore <- read.table(args[2], sep = "\t", skip = 1, col.names = c("Chr", "Start", "End", "Cscore"))
res <- args[3]

#create LUT of cscores for given chrs and positions
getcscore <- cscore$Cscore
names(getcscore) <- paste(cscore$Chr, cscore$End, sep = "-")

#assign contacts' cscores, order by cscore, reestablish .summary format
contacts$tomatch <- paste(contacts$Chr1, roundup(contacts$Pos1, 10000), sep = "-")
contacts$Cscore <- unname(getcscore[contacts$tomatch])
ordcontacts <- contacts[order(contacts$Cscore),]
saddledata <- ordcontacts[,1:8]

write.table(saddledata, sep = "\t", file = paste0("saddle-", args[1]), quote = FALSE, row.names = FALSE, col.names = FALSE) 
