library(ggplot2)
library(stringr)
#library(sjmisc)

maxes <- function(d){
  i <- which(diff(sign(diff(d$y)))<0)+1
  data.frame(X = d$x[i], Y = d$y[i], Orientation = 0, Sample = 0)
}

samplelist <- c("31_q1_gene-mESC_low_CKDL210025175-1a-4_HVWV3DSX2_L1", "18_q1_gene-mESC_WT_CKDL210025175-1a-2_HVWV3DSX2_L1", "18_q4_gene-mESC_WT_CKDL210025175-1a-2_HVWV3DSX2_L1", "31_q4_gene-mESC_low_CKDL210025175-1a-4_HVWV3DSX2_L1")
datadir <- "contactProb/"
setwd(paste0(datadir))

alldf <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Distance", "Orientation","Sample", "CellType", "Region"))
allmaxes <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("X", "Y", "Orientation", "Peak", "CellType", "Region"))
allnrls <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("NRL", "Orientation", "CellType", "Region"))

for (sample in samplelist){
  #name <- sample
  namebits <- str_split(sample, "-", 2, simplify = TRUE)
  name <- namebits[1]
  #name <- paste(namebits[1], namebits[2], sep = "-")
  #name <- paste(str_split(namebits[1], "_", 3, simplify = TRUE)[1], str_split(namebits[1], "_", 3, simplify = TRUE)[2], sep = "_")
  region <- str_split(name, "_", 3, simplify = TRUE)[2]
  #region <- "Genome-wide"
  celltype <- str_split(name, "_", 3, simplify = TRUE)[1]
 # celltype <- name
  #replicate <- str_split(name, "_", 3, simplify = TRUE)[3]
 
  Infile <- paste0(sample, "-mapped-filt-inward-dist.txt")
  Outfile <- paste0(sample, "-mapped-filt-outward-dist.txt")
  Tandemfile <- paste0(sample, "-mapped-filt-tandemplus-dist.txt")
  
  orientlist <- c("Tandem", "In", "Out")
  
  for (orient in orientlist){
    orfile <- paste0(orient, "file")
    if(orient=="In"){toload <- -1} #imports all
    else{toload <- 5000000}
    df <- read.table(get(orfile), sep="\t", col.names = "Distance", nrows = toload)
    df$Orientation <- orient
    df$Sample <- name
    df$Region <- region
    df$CellType <- celltype
  #  df$Replicate <- replicate

   firstzero <- which.min(abs(df$Distance - 0))
   last1500 <- which(df$Distance == 1500)[length(which(df$Distance == 1500))]
#   if (is_empty(last1500) == TRUE){
#     last1500 <- which.min(abs(df$Distance - 1500))
#   }
   d <- density(df$Distance[firstzero:last1500])
   samplemaxes <- maxes(d)
   samplemaxes$Orientation <- orient
   samplemaxes$Sample <- name
   samplemaxes$Region <- region
   samplemaxes$CellType <- celltype
   samplemaxes$Peak <- as.numeric(row.names(samplemaxes))
  # samplemaxes$Replicate <- replicate
 
   nrl <- mean(diff(samplemaxes$X, 1))
   samplenrl <- nrl
   samplenrl$Orientation <- orient
   samplenrl$CellType <- celltype
   samplenrl$Region <- region
   samplenrl <- setNames(samplenrl, c("NRL", "Orientation", "CellType", "Region"))
  
   maxfilename <- paste(name, orient, "maxes.txt", sep = "-")
   write.table(samplemaxes, file=maxfilename)
    
    if (str_split(name,"_",2,simplify=TRUE)[2] == "WT"){
      mycolor <- "#00BFC4"
    } else {
      mycolor <- "#F8766D"
    }
    
#    fig1name <- paste(name, orient, "contactprob.pdf", sep = "-")
#    pdf(fig1name, width = 5, height = 3)
#    print(ggplot(df, aes(x=Distance)) +
#            xlim(c(0,1500)) +
#            geom_histogram(aes(y=..density..), binwidth=10) +
##            geom_histogram(aes(y=..density..), bins = 500) +
##            scale_x_log10(limits= c(1,150000)) +
#            geom_line(aes(y= ..density..),stat = 'density', size=1, color=mycolor) +
#            xlab("Genomic distance (bp)")+
#            geom_point(data=samplemaxes, mapping= aes(x=X, y=Y), color=mycolor, size=3, alpha = 0.9) +
#            ylab("Contact probability") +
#            theme_classic(base_size = 20))
#    dev.off()
    
    alldf <- rbind(alldf, df)
    allmaxes <- rbind(allmaxes, samplemaxes)
    allnrls <- rbind(allnrls, samplenrl)   
  }
  
# fig2name <- paste0(name, "-all.pdf")
# pdf(fig2name, width = 10, height = 6)
# print(ggplot(subset(alldf, Sample=name), aes(x=Distance, color=Orientation)) +
#         geom_density() +
#         xlim(0,1500) +
#         ggtitle(name))
# dev.off()
 
}

for (region in unique(alldf$Region)){
  for (orient in unique(alldf$Orientation)){
    fig3name <- paste(orient, region, "18vs31-CP.pdf", sep = '-')
    pdf(fig3name, width = 3.5, height = 1.5)
    print(ggplot(subset(alldf, Orientation == orient & Region == region), aes(x=Distance, color=Sample)) +
      geom_density() +
      xlim(0,1500) +
      xlab("Distance (bp)") +
      ylab("Density"))
      scale_x_continuous(trans = 'log10') +
    dev.off()
    
    fig4name <- paste(orient, region, "18vs31-maxes.pdf", sep = '-')
    pdf(fig4name, width = 3.5, height = 1.5)
    print(ggplot(subset(allmaxes, Orientation == orient && Region == region)) +
      geom_line(mapping= aes(x=Peak, y=Y, color=Sample)) +
      geom_point(mapping= aes(x=Peak, y=Y, color=Sample), size = 2, alpha = 0.9) +
      ylim(0,0.002) +
      ylab("Peak probability") +
      xlab("Peak number") +
      theme_classic())
    dev.off()
  }
}

