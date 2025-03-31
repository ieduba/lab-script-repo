library(ggplot2)
library(stringr)
library(sjmisc)

## functions for finding max and min of contact prob curves, used for NRL and odd/even calculations
maxes <- function(d){
  i <- which(diff(sign(diff(d$y)))<0)+1
  data.frame(X = d$x[i], Y = d$y[i], Orientation = 0, Sample = 0)
}
mins <- function(d){
  i <- which(diff(sign(diff(d$y)))>0)+1
  data.frame(X = d$x[i], Y = d$y[i])
}


## list of samples you want to plot, without the suffix (specified below)
fignum <- "hg19" #for naming the plots that are generated
#samplelist <- c("dH1-12-combo", "dH1-34-combo", "scr-12-combo", "scr-34-combo")
#samplelist <- c("scr-12-combo-Het","scr-34-combo-Het","dH1-12-combo-Het","dH1-34-combo-Het","scr-12-combo-EnhA","scr-34-combo-EnhA","dH1-12-combo-EnhA","dH1-34-combo-EnhA","scr-12-combo-Tx","scr-34-combo-Tx","dH1-12-combo-Tx","dH1-34-combo-Tx","scr-12-combo-TssA","scr-34-combo-TssA","dH1-12-combo-TssA","dH1-34-combo-TssA", "scr-12-combo-ReprPC","scr-34-combo-ReprPC","dH1-12-combo-ReprPC","dH1-34-combo-ReprPC")
samplelist <- c("scr-12-combo-K562-H3K27ac", "scr-34-combo-K562-H3K27ac", "dH1-12-combo-K562-H3K27ac", "dH1-34-combo-K562-H3K27ac","scr-12-combo-K562-H3K9me3", "scr-34-combo-K562-H3K9me3", "dH1-12-combo-K562-H3K9me3", "dH1-34-combo-K562-H3K9me3","scr-12-combo-K562-H3K27me3", "scr-34-combo-K562-H3K27me3", "dH1-12-combo-K562-H3K27me3", "dH1-34-combo-K562-H3K27me3")
#samplelist <- c("dH1-12-combo-WTvslowup-mean", "dH1-34-combo-WTvslowup-mean", "scr-12-combo-WTvslowup-mean", "scr-34-combo-WTvslowup-mean", "dH1-12-combo-WTvslow-nochange-basemeanLFC-meandownsample", "dH1-34-combo-WTvslow-nochange-basemeanLFC-meandownsample", "scr-12-combo-WTvslow-nochange-basemeanLFC-meandownsample", "scr-34-combo-WTvslow-nochange-basemeanLFC-meandownsample")
#samplelist <- c("scr-12-combo-K562_scr1_q1", "scr-34-combo-K562_scr1_q1", "dH1-12-combo-K562_low1_q1", "dH1-34-combo-K562_low1_q1","scr-12-combo-K562_scr1_q4sample", "scr-34-combo-K562_scr1_q4sample", "dH1-12-combo-K562_low1_q4sample", "dH1-34-combo-K562_low1_q4sample")
#samplelist <- c("mESC-WT-ctl1", "mESC-WT-ctl2", "mESC-WT-slowed", "mESC-low-ctl1", "mESC-low-ctl2", "mESC-low-slowed")

datadir <- "/ru-auth/local/home/iduba/linker-histone/Micro-C/in-situ/old-new-combo/sep-bio-reps/contactProb"
#datadir <- "/ru-auth/local/home/iduba/linker-histone/RNA-seq/mouse/enlongation-genes-microc/contactProb"
setwd(datadir)

## initialize dataframes
alldf <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Distance", "Orientation","Sample", "CellType", "Region", "Replicate"))
allcounts <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("Breaks", "Density", "Orientation","Sample", "CellType", "Region", "Replicate"))
allmaxes <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("X", "Y", "Orientation", "Sample", "Peak", "CellType", "Region", "Replicate"))
allslopes <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("N", "dY", "Orientation", "Sample", "CellType", "Region", "Replicate"))
allpeakinfo <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("NRL", "Orientation", "Sample", "CellType", "Region", "Replicate", "LongShort", "OddEven"))

for (sample in samplelist){
  ## parse the sample names defined above. you will likely need to edit this so that all the variables are correctly assigned
  namebits <- str_split(sample, "-", 5, simplify=TRUE)
  region <- namebits[5]
  celltype <- namebits[1]
  replicate <- namebits[2]
  name <- paste(region, celltype, replicate, sep='-')
 
  ## suffix for files, might change. item in sample list + this should be your file name.
  Infile <- paste0(sample, "-mapped-filt-inward-dist.txt")
  Outfile <- paste0(sample, "-mapped-filt-outward-dist.txt")
  Tandemfile <- paste0(sample, "-mapped-filt-tandemplus-dist.txt")
  orientlist <- c("Out")
  
  for (orient in orientlist){
    orfile <- paste0(orient, "file")
    if(orient=="In"){toload <- -1} #imports all
    #else{toload <- 8000000} #use this if importing Out/Tandem takes too long (cuts off the data after 8M rows)
    else{toload <- -1}
    df <- read.table(get(orfile), sep="\t", col.names = "Distance", nrows = toload)
    df$Orientation <- orient
    df$Sample <- name
    df$Region <- region
    df$CellType <- celltype
    df$Replicate <- replicate

   ## find the first row number with value 10 and last row number with value 1500 to set the range of rows for plotting distances 10-1500 bp
   firstten <- which.min(abs(df$Distance - 10))
   last1500 <- which(df$Distance == 1500)[length(which(df$Distance == 1500))]
   if (is_empty(last1500) == TRUE){
     last1500 <- which.min(abs(df$Distance - 1500))
   }
   
   ## make smoothed histogram of distances from 10-1500 bp for calculating maxes, fill out max dataframe
   d <- density(df$Distance[firstten:last1500])
   samplemaxes <- maxes(d)
   samplemaxes$Orientation <- orient
   samplemaxes$Sample <- name
   samplemaxes$Region <- region
   samplemaxes$CellType <- celltype
   samplemaxes$Peak <- as.numeric(row.names(samplemaxes))
   samplemaxes$Replicate <- replicate
   
   ## calculate slopes between maxes, fill out slope dataframe
   sampleslopes <- setNames(data.frame(matrix(ncol = 2, nrow = length(samplemaxes$Y)-1)), c("N","dY"))
   for (i in 2:length(samplemaxes$Peak)){
     sampleslopes$N[i-1] <- i
     sampleslopes$dY[i-1] <- samplemaxes$Y[i]-samplemaxes$Y[i-1]
   }
   sampleslopes$Orientation <- orient
   sampleslopes$Sample <- name
   sampleslopes$Region <- region
   sampleslopes$CellType <- celltype
   sampleslopes$Replicate <- replicate   
   
   ## calculate NRL as mean difference between maxes, fill out NRL dataframe
   nrl <- mean(diff(samplemaxes$X, 1))
   samplepeakinfo <- data.frame(NRL=nrl)
   samplepeakinfo$Orientation <- orient
   samplepeakinfo$Sample <- name
   samplepeakinfo$CellType <- celltype
   samplepeakinfo$Region <- region
   samplepeakinfo$Replicate <- replicate
   
   ## save max locations as txt file
   maxfilename <- paste(name, orient, "maxes.txt", sep = "-")
   write.table(samplemaxes, file=maxfilename)

   ## calculate mins of of smoothed contact prob curve. use mins as dividers between peaks for calculating peak sums
   samplemins <- mins(d)
   sums <- setNames(data.frame(matrix(ncol = 3, nrow = length(samplemins$X))), c("start", "end","sum"))
   for (r in seq(length(samplemins$X))){
     #samplemins$Ind[r] <- which.min(abs(df$Distance - samplemins$X[r]))
     if (r == 1){
       sums$start[r] <- 0
     } else {
         sums$start[r] <- which(d$x == samplemins$X[r-1])
     }
     sums$end[r] <- which(d$x == samplemins$X[r])
     sums$sum[r] <- sum(d$y[sums$start[r]:sums$end[r]]) # sums up the smoothed contact prob between peak minima
   }
   longshort <- sum(sums$sum[4:6])/sum(sums$sum[1:3]) # ratio of long:short peak sums, defined as peaks 4-6/1-3
   oddeven <- sum(sums$sum[3],sums$sum[5])/sum(sums$sum[2],sums$sum[4]) # ratio of odd:even peak sums, defined as peaks 3,5/2,4
   
   ## fill out peak info dataframe
   samplepeakinfo$LongShort <- longshort
   samplepeakinfo$OddEven <- oddeven
   samplepeakinfo <- setNames(samplepeakinfo, c("NRL", "Orientation", "Sample", "CellType", "Region", "Replicate", "LongShort", "OddEven"))
   
   histdata <- hist(df$Distance[firstten:last1500], breaks = 500, plot = FALSE) #raw counts data for histogram plotting with 500 bins
   countsdf <- data.frame(Breaks=histdata$breaks[1:length(histdata$density)], Density=histdata$density) #convert from counts to density
   countsdf$Orientation = orient
   countsdf$Sample = name
   countsdf$CellType = celltype
   countsdf$Region = region
   countsdf$Replicate = replicate
   
    if (celltype == "scr"){
      mycolor <- "#00BFC4"
    } else {
      mycolor <- "#F8766D"
    }
    
   ## set colors based on whatever variable you're interested in (usually region)
   if (region == "H3K27ac"){
     mycolor <- "#E68613"
   } else if (region == "H3K27me3"){
     mycolor <- "#7CAE00"
   } else if (region == "H3K9me3"){
     mycolor <- "#C77CFF"
   }
   
#   if (region == "Enh"){
#     mycolor <- "#FFC425"
#   } else if (region == "Het"){
#     mycolor <- "#FF61CC"
#   } else if (region == "TssA"){
#     mycolor <- "#8494FF"
#   } else if (region == "Tx"){
#     mycolor <- "#66A182"
#   }
   
     ## plot of individual sample histogram (usually not that useful)
#    fig1name <- paste(sample, orient, "contactprob.pdf", sep = "-")
#    pdf(fig1name, width = 5, height = 3)
#    print(ggplot(df, aes(x=Distance)) +
#            xlim(c(0,1500)) +
#            geom_histogram(aes(y=..density..), binwidth=10) +
#            geom_line(aes(y= ..density..),stat = 'density', size=1, color=mycolor) +
#            xlab("Genomic distance (bp)")+
#            geom_point(data=samplemaxes, mapping= aes(x=X, y=Y), color=mycolor, size=3, alpha = 0.9) +
#            ylab("Contact probability") +
#            theme_classic(base_size = 20))
#    dev.off()
#    
    ## add all the sample data to the master dataframes for summary plots
    alldf <- rbind(alldf, df)
    allcounts <- rbind(allcounts, countsdf)
    allmaxes <- rbind(allmaxes, samplemaxes)
    allslopes <- rbind(allslopes, sampleslopes)
    allpeakinfo <- rbind(allpeakinfo, as.data.frame(samplepeakinfo))   
  }
}

## SUMMARY PLOTS -- MAY NEED TO CHANGE OUTPUT FILE NAMES
  ## plots replicates of WT and low histograms, maxes, slopes for each region
  for (orient in "Out"){
    for (region in unique(alldf$Region)){
        
        pdf(paste(orient, region, "WTvslow-maxes.pdf", sep = '-'), width = 6, height = 3)
        print(ggplot(subset(allmaxes, Orientation == orient & Region == region)) +
                geom_line(mapping= aes(x=Peak, y=Y, color=Sample)) +
                geom_point(mapping= aes(x=Peak, y=Y, color=Sample), size = 1, alpha = 0.9) +
                ylim(0.0003,0.0028) +
                ylab("Peak probability") +
                xlab("Peak number") +
                scale_x_continuous(breaks=c(2,4,6,8)) +
                scale_color_manual(values=c("#F8766D","#F8766D","#00BFC4","#00BFC4")) +
                theme_classic())
        dev.off()
        
        pdf(paste(orient, region, "WTvslow-slopes.pdf", sep = '-'), width = 6, height = 3)
        print(ggplot(subset(allslopes, Orientation == orient & Region == region)) +
                geom_line(mapping= aes(x=N, y=dY, color=Sample)) +
                geom_point(mapping= aes(x=N, y=dY, color=Sample), size = 1) +
                ylab("Peak slope") +
                xlab("Peak number") +
                scale_x_continuous(breaks=c(2,4,6,8)) +
                scale_color_manual(values=c("#F8766D","#F8766D","#00BFC4","#00BFC4")) +
                theme_classic())
        dev.off()
        
        pdf(paste(orient, region, "WTvslow-hist.pdf", sep = '-'), width = 6, height = 3)
        print(ggplot(subset(allcounts, Orientation == orient & Region == region)) +
                geom_line(mapping= aes(x=Breaks, y=Density, color=Sample), alpha=0.8) +
                ylab("Density") +
                xlab("Distance (bp)") +
                scale_color_manual(values=c("#F8766D","#F8766D","#00BFC4","#00BFC4")) +
                theme_classic())
        dev.off()
    }
    
    ## plots all regions on the same plot, for WT and H1-low
    for (celltype in unique(alldf$CellType)){
     
      pdf(paste(orient, celltype, "hg19-marks-maxes.pdf", sep = '-'), width = 6, height = 3)
      print(ggplot(subset(allmaxes, Orientation == orient & CellType == celltype)) +
        geom_line(mapping= aes(x=Peak, y=Y, color=Sample)) +
        geom_point(mapping= aes(x=Peak, y=Y, color=Sample), size = 1, alpha = 0.9) +
        ylim(0.0003,0.0028) +
        ylab("Peak probability") +
        xlab("Peak number") +
        scale_x_continuous(breaks=c(2,4,6,8)) +
        theme_classic() +
        #scale_color_manual(values=c("#FFC425","#FFC425","#FF61CC","#FF61CC","#8494FF","#8494FF","#66A182","#66A182","#E68000","#E68000")))
        scale_color_manual(values=c("#E68613","#E68613","#7CAE00","#7CAE00","#C77CFF","#C77CFF")))
      dev.off()
      
      pdf(paste(orient, celltype, "hg19-marks-slopes.pdf", sep = '-'), width = 6, height = 3)
      print(ggplot(subset(allslopes, Orientation == orient & CellType == celltype)) +
              geom_line(mapping= aes(x=N, y=dY, color=Sample)) +
              geom_point(mapping= aes(x=N, y=dY, color=Sample), size = 1, alpha = 0.9) +
              ylab("Peak slope") +
              xlab("Peak number") +
              scale_x_continuous(breaks=c(2,4,6,8)) +
              theme_classic() +
              #scale_color_manual(values=c("#FFC425","#FFC425","#FF61CC","#FF61CC","#8494FF","#8494FF","#66A182","#66A182","#E68000","#E68000")))
              scale_color_manual(values=c("#E68613","#E68613","#7CAE00","#7CAE00","#C77CFF","#C77CFF")))
      dev.off()
      
      pdf(paste(orient, celltype, "hg19-marks-hist.pdf", sep = '-'), width = 6, height = 3)
      print(ggplot(subset(allcounts, Orientation == orient & CellType == celltype)) +
        geom_line(mapping= aes(x=Breaks, y=Density, color=Sample), alpha=0.8) +
        ylab("Density") +
        xlab("Distance (bp)") +
        theme_classic() +
        #scale_color_manual(values=c("#FFC425","#FFC425","#FF61CC","#FF61CC","#8494FF","#8494FF","#66A182","#66A182","#E68000","#E68000")))
        scale_color_manual(values=c("#E68613","#E68613","#7CAE00","#7CAE00","#C77CFF","#C77CFF")))
      dev.off()

   }
  }

## save peak info as summary txt file
summaryfilename <- paste(fignum, "summary.txt", sep = "-")
write.table(allpeakinfo, file=summaryfilename)

## prep for bar plot of summary quantifications
formean <- subset(allpeakinfo, Orientation =='Out') #use only outward ligations for summary stats
formean$ctreg <- paste(formean$CellType, formean$Region, sep='-') #combine cell type and region as a new variable so it can be aggregated over

## find average and standard devation of NRL, short vs long, and odd vs even for outward ligation averaged across replicates for each cell type/region combo
meansdnrl <- aggregate(formean$NRL, by=list(formean$ctreg), FUN=function(y) c(M=mean(y,na.rm=TRUE), SD=sd(y)))
meansdls <- aggregate(formean$LongShort, by=list(formean$ctreg), FUN=function(y) c(M=mean(y,na.rm=TRUE), SD=sd(y)))
meansdoe <- aggregate(formean$OddEven, by=list(formean$ctreg), FUN=function(y) c(M=mean(y,na.rm=TRUE), SD=sd(y)))

## bar charts with average quantifications, error bars are +/- standard deviation
fig5name <- paste(fignum,'allNRL.pdf',sep='-')
pdf(fig5name, width = 3, height = 3)
print(ggplot(data = meansdnrl, aes(x=Group.1, y=x[,"M"])) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(aes(ymin=x[,"M"]-x[,"SD"], ymax=x[,"M"]+x[,"SD"]), width = 0.3) +
  coord_cartesian(ylim=c(160,200)) +
  ylab('NRL (bp)') +
  xlab('Sample') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8)))
dev.off()

fig6name <- paste(fignum,'alllongshort.pdf',sep='-')
pdf(fig6name,width = 3, height = 3)
print(ggplot(data = meansdls, aes(x=Group.1, y=x[,"M"])) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(aes(ymin=x[,"M"]-x[,"SD"], ymax=x[,"M"]+x[,"SD"]), width = 0.3) +
  coord_cartesian(ylim=c(0.3,0.7)) +
  ylab('Long/Short contacts') +
  xlab('Sample')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8)))
dev.off()

fig7name <- paste(fignum,'alloddeven.pdf',sep='-')
pdf(fig7name,width = 3, height = 3)
print(ggplot(data = meansdoe, aes(x=Group.1, y=x[,"M"])) +
        geom_bar(stat = "identity", width = 0.8) +
        geom_errorbar(aes(ymin=x[,"M"]-x[,"SD"], ymax=x[,"M"]+x[,"SD"]), width = 0.3) +
        coord_cartesian(ylim=c(0.85, 1.15)) +
        ylab('Odd/Even contacts') +
        xlab('Sample') +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8)))
dev.off()
