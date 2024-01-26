library(ggplot2)
library(stringr)
library(sjmisc)

maxes <- function(d){
  i <- which(diff(sign(diff(d$y)))<0)+1
  data.frame(X = d$x[i], Y = d$y[i], Orientation = 0, Sample = 0)
}

mins <- function(d){
  i <- which(diff(sign(diff(d$y)))>0)+1
  data.frame(X = d$x[i], Y = d$y[i])
}

fignum <- "gp-chip"
#samplelist <- c("scr-12-Het","scr-34-Het","dH1-12-Het","dH1-34-Het","scr-12-Enh","scr-34-Enh","dH1-12-Enh","dH1-34-Enh","scr-12-Tx","scr-34-Tx","dH1-12-Tx","dH1-34-Tx","scr-12-TssA","scr-34-TssA","dH1-12-TssA","dH1-34-TssA")
samplelist <- c("scr-12-K562-H3K27ac", "scr-34-K562-H3K27ac", "dH1-12-K562-H3K27ac", "dH1-34-K562-H3K27ac","scr-12-K562-H3K9me3", "scr-34-K562-H3K9me3", "dH1-12-K562-H3K9me3", "dH1-34-K562-H3K9me3","scr-12-K562-H3K27me3", "scr-34-K562-H3K27me3", "dH1-12-K562-H3K27me3", "dH1-34-K562-H3K27me3")
datadir <- "/lustre/fs4/risc_lab/scratch/iduba/linker-histone/Micro-C/in-situ/contact-prob-compare"
setwd(datadir)

alldf <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Distance", "Orientation","Sample", "CellType", "Region", "Replicate"))
allcounts <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("Breaks", "Density", "Orientation","Sample", "CellType", "Region", "Replicate"))
allmaxes <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("X", "Y", "Orientation", "Peak", "CellType", "Region", "Replicate"))
allpeakinfo <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("NRL", "Orientation", "Sample", "CellType", "Region", "Replicate", "LongShort", "OddEven"))

for (sample in samplelist){
  namebits <- str_split(sample, "-", 4, simplify=TRUE)
  region <- namebits[4]
  celltype <- namebits[1]
  replicate <- namebits[2]
  name <- paste(region, celltype, replicate, sep='-')
 
  Infile <- paste0(sample, "-gp-mapped-filt-inward-dist.txt")
  Outfile <- paste0(sample, "-gp-mapped-filt-outward-dist.txt")
  Tandemfile <- paste0(sample, "-gp-mapped-filt-tandemplus-dist.txt")
  orientlist <- c("Tandem", "Out")
  
  for (orient in orientlist){
    orfile <- paste0(orient, "file")
    if(orient=="In"){toload <- -1} #imports all
    #else{toload <- 8000000}
    else{toload <- -1}
    df <- read.table(get(orfile), sep="\t", col.names = "Distance", nrows = toload)
    df$Orientation <- orient
    df$Sample <- name
    df$Region <- region
    df$CellType <- celltype
    df$Replicate <- replicate

   firstten <- which.min(abs(df$Distance - 10))
   last1500 <- which(df$Distance == 1500)[length(which(df$Distance == 1500))]
   if (is_empty(last1500) == TRUE){
     last1500 <- which.min(abs(df$Distance - 1500))
   }
   d <- density(df$Distance[firstten:last1500])
   samplemaxes <- maxes(d)
   samplemaxes$Orientation <- orient
   samplemaxes$Sample <- name
   samplemaxes$Region <- region
   samplemaxes$CellType <- celltype
   samplemaxes$Peak <- as.numeric(row.names(samplemaxes))
   samplemaxes$Replicate <- replicate
   
   nrl <- mean(diff(samplemaxes$X, 1))
   samplepeakinfo <- data.frame(NRL=nrl)
   samplepeakinfo$Orientation <- orient
   samplepeakinfo$Sample <- name
   samplepeakinfo$CellType <- celltype
   samplepeakinfo$Region <- region
   samplepeakinfo$Replicate <- replicate
   
   maxfilename <- paste(name, orient, "maxes.txt", sep = "-")
   write.table(samplemaxes, file=maxfilename)

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
     sums$sum[r] <- sum(d$y[sums$start[r]:sums$end[r]])
   }
   longshort <- sum(sums$sum[4:6])/sum(sums$sum[1:3])
   oddeven <- sum(sums$sum[3],sums$sum[5])/sum(sums$sum[2],sums$sum[4])
   samplepeakinfo$LongShort <- longshort
   samplepeakinfo$OddEven <- oddeven
   samplepeakinfo <- setNames(samplepeakinfo, c("NRL", "Orientation", "Sample", "CellType", "Region", "Replicate", "LongShort", "OddEven"))
   
   histdata <- hist(df$Distance[firstten:last1500], breaks = 500, plot = FALSE)
   countsdf <- data.frame(Breaks=histdata$breaks[1:length(histdata$density)], Density=histdata$density)
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
   
    fig1name <- paste(sample, orient, "contactprob.pdf", sep = "-")
    pdf(fig1name, width = 5, height = 3)
    print(ggplot(df, aes(x=Distance)) +
            xlim(c(0,1500)) +
            geom_histogram(aes(y=..density..), binwidth=10) +
            geom_line(aes(y= ..density..),stat = 'density', size=1, color=mycolor) +
            xlab("Genomic distance (bp)")+
            geom_point(data=samplemaxes, mapping= aes(x=X, y=Y), color=mycolor, size=3, alpha = 0.9) +
            ylab("Contact probability") +
            theme_classic(base_size = 20))
    dev.off()
    
    alldf <- rbind(alldf, df)
    allcounts <- rbind(allcounts, countsdf)
    allmaxes <- rbind(allmaxes, samplemaxes)
    allpeakinfo <- rbind(allpeakinfo, as.data.frame(samplepeakinfo))   
  }
}

#for (replicate in unique(alldf$Replicate)){
  for (orient in "Out"){
    for (region in unique(alldf$Region)){
        
#        fig3name <- paste(orient, region, "WTvslow-CP.pdf", sep = '-')
#        pdf(fig3name, width = 6, height = 2.5)
#        print(ggplot(subset(alldf, Orientation == orient & Region == region), aes(x=Distance, color=Sample)) +
#                geom_density() +
#                xlim(0,1500) +
#                xlab("Distance (bp)") +
#                ylab("Density") +
#                scale_color_manual(values=c("#F8766D","#F8766D","#00BFC4","#00BFC4")) +
#                theme_classic())
#        dev.off()
        
        fig4name <- paste(orient, region, "WTvslow-maxes.pdf", sep = '-')
        pdf(fig4name, width = 6, height = 2.5)
        print(ggplot(subset(allmaxes, Orientation == orient & Region == region)) +
                geom_line(mapping= aes(x=Peak, y=Y, color=Sample)) +
                geom_point(mapping= aes(x=Peak, y=Y, color=Sample), size = 1, alpha = 0.9) +
                ylim(0.0003,0.0028) +
                ylab("Peak probability") +
                xlab("Peak number") +
                scale_x_continuous(breaks=c(2,4,6)) +
                scale_color_manual(values=c("#F8766D","#F8766D","#00BFC4","#00BFC4")) +
                theme_classic())
        dev.off()
        
        fig4name <- paste(orient, region, "WTvslow-hist.pdf", sep = '-')
        pdf(fig4name, width = 6, height = 2.5)
        print(ggplot(subset(allcounts, Orientation == orient & Region == region)) +
                geom_line(mapping= aes(x=Breaks, y=Density, color=Sample), alpha=0.8) +
                ylab("Density") +
                xlab("Distance (bp)") +
                scale_color_manual(values=c("#F8766D","#F8766D","#00BFC4","#00BFC4")) +
                theme_classic())
        dev.off()
    }
    for (celltype in unique(alldf$CellType)){
 #     fig3name <- paste(orient, celltype, "hmmregions-CP.pdf", sep = '-')
 #     pdf(fig3name, width = 6, height = 3)
 #     print(ggplot(subset(alldf, Orientation == orient & CellType == celltype), aes(x=Distance, color=Sample)) +
 #       geom_density() +
 #       xlim(0,1500) +
 #       xlab("Distance (bp)") +
 #       ylab("Density") +
 #       scale_color_manual(values=c("#E68613","#E68613","#7CAE00","#7CAE00","#C77CFF","#C77CFF")) +
 #       theme_classic())
 #     dev.off()
      
      fig4name <- paste(orient, celltype, "regions-maxes.pdf", sep = '-')
      pdf(fig4name, width = 6, height = 3)
      print(ggplot(subset(allmaxes, Orientation == orient & CellType == celltype)) +
        geom_line(mapping= aes(x=Peak, y=Y, color=Sample)) +
        geom_point(mapping= aes(x=Peak, y=Y, color=Sample), size = 1, alpha = 0.9) +
        ylim(0.0003,0.0025) +
        ylab("Peak probability") +
        xlab("Peak number") +
        scale_x_continuous(breaks=c(2,4,6,8)) +
        theme_classic() +
        #scale_color_manual(values=c("#FFC425","#FFC425","#FF61CC","#FF61CC","#8494FF","#8494FF","#66A182","#66A182")))
        scale_color_manual(values=c("#E68613","#E68613","#7CAE00","#7CAE00","#C77CFF","#C77CFF")))
      dev.off()
      
      fig4name <- paste(orient, celltype, "regions-hist.pdf", sep = '-')
      pdf(fig4name, width = 6, height = 3)
      print(ggplot(subset(allcounts, Orientation == orient & CellType == celltype)) +
        geom_line(mapping= aes(x=Breaks, y=Density, color=Sample), alpha=0.8) +
        ylab("Density") +
        xlab("Distance (bp)") +
        theme_classic() +
        #scale_color_manual(values=c("#FFC425","#FFC425","#FF61CC","#FF61CC","#8494FF","#8494FF","#66A182","#66A182")))
        scale_color_manual(values=c("#E68613","#E68613","#7CAE00","#7CAE00","#C77CFF","#C77CFF")))
      dev.off()

   }
  }
#}

summaryfilename <- paste(fignum, "summary.txt", sep = "-")
write.table(allpeakinfo, file=summaryfilename)

formean <- subset(allpeakinfo, Orientation =='Out')
formean$ctreg <- paste(formean$CellType, formean$Region, sep='-')

#Averages across orientations, then averages and SDs across reps. should do this all in one step to calculate SD for all
meansdnrl <- aggregate(formean$NRL, by=list(formean$ctreg), FUN=function(y) c(M=mean(y,na.rm=TRUE), SD=sd(y)))
meansdls <- aggregate(formean$LongShort, by=list(formean$ctreg), FUN=function(y) c(M=mean(y,na.rm=TRUE), SD=sd(y)))
meansdoe <- aggregate(formean$OddEven, by=list(formean$ctreg), FUN=function(y) c(M=mean(y,na.rm=TRUE), SD=sd(y)))

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
