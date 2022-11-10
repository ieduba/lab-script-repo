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

fignum <- "mito-opt1-QC"
samplelist <- c("WT-mito_opt1", "WT-async_opt1", "low-mito_opt1", "low-async_opt1", "WT-mito_opt2", "WT-async_opt2", "low-mito_opt2", "low-async_opt2")
#samplelist <- c("K562_WT_CKDL220004209-1a-12_HG5T2DSX3_L2", "K562_low_CKDL220004209-1a-19_HG5T2DSX3_L2", "K562_WT_CKDL210025175-1a-5_HVWV3DSX2_L1", "K562_low_CKDL210025175-1a-6_HVWV3DSX2_L1")
#samplelist <- c("mESC_low_CKDL210025175-1a-4_HVWV3DSX2_L1", "mESC_WT_CKDL210025175-1a-2_HVWV3DSX2_L1", "mESC_low_CKDL220004209-1a-8_HG5T2DSX3_L2", "mESC_WT_CKDL220004209-1a-7_HG5T2DSX3_L2")
#samplelist <- c("eu-ENCFF163QUM-hg38-K562_low_CKDL210025175-1a-6_HVWV3DSX2_L1", "eu-ENCFF163QUM-hg38-K562_low_CKDL220004209-1a-19_HG5T2DSX3_L2","eu-ENCFF163QUM-hg38-K562_WT_CKDL210025175-1a-5_HVWV3DSX2_L1","eu-ENCFF163QUM-hg38-K562_WT_CKDL220004209-1a-12_HG5T2DSX3_L2", "het-ENCFF163QUM-hg38-K562_low_CKDL210025175-1a-6_HVWV3DSX2_L1", "het-ENCFF163QUM-hg38-K562_low_CKDL220004209-1a-19_HG5T2DSX3_L2", "het-ENCFF163QUM-hg38-K562_WT_CKDL210025175-1a-5_HVWV3DSX2_L1", "het-ENCFF163QUM-hg38-K562_WT_CKDL220004209-1a-12_HG5T2DSX3_L2")
#samplelist <- c("eu-K562_low_CKDL210025175-1a-6_HVWV3DSX2_L1", "eu-K562_WT_CKDL210025175-1a-5_HVWV3DSX2_L1","het-K562_low_CKDL210025175-1a-6_HVWV3DSX2_L1", "het-K562_WT_CKDL210025175-1a-5_HVWV3DSX2_L1")
#samplelist <- c("eu-mESC_WT_CKDL210025175-1a-2_HVWV3DSX2_L1","het-mESC_WT_CKDL210025175-1a-2_HVWV3DSX2_L1","eu-mESC_low_CKDL210025175-1a-4_HVWV3DSX2_L1", "het-mESC_low_CKDL210025175-1a-4_HVWV3DSX2_L1","eu-mESC_WT_CKDL220004209-1a-7_HG5T2DSX3_L2","het-mESC_WT_CKDL220004209-1a-7_HG5T2DSX3_L2", "eu-mESC_low_CKDL220004209-1a-8_HG5T2DSX3_L2", "het-mESC_low_CKDL220004209-1a-8_HG5T2DSX3_L2" )
#samplelist <- c("Enh-hg38-K562_low_CKDL210025175-1a-6_HVWV3DSX2_L1", "Enh-hg38-K562_WT_CKDL210025175-1a-5_HVWV3DSX2_L1", "Het-hg38-K562_low_CKDL210025175-1a-6_HVWV3DSX2_L1", "Het-hg38-K562_WT_CKDL210025175-1a-5_HVWV3DSX2_L1", "TssA-hg38-K562_low_CKDL210025175-1a-6_HVWV3DSX2_L1", "TssA-hg38-K562_WT_CKDL210025175-1a-5_HVWV3DSX2_L1", "Enh-hg38-K562_low_CKDL220004209-1a-19_HG5T2DSX3_L2", "Enh-hg38-K562_WT_CKDL220004209-1a-12_HG5T2DSX3_L2", "Het-hg38-K562_low_CKDL220004209-1a-19_HG5T2DSX3_L2", "Het-hg38-K562_WT_CKDL220004209-1a-12_HG5T2DSX3_L2", "TssA-hg38-K562_low_CKDL220004209-1a-19_HG5T2DSX3_L2", "TssA-hg38-K562_WT_CKDL220004209-1a-12_HG5T2DSX3_L2")
#samplelist <- c("Enh-mm10-mESC_low_CKDL210025175-1a-4_HVWV3DSX2_L1", "Enh-mm10-mESC_WT_CKDL210025175-1a-2_HVWV3DSX2_L1", "Het-mm10-mESC_low_CKDL210025175-1a-4_HVWV3DSX2_L1", "Het-mm10-mESC_WT_CKDL210025175-1a-2_HVWV3DSX2_L1", "Pro-mm10-mESC_low_CKDL210025175-1a-4_HVWV3DSX2_L1", "Pro-mm10-mESC_WT_CKDL210025175-1a-2_HVWV3DSX2_L1", "Enh-mm10-mESC_low_CKDL220004209-1a-8_HG5T2DSX3_L2", "Enh-mm10-mESC_WT_CKDL220004209-1a-7_HG5T2DSX3_L2", "Het-mm10-mESC_low_CKDL220004209-1a-8_HG5T2DSX3_L2","Het-mm10-mESC_WT_CKDL220004209-1a-7_HG5T2DSX3_L2","Pro-mm10-mESC_low_CKDL220004209-1a-8_HG5T2DSX3_L2","Pro-mm10-mESC_WT_CKDL220004209-1a-7_HG5T2DSX3_L2")
#samplelist <- c("18_q1_gene-mESC_WT_CKDL220004209-1a-7_HG5T2DSX3_L2", "18_q4_gene-mESC_WT_CKDL220004209-1a-7_HG5T2DSX3_L2", "18_q1_gene-mESC_WT_CKDL210025175-1a-2_HVWV3DSX2_L1", "18_q4_gene-mESC_WT_CKDL210025175-1a-2_HVWV3DSX2_L1", "31_q1_gene-mESC_low_CKDL220004209-1a-8_HG5T2DSX3_L2", "31_q4_gene-mESC_low_CKDL220004209-1a-8_HG5T2DSX3_L2", "31_q1_gene-mESC_low_CKDL210025175-1a-4_HVWV3DSX2_L1","31_q4_gene-mESC_low_CKDL210025175-1a-4_HVWV3DSX2_L1")
#samplelist <- c("18to31up_gene-mESC_WT_CKDL210025175-1a-2_HVWV3DSX2_L1", "18to31up_gene-mESC_low_CKDL210025175-1a-4_HVWV3DSX2_L1", "18to31down_gene-mESC_WT_CKDL210025175-1a-2_HVWV3DSX2_L1","18to31down_gene-mESC_low_CKDL210025175-1a-4_HVWV3DSX2_L1", "18to31up_gene-mESC_WT_CKDL220004209-1a-7_HG5T2DSX3_L2", "18to31up_gene-mESC_low_CKDL220004209-1a-8_HG5T2DSX3_L2", "18to31down_gene-mESC_WT_CKDL220004209-1a-7_HG5T2DSX3_L2", "18to31down_gene-mESC_low_CKDL220004209-1a-8_HG5T2DSX3_L2" )

#datadir <- "/Users/ireneduba/Dropbox (Dropbox @RU)/Risca Laboratory/Risca Laboratory/Users/Irene/_Work_In_Progress/1 - Linker histone project/Micro-C/dovetail/tech rep 1/deep sequencing/mESC/eu-het"
datadir <- "/Users/ireneduba/Dropbox (Dropbox @RU)/Risca Laboratory/Risca Laboratory/Users/Irene/_Work_In_Progress/1 - Linker histone project/Micro-C/mitosis/QC/compare"
setwd(paste0(datadir))

alldf <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Distance", "Orientation","Sample", "CellType", "CellCycle", "Replicate"))
allmaxes <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("X", "Y", "Orientation", "Peak", "CellType", "CellCycle", "Replicate"))
allpeakinfo <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("NRL", "Orientation", "Sample", "CellType", "CellCycle", "Replicate", "LongShort", "OddEven"))

for (sample in samplelist){
  namebits <- str_split(sample, "_",3,simplify=TRUE)
  name <- namebits[1]
  celltype <- str_split(name,"-",2,simplify=TRUE)[1]
  cellcycle <- str_split(name,"-",2,simplify=TRUE)[2]
  replicate <- namebits[2]
  
  #name <- sample
  #namebits <- str_split(sample, "-", 4, simplify = TRUE)
  #name <- namebits[1]
  #name <- paste(namebits[1], namebits[2], sep = "-")
  #name <- paste(str_split(namebits[2], '_', 3, simplify = TRUE)[1], str_split(namebits[2], '_', 3, simplify = TRUE)[2], namebits[1], sep = "_")
  #name <- paste(str_split(namebits[1], "_", 3, simplify = TRUE)[1], str_split(namebits[1], "_", 3, simplify = TRUE)[2], sep = "_")
  #name <- paste(str_split(namebits[1],'_',2,simplify=TRUE)[1],str_split(namebits[2],'_',3,simplify=TRUE)[1],str_split(namebits[2],'_',3,simplify=TRUE)[2],sep='_')
  #region <- str_split(name, "_", 2, simplify = TRUE)[2]
  #region <- "Genome-wide"
  #region <- namebits[1]
  #celltype <- paste(str_split(namebits[2], "_", 3, simplify = TRUE)[1],str_split(namebits[2], "_", 3, simplify = TRUE)[2],sep='_')
  #celltype <- str_split(name, "_", 2, simplify = TRUE)[2]
  #celltype <- name
  #replicate <- str_split(namebits[4], "_", 3, simplify = TRUE)[3]
 
  Infile <- paste0(sample, "-mapped-filt-inward-dist.txt")
  Outfile <- paste0(sample, "-mapped-filt-outward-dist.txt")
  Tandemfile <- paste0(sample, "-mapped-filt-tandemplus-dist.txt")
  #Tandemfile <- paste0(sample, "-mapped.pairs-tandem-adj-dist.txt")
  #Infile <- paste0(sample, "-mapped.pairs-in-adj-dist.txt")
  #Outfile <- paste0(sample, "-mapped.pairs-out-adj-dist.txt")  
  orientlist <- c("Tandem", "Out")
  
  for (orient in orientlist){
    orfile <- paste0(orient, "file")
    #if(orient=="In"){toload <- -1} #imports all
   # else{toload <- 3000000}
    toload <- 8000000
    df <- read.table(get(orfile), sep="\t", col.names = "Distance", nrows = toload)
    df$Orientation <- orient
    df$Sample <- name
    df$CellCycle <- cellcycle
    df$CellType <- celltype
    df$Replicate <- replicate

   firstzero <- which.min(abs(df$Distance - 0))
   last1500 <- which(df$Distance == 1500)[length(which(df$Distance == 1500))]
   if (is_empty(last1500) == TRUE){
     last1500 <- which.min(abs(df$Distance - 1500))
   }
   d <- density(df$Distance[firstzero:last1500])
   samplemaxes <- maxes(d)
   samplemaxes$Orientation <- orient
   samplemaxes$Sample <- name
   samplemaxes$CellCycle <- cellcycle
   samplemaxes$CellType <- celltype
   samplemaxes$Peak <- as.numeric(row.names(samplemaxes))
   samplemaxes$Replicate <- replicate
   
   nrl <- mean(diff(samplemaxes$X, 1))
   samplepeakinfo <- nrl
   samplepeakinfo$Orientation <- orient
   samplepeakinfo$Sample <- name
   samplepeakinfo$CellType <- celltype
   samplepeakinfo$CellCycle <- cellcycle
   samplepeakinfo$Replicate <- replicate
   
   maxfilename <- paste(name, orient, replicate, "maxes.txt", sep = "-")
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
   
    if (str_split(name,"_",2,simplify=TRUE)[2] == "WT"){
      mycolor <- "#00BFC4"
    } else {
      mycolor <- "#F8766D"
    }
    
    fig1name <- paste(name, orient, replicate, "contactprob.pdf", sep = "-")
    pdf(fig1name, width = 5, height = 3)
    print(ggplot(df, aes(x=Distance)) +
            xlim(c(0,1500)) +
            geom_histogram(aes(y=..density..), binwidth=10) +
#            geom_histogram(aes(y=..density..), bins = 500) +
#            scale_x_log10(limits= c(1,150000)) +
            geom_line(aes(y= ..density..),stat = 'density', size=1, color=mycolor) +
            xlab("Genomic distance (bp)")+
            geom_point(data=samplemaxes, mapping= aes(x=X, y=Y), color=mycolor, size=3, alpha = 0.9) +
            ylab("Contact probability") +
            theme_classic(base_size = 20))
    dev.off()
    
    alldf <- rbind(alldf, df)
    allmaxes <- rbind(allmaxes, samplemaxes)
    allpeakinfo <- rbind(allpeakinfo, samplepeakinfo)   
  }
  
# fig2name <- paste0(name, "-all.pdf")
# pdf(fig2name, width = 10, height = 6)
# print(ggplot(subset(alldf, Sample=name), aes(x=Distance, color=Orientation)) +
#         geom_density() +
#         xlim(0,1500) +
#         ggtitle(name))
# dev.off()
 
}

for (celltype in unique(alldf$CellType)){
  for (orient in unique(alldf$Orientation)){
    for (cellcycle in unique(alldf$CellCycle)){
      fig3name <- paste(orient, celltype, cellcycle, "digestion-hist.pdf", sep = '-')
      pdf(fig3name, width = 4, height = 2.5)
      print(ggplot(subset(alldf, Orientation == orient & CellType == celltype & CellCycle == cellcycle), aes(x=Distance, fill=Replicate, color=Replicate)) +
              geom_histogram(aes(y=..density..), binwidth=20, position = "identity", alpha=0.5) +
              xlim(0,1500) +
              xlab("Distance (bp)") +
              ylab("Density") +
              theme_classic())
      dev.off()
      
      fig4name <- paste(orient, celltype, cellcycle, "digestion-maxes.pdf", sep = '-')
      pdf(fig4name, width = 4, height = 2.5)
      print(ggplot(subset(allmaxes, Orientation == orient & CellType == celltype & CellCycle == cellcycle)) +
              geom_line(mapping= aes(x=Peak, y=Y, color=Replicate)) +
              geom_point(mapping= aes(x=Peak, y=Y, color=Replicate), size = 2, alpha = 0.9) +
              ylim(0,0.002) +
              ylab("Peak probability") +
              xlab("Peak number") +
              theme_classic())
      dev.off()
    }
    for (replicate in unique(alldf$Replicate)){
      fig3name <- paste(orient, celltype, replicate, "cellcycle-hist.pdf", sep = '-')
      pdf(fig3name, width = 4, height = 2.5)
      print(ggplot(subset(alldf, Orientation == orient & CellType == celltype & Replicate == replicate), aes(x=Distance, fill=CellType, color=CellCycle)) +
              geom_histogram(aes(y=..density..), binwidth=20, position = "dodge", alpha=0.5) +
              xlim(0,1500) +
              xlab("Distance (bp)") +
              ylab("Density") +
              theme_classic())
      dev.off()
      
      fig4name <- paste(orient, celltype, replicate, "cellcycle-maxes.pdf", sep = '-')
      pdf(fig4name, width = 4, height = 2.5)
      print(ggplot(subset(allmaxes, Orientation == orient & CellType == celltype & Replicate == replicate)) +
              geom_line(mapping= aes(x=Peak, y=Y, color=CellCycle)) +
              geom_point(mapping= aes(x=Peak, y=Y, color=CellCycle), size = 2, alpha = 0.9) +
              ylim(0,0.002) +
              ylab("Peak probability") +
              xlab("Peak number") +
              theme_classic())
      dev.off()
    }
  }
}





for (replicate in unique(alldf$Replicate)){
  for (orient in unique(alldf$Orientation)){
    for (celltype in unique(alldf$CellType)){
      fig3name <- paste(orient, celltype, replicate, "regions-CP.pdf", sep = '-')
      pdf(fig3name, width = 4, height = 2.5)
      print(ggplot(subset(alldf, Orientation == orient & CellType == celltype & Replicate == replicate), aes(x=Distance, color=Region)) +
        geom_density() +
        xlim(0,1500) +
        xlab("Distance (bp)") +
        ylab("Density") +
        theme_classic())
      dev.off()
      
      fig4name <- paste(orient, celltype, replicate, "regions-maxes.pdf", sep = '-')
      pdf(fig4name, width = 4, height = 2.5)
      print(ggplot(subset(allmaxes, Orientation == orient & CellType == celltype & Replicate == replicate)) +
        geom_line(mapping= aes(x=Peak, y=Y, color=Region)) +
        geom_point(mapping= aes(x=Peak, y=Y, color=Region), size = 2, alpha = 0.9) +
        ylim(0,0.002) +
        ylab("Peak probability") +
        xlab("Peak number") +
        theme_classic())
      dev.off()
    }
   for (region in unique(alldf$Region)){
     fig3name <- paste(orient, region, replicate, "WTvsLow-CP.pdf", sep = '-')
     pdf(fig3name, width = 4, height = 2.5)
     print(ggplot(subset(alldf, Orientation == orient & Region == region & Replicate == replicate), aes(x=Distance, color=CellType)) +
             geom_density() +
             xlim(0,1500) +
             xlab("Distance (bp)") +
             ylab("Density") +
             theme_classic())
     dev.off()
     
     fig4name <- paste(orient, region, replicate, "WTvsLow-maxes.pdf", sep = '-')
     pdf(fig4name, width = 4, height = 2.5)
     print(ggplot(subset(allmaxes, Orientation == orient & Region == region & Replicate == replicate)) +
             geom_line(mapping= aes(x=Peak, y=Y, color=CellType)) +
             geom_point(mapping= aes(x=Peak, y=Y, color=CellType), size = 2, alpha = 0.9) +
             ylim(0,0.002) +
             ylab("Peak probability") +
             xlab("Peak number") +
             theme_classic())
     dev.off()
   }
  }
}

formean <- subset(allpeakinfo, Orientation != 'In')

meannrl <- aggregate(formean$NRL,by=list(formean$Sample, formean$CellType, formean$Region, formean$Replicate), FUN=mean)
meansdnrl <- aggregate(meannrl$x, by=list(meannrl$Group.1, meannrl$Group.2, meannrl$Group.3), FUN=function(y) c(M=mean(y), SD=sd(y)))

meanlongshort <- aggregate(formean$LongShort,by=list(formean$Sample, formean$CellType, formean$Region, formean$Replicate), FUN=mean)
meansdls <- aggregate(meanlongshort$x, by=list(meanlongshort$Group.1, meanlongshort$Group.2, meanlongshort$Group.3), FUN=function(y) c(M=mean(y), SD=sd(y)))

meanoddeven <- aggregate(formean$OddEven,by=list(formean$Sample, formean$CellType, formean$Region, formean$Replicate), FUN=mean)
meansdoe <- aggregate(meanoddeven$x, by=list(meanoddeven$Group.1, meanoddeven$Group.2, meanoddeven$Group.3), FUN=function(y) c(M=mean(y), SD=sd(y)))

fig5name <- paste(fignum,'allNRL.pdf',sep='-')
pdf(fig5name, width = 3, height = 2)
print(ggplot(data = meansdnrl, aes(x=Group.1, y=x[,"M"])) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(aes(ymin=x[,"M"]-x[,"SD"], ymax=x[,"M"]+x[,"SD"]), width = 0.3) +
  coord_cartesian(ylim=c(160,200)) +
  ylab('NRL (bp)') +
  xlab('Sample') +
  theme_classic())
dev.off()

fig6name <- paste(fignum,'alllongshort.pdf',sep='-')
pdf(fig6name,width = 3, height = 2)
print(ggplot(data = meansdls, aes(x=Group.1, y=x[,"M"])) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(aes(ymin=x[,"M"]-x[,"SD"], ymax=x[,"M"]+x[,"SD"]), width = 0.3) +
  coord_cartesian(ylim=c(0.6,1)) +
  ylab('Long/Short contacts') +
  xlab('Sample')+
  theme_classic())
dev.off()

fig7name <- paste(fignum,'alloddeven.pdf',sep='-')
pdf(fig7name,width = 3, height = 2)
print(ggplot(data = meansdoe, aes(x=Group.1, y=x[,"M"])) +
        geom_bar(stat = "identity", width = 0.8) +
        geom_errorbar(aes(ymin=x[,"M"]-x[,"SD"], ymax=x[,"M"]+x[,"SD"]), width = 0.3) +
        coord_cartesian(ylim=c(0.8, 1.2)) +
        ylab('Odd/Even contacts') +
        xlab('Sample') +
        theme_classic())
dev.off()
