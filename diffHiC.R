library(multiHiCcompare)
setwd("~/linker-histone/HiC/KH3/deepseq/multiHiCcompare")

reformat <- function(x){
  colorder <- c("chr","region1","region2","IF")
  chrlists = list()
  for (chrom in c(1:22)){
    chromname <- paste0("chr",chrom)
    x[[chromname]]$chr <- chrom
    chrlists <- append(chrlists,x[chromname])
  }
  x[["chrX"]]$chr <- 23
  chrlists <- append(chrlists,x["chrX"])
  flat <- do.call(rbind,chrlists)[,..colorder]
  flat
}

scr100kb1 <- read.table("scr-1-mapped-filt_100kb.txt", header = FALSE)
low100kb1 <- read.table("dH1-1-mapped-filt_100kb.txt", header = FALSE)
scr100kb2 <- read.table("scr-2-mapped-filt_100kb.txt", header = FALSE)
low100kb2 <- read.table("dH1-2-mapped-filt_100kb.txt", header = FALSE)

scrSparse1 <- HiCcompare::cooler2sparse(scr100kb1)
lowSparse1 <- HiCcompare::cooler2sparse(low100kb1)
scrSparse2 <- HiCcompare::cooler2sparse(scr100kb2)
lowSparse2 <- HiCcompare::cooler2sparse(low100kb2)

scrSparseFlat1 <- reformat(scrSparse1)
lowSparseFlat1 <- reformat(lowSparse1)
scrSparseFlat2 <- reformat(scrSparse2)
lowSparseFlat2 <- reformat(lowSparse2)

hicexp <- make_hicexp(scrSparseFlat1, lowSparseFlat1, scrSparseFlat2, lowSparseFlat2, groups = c(0,1,0,1), 
                      zero.p = 0.8, A.min = 5, filter = TRUE, remove.regions = hg38_cyto)

MD_hicexp(hicexp, plot.chr=1, plot.loess = TRUE)
hicexp1 <- cyclic_loess(hicexp, verbose = TRUE, parallel = TRUE, span = 0.2)
MD_hicexp(hicexp1, plot.chr=1, plot.loess=TRUE)

hicexp1 <- hic_exactTest(hicexp1, p.method = 'fdr', parallel = TRUE)
manhattan_hicexp(hicexp1)
