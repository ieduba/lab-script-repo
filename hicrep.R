library(hicrep)
library(ggplot2)

datadir <- "/lustre/fs4/home/iduba/linker-histone/HiC/KH3/deepseq/hicrep"
setwd(datadir)

test1 <- read.table("scr-1-chr1-500kb.matrix")
test2 <- read.table("dH1-1-chr1-500kb.matrix")
h_value <- htrain(mat1, mat2, resol = 500000, lbr = 0, ubr = 5000000, range = 0:10)

scr.scc <- list()
low.scc <- list()
compare1.scc <- list()
compare2.scc <- list()
compare3.scc <- list()
compare4.scc <- list()
for (i in paste0("chr", c(as.character(1:22), 'X'))){
  scr1 = read.table(paste0("scr-1-", i, "-500kb.matrix"))
  scr2 = read.table(paste0("scr-2-", i, "-500kb.matrix"))
  low1 = read.table(paste0("dH1-1-", i, "-500kb.matrix"))
  low2 = read.table(paste0("dH1-2-", i, "-500kb.matrix"))
  scr.scc[[i]] = get.scc(scr1, scr2, 500000, h_value, lbr = 0, ubr = 5000000)
  low.scc[[i]] = get.scc(low1, low2, 500000, h_value, lbr = 0, ubr = 5000000)
  compare1.scc[[i]] = get.scc(scr1, low1, 500000, h_value, lbr = 0, ubr = 5000000)
  compare2.scc[[i]] = get.scc(scr2, low2, 500000, h_value, lbr = 0, ubr = 5000000)
  compare3.scc[[i]] = get.scc(scr1, low2, 500000, h_value, lbr = 0, ubr = 5000000)
  compare4.scc[[i]] = get.scc(scr2, low1, 500000, h_value, lbr = 0, ubr = 5000000)
}

scrsum = 0
lowsum = 0
comp1sum = 0
comp2sum = 0
comp3sum = 0
comp4sum = 0
for (i in paste0("chr", c(as.character(1:22), 'X'))){
  scrsum <- scrsum + scr.scc[[i]]$scc[1,1]
  lowsum <- lowsum + low.scc[[i]]$scc[1,1]
  comp1sum <- comp1sum + compare1.scc[[i]]$scc[1,1]
  comp2sum <- comp2sum + compare2.scc[[i]]$scc[1,1]
  comp3sum <- comp3sum + compare3.scc[[i]]$scc[1,1]
  comp4sum <- comp4sum + compare4.scc[[i]]$scc[1,1]
}
scrmean = scrsum/23
lowmean = lowsum/23
comp1mean = comp1sum/23
comp2mean = comp2sum/23
comp3mean = comp3sum/23
comp4mean = comp4sum/23

toplot <- data.frame(Xvals = c(0.1, 0.2, 0.45, 0.5, 0.55, 0.6), SCC = c(scrmean, lowmean, comp1mean, comp2mean, comp3mean, comp4mean), group = c("scr", "low", "comp", "comp", "comp", "comp"))

pdf("SCC-scatter.pdf", width = 4, height = 4)
print(ggplot(toplot) +
        geom_point(aes(x=Xvals, y=SCC, color=group), size = 2) +
        ylim(0.965, 0.98) +
        ylab("SCC") +
        xlab(NULL) +
        scale_x_continuous(labels = NULL, breaks = NULL) +
        theme_classic())
dev.off()