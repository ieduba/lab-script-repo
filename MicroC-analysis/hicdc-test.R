library('HiCDCPlus')
library('data.table')
library('BSgenome.Hsapiens.UCSC.hg38')
library('DESeq2')

setDTthreads(threads = 1)
setwd('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/Micro-C/in-situ/hicdc')

#add .hic counts
hicfile_paths<-c('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/Micro-C/in-situ/KIM1-2-combo/dH1-12/dH1-12-mapped-filt-mq30.hic',
                 '/lustre/fs4/risc_lab/scratch/iduba/linker-histone/Micro-C/in-situ/KIM1-2-combo/scr-12/scr-12-mapped-filt-mq30.hic',
                 '/lustre/fs4/risc_lab/scratch/iduba/linker-histone/Micro-C/in-situ/KIM3-4-combo/dH1-34/dH1-34-mapped-filt-mq30.hic',
                 '/lustre/fs4/risc_lab/scratch/iduba/linker-histone/Micro-C/in-situ/KIM3-4-combo/scr-34/scr-34-mapped-filt-mq30.hic'
                 )

indexfile<-data.frame()
construct_features(output_path=getwd(), gen="Hsapiens", gen_ver="hg38", sig="NULL", bin_type="Bins-uniform", binsize=50000, chrs=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"), feature_type = "RE-agnostic")

for(hicfile_path in hicfile_paths){
  output_path<-paste0(getwd(),'/', gsub("^(.*[\\/])", "",gsub('.hic','.txt.gz',hicfile_path)))
  #generate gi_list instance
  gi_list<-generate_bintolen_gi_list(paste0(getwd(),'_bintolen.txt.gz'),gen_ver="hg38")
  gi_list<-add_hic_counts(gi_list,hic_path = hicfile_path)
  #expand features for modeling
  gi_list<-expand_1D_features(gi_list)
  #run HiC-DC+
  set.seed(1010) #HiC-DC downsamples rows for modeling
  gi_list<-HiCDCPlus(gi_list,ssize=0.1)
  for (i in seq(length(gi_list))){
    indexfile<-unique(rbind(indexfile,
                            as.data.frame(gi_list[[i]][gi_list[[i]]$qvalue<=0.05])[c('seqnames1',
                                                                                     'start1','start2')]))
  }
  #write results to a text file
  gi_list_write(gi_list,fname=output_path,chrs=c('chr1','chr2'))
}

#save index file---union of significants at 50kb
colnames(indexfile)<-c('chr','startI','startJ')
data.table::fwrite(indexfile,
                   paste0(getwd(),'/hicdc_significant_indices.txt.gz'),
                   sep='\t',row.names=FALSE,quote=FALSE)

#Differential analysis using modified DESeq2 
hicdcdiff(input_paths=list(scr=c(paste0(getwd(),'/scr-12-mapped-filt-mq30.txt.gz'),
                                  paste0(getwd(),'/scr-34-mapped-filt-mq30.txt.gz')),
                           dH1=c(paste0(getwd(),'/dH1-12-mapped-filt-mq30.txt.gz'),
                                 paste0(getwd(),'/dH1-34-mapped-filt-mq30.txt.gz'))),
          filter_file=paste0(getwd(),'/hicdc_significant_indices.txt.gz'),
          output_path=paste0(getwd(),'/diff_analysis_example/'),
          fitType = 'mean',
          binsize=50000,
          diagnostics=TRUE)

