library(stringr)
library(pheatmap)
library(ggplot2)

#### FUNCTIONS ####

## Takes output of MEME AME and make it into a neat dataframe
process_motif_data <- function(base_dir, expressed_genes = NULL) {
  library(dplyr)
  library(readr)
  
  # Initialize an empty list to store the data frames
  Peaks_AME <- list()
  processed_dirs <- vector()
  setwd(base_dir)
  file_paths=dir()
  
  # Loop through each directory and process the files
  for (i in seq_along(file_paths)) {
    # Construct the path to the ame.tsv file
    ame_file_path <- paste0(file_paths[i], "/ame.tsv")
    
    # Check the file size
    file_info <- file.info(ame_file_path)
    if (file_info$size == 1.1 || is.na(file_info$size)) {
      # If the file size is 0 or file does not exist, skip to the next iteration
      next
    }
    
    # Read the ame.tsv file
    Peaks_AME[[length(Peaks_AME) + 1]] <- as.data.frame(read_tsv(file = ame_file_path, comment = "#"))
    
    # Add the current directory to the processed_dirs vector
    processed_dirs <- c(processed_dirs, file_paths[i])
  }
  
  # Set names for the list elements
  names(Peaks_AME) <- processed_dirs
  
  # Post-processing each dataframe
  for (i in seq_along(Peaks_AME)) {
    
    # Calculate TPvsFP and add as a new column
    Peaks_AME[[i]]$TPvsFP <- Peaks_AME[[i]]$`%TP` - Peaks_AME[[i]]$`%FP`
    
    # Remove everything after "_" in the motif_ID column
    Peaks_AME[[i]]$motif_ID <- gsub("\\_.*", "", Peaks_AME[[i]]$motif_ID)
  }
  
  # Convert the list of MEME-AME results into a single dataframe
  Peaks_AME <- Peaks_AME[!sapply(Peaks_AME, is.null)]
  Peaks_AME_df <- do.call(rbind, Peaks_AME)
  
  # Filter for motifs that are expressed, if provided
  if (!is.null(expressed_genes)) {
    Peaks_AME_df <- Peaks_AME_df %>%
      filter(motif_ID %in% expressed_genes)
  }
  colnames(Peaks_AME_df)[c(7,15,17)] <- c("adj_pvalue", "percentTP", "percentFP") # change column names without special symbols 
  Peaks_AME_df$neglog10padj <- -log10(Peaks_AME_df$adj_pvalue) # add in column for -log10(adj pvalue)
  return(Peaks_AME_df)
}

## Makes matrix of most significant peaks per group from df of all motif data
create_AME_peaks_matrix_byGroup <- function(sig_peaks_df) {
  # Create an empty matrix for significant peaks
  unique_motifs <- unique(sig_peaks_df$motif_ID)
  group <- unique(sig_peaks_df$updown)
  sig_peaks_mat <- matrix(nrow = length(unique_motifs), ncol = length(group))
  
  rownames(sig_peaks_mat) <- unique_motifs # Row names are motif IDs
  colnames(sig_peaks_mat) <- group   # Column names are DESeq comparisons
  
  # Populate the matrix with neglog10padj values
  for (i in 1:ncol(sig_peaks_mat)) {
    temp_comparison <- sig_peaks_df[sig_peaks_df$updown == group[i], ]
    sig_peaks_mat[, i] <- temp_comparison[match(rownames(sig_peaks_mat), temp_comparison$motif_ID), ]$neglog10padj
    sig_peaks_mat[, i] <- ifelse(is.na(sig_peaks_mat[, i]), 0, sig_peaks_mat[, i])
  }
  
  return(sig_peaks_mat)
}

## Makes dotplots of motif enrichment results
create_AME_dotplots <- function(AME_df=AME_padj_dots, Title=Title){
  ggplot(AME_df, aes(x=updown, y=motif_ID, size=TPvsFP)) +
    geom_point(aes(fill=neglog10padj), colour="black", pch=21, alpha=0.8) +
    scale_fill_gradient(name="-log10(padj)", low="light yellow", high="dark orange") +
    scale_size_continuous(name="%TP - % FP") +
    theme_minimal() +
    ggtitle(Title) +
    theme(
      text = element_text(size=8), 
      axis.text.x = element_text(size=8), 
      strip.text.x = element_text(size=8, face="bold"),
      strip.text.y = element_text(size=8, face="bold")
    ) +
    theme_classic()+
    facet_wrap(~updown, scales="free", nrow=1)+
    labs(x="peak set", y="motif ID")
}


#### MAIN CODE ####

## PROCESSING MEME OUTPUT
Peaks_AME_updown <- process_motif_data("/lustre/fs4/risc_lab/scratch/iduba/linker-histone/multi-results/ATACxRNA/newRNA/MEME-AME-out")
Peaks_AME_updown$updown <- gsub('\\.[0-9]*$', '', row.names(Peaks_AME_updown)) 

## PLOTTING
AME_padj_mat <- list()
comparisons <- c("down-peaks-down-genes", "down-peaks-up-genes", "up-peaks-down-genes", "up-peaks-up-genes")
for(i in 1:length(comparisons)){
  AME_padj_mat[[i]] <- create_AME_peaks_matrix_byGroup(Peaks_AME_updown[Peaks_AME_updown$updown %in% comparisons[i], ])
}
names(AME_padj_mat) <- c("down-peaks-down-genes", "down-peaks-up-genes", "up-peaks-down-genes", "up-peaks-up-genes")

pdf("/lustre/fs4/risc_lab/scratch/iduba/linker-histone/multi-results/ATACxRNA/newRNA/MEME-AME-out/Up.Down.peaks.near.up.down.genes.enriched.motifs.pdf", width=3, height=9)
  pheatmap(AME_padj_mat[[1]][order(rowSums(AME_padj_mat[[1]]), decreasing=T), ], cluster_rows = F, cluster_cols=F, color = colorRampPalette(c("white",  "#FDBB84", "#B30000"))(100), breaks = seq(0, 150,length.out = 100), main="down peaks near down genes")
  pheatmap(AME_padj_mat[[2]][order(rowSums(AME_padj_mat[[2]]), decreasing=T), ], cluster_rows = F, cluster_cols=F, color = colorRampPalette(c("white",  "#FDBB84", "#B30000"))(100), breaks = seq(0, 30,length.out = 100), main="down peaks near up genes")
  pheatmap(AME_padj_mat[[3]][order(rowSums(AME_padj_mat[[3]]), decreasing=T), ], cluster_rows = F, cluster_cols=F, color = colorRampPalette(c("white",  "#FDBB84", "#B30000"))(100), breaks = seq(0, 40,length.out = 100), main="up peaks near down genes") 
  pheatmap(AME_padj_mat[[4]][order(rowSums(AME_padj_mat[[4]]), decreasing=T), ], cluster_rows = F, cluster_cols=F, color = colorRampPalette(c("white",  "#FDBB84", "#B30000"))(100), breaks = seq(0, 20,length.out = 100), main="up peaks near up genes")
dev.off()

AME_padj_dots <- list()
for(i in 1:length(comparisons)){
  AME_padj_dots[[i]] <- Peaks_AME_updown[Peaks_AME_updown$updown %in% comparisons[i], ]
  AME_padj_dots[[i]] <- AME_padj_dots[[i]][AME_padj_dots[[i]]$motif_ID %in% row.names(AME_padj_mat[[i]])[1:15], ] # show only top 15 most significant padj values motif IDs
  AME_padj_dots[[i]]$motif_ID <- factor(AME_padj_dots[[i]]$motif_ID, levels=rev(unique(AME_padj_dots[[i]]$motif_ID)))
}
names(AME_padj_dots) <- c("down-peaks-down-genes", "down-peaks-up-genes", "up-peaks-down-genes", "up-peaks-up-genes")

pdf("/lustre/fs4/risc_lab/scratch/iduba/linker-histone/multi-results/ATACxRNA/newRNA/MEME-AME-out/Up.Down.peaks.near.RNAclusters.enriched.motifs.dotplots.pdf", width=8, height=5)
for(i in 1:length(AME_padj_dots)){
  p <- create_AME_dotplots(AME_padj_dots[[i]], Title=names(AME_padj_dots)[i]) 
  print(p)
}
dev.off()
