# 11 January 2022

# Script for the differential expression computation of metacells, as the memory
# required is high

# ---------------------------------------------------------
# Libraries and dependencies
# ---------------------------------------------------------

library(Seurat)
library(ggplot2)
library(data.table)
library(DESeq2)
library(edgeR)
library(tidyr)
library(limma)
library(ggpubr)
library(ggrepel)
library(stringr)
library(tidyseurat)
library(ggExtra)
library(dplyr)
library(Matrix)
library(ggplot2)
library(weights)
library(zoo)
library(metacell)
library(DropletUtils)

source('src/utility.R')
source('src/supercells.R')
source('src/analysis.R')

# ---------------------------------------------------------
# Meta parameters
# ---------------------------------------------------------
config <- config::get(file = 'configs/hagai_mouse_lps_config.yml')

filename <- config$filename
data_folder <- file.path("data", config$intermediaryDataFile)
results_folder <- file.path("data", config$resultsFile)
dir.create(results_folder, showWarnings = F, recursive = T)


# ---------------------------------------------------------
# Data loadings
# ---------------------------------------------------------
sc_filtered_data <- readRDS(file = file.path(data_folder, "singleCellFiltered.rds"))
samples <- unique(sc_filtered_data$sample)

# data should be in 10x format for metacell to work
for(sample in samples){
  path_per_sample <- file.path(data_folder, paste0('sc10x', sample))
  write10xCounts(path = path_per_sample, 
                 x = GetAssayData(sc_filtered_data)[, which(sc_filtered_data$sample == sample)], 
                 overwrite = T)
}


# ---------------------------------------------------------
# Metacells creation
# ---------------------------------------------------------
samples_mc <- list()
for(sample in samples){
  path_per_sample <- file.path(data_folder, paste0('sc10x', sample))
  scdb_init(path_per_sample, force_reinit=T)
  mcell_import_scmat_10x("meta", base_dir = path_per_sample)
  
  mcell_add_gene_stat(gstat_id="meta_gs", mat_id="meta", force=T)
  mcell_gset_filter_cov(gset_id = "meta_feats", gstat_id="meta_gs", T_tot=100, T_top3=10)
  
  sizes <- c(1, 2, 5, 10, 50, 100)
  sample_mc <- list()
  for(size in sizes){
    
    mcell_add_cgraph_from_mat_bknn(mat_id="meta", 
                                   gset_id = "meta_feats", 
                                   graph_id="meta_graph",
                                   K=100,
                                   dsamp=T)
    mcell_coclust_from_graph_resamp(
      coc_id="meta_coc500", 
      graph_id="meta_graph",
      min_mc_size=20, 
      p_resamp=0.75, n_resamp=500)
    
    mcell_mc_from_coclust_balanced(
      coc_id="meta_coc500", 
      mat_id= "meta",
      mc_id= "meta_mc", 
      K=30, min_mc_size=size, alpha=2)
    
    mat <- scdb_mc("meta_mc")
    sample_mc[[size]] <- mat@mc
    print(length(mat@mc))
  }
  samples_mc[[sample]] <- sample_mc
}
saveRDS(samples_mc, file.path(results_folder, "mcComposition.rds"))
