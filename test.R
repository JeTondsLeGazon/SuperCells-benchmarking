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


source('utility.R')
source('supercells.R')
source('analysis.R')


# ---------------------------------------------------------
# Meta parameters
# ---------------------------------------------------------
config <- config::get(file = 'configs/hagai_mouse_lps_config.yml')

filename <- config$filename
data_folder <- file.path("data", config$intermediaryDataFile)
results_folder <- file.path("data", config$resultsFile)
dir.create(results_folder, showWarnings = F, recursive = T)

stat.test <- config$statTest

# Should we compute these DE ?
computeSingle <- config$DE$computeSingle
computeSuper <- config$DE$computeSuper
computeSingleManual <- config$DE$computeSingleManual
computeSuperDes <- config$DE$computeSuperDes
computePseudo <- config$DE$computePseudo
computePseudoManual <- config$DE$computePseudoManual
computeBulk <- config$DE$computeBulk
computeBulkManual <- config$DE$computeBulkManual
computeMeta <- config$DE$computeMeta

gammas <- config$gammas

library(metacell)
library(DropletUtils)
data_folder <- file.path("data", config$intermediaryDataFile)

sc_filtered_data <- readRDS(file = file.path(data_folder, "singleCellFiltered.rds"))
write10xCounts(path = 'data/hagai_mouse_lps_data/sc10x', x = GetAssayData(sc_filtered_data)[, 1:15000], overwrite = T)

scdb_init("data/hagai_mouse_lps_data/sc10x", force_reinit=T)
mcell_import_scmat_10x("meta", base_dir="data/hagai_mouse_lps_data/sc10x")

mcell_add_gene_stat(gstat_id="meta_gs", mat_id="meta", force=T)
mcell_gset_filter_cov(gset_id = "meta_feats", gstat_id="meta_gs", T_tot=100, T_top3=2)
l <- c()
timings <- c()
for(k in c(20, 10))
{
  for(m in c(10, 100))
  {
  start <- Sys.time()
  mcell_add_cgraph_from_mat_bknn(mat_id="meta", 
                                 gset_id = "meta_feats", 
                                 graph_id="meta_graph",
                                 K=k,
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
    K=k, min_mc_size=m, alpha=2)
  
  mat <- scdb_mc("meta_mc")
  l <- c(l, length(unique(mat@mc)))
  timings <- c(timings, difftime(Sys.time(), start, units = "secs"))
  }
}
print(l)
print(timings)
