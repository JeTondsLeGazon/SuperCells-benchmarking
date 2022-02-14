# 11 January 2022

# Script for the computation of metacells, as the memory requirement is high
# Can only be run on linux and macOS

# ---------------------------------------------------------
#  Header
# ---------------------------------------------------------
if(.Platform$OS.type == 'windows'){
    stop('Cannot run this script on Windows... Try on Linux or Mac!')
}

args <- commandArgs(trailingOnly = TRUE)
args <- 'configs/hagai_mouse_lps_config.yml'
if (length(args) == 0){
    stop('You must provide a configuration file', call. = FALSE)
}

# SHOULD BE CHANGED ACCORDINGLY TO LOCATIONS OF R LIBRARIES
#.libPaths("C:/Users/miche/OneDrive/Documents/R/win-library/4.1")


# ---------------------------------------------------------
# Libraries and dependencies
# ---------------------------------------------------------
library(Seurat)
library(ggplot2)
library(data.table)
library(DESeq2)
library(edgeR)
library(tidyr)
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
library(SuperCellBM)
library(SuperCell)

source('src/utility.R')
source('src/supercells.R')
source('src/analysis.R')


# ---------------------------------------------------------
# Functions
# ---------------------------------------------------------
# As metacells are not saved within each group (sample or condition), we need
# to find common genes for the same gamma
find_common_genes <- function(mc, mc.type){
    genes.per.gamma <- list()
    for(group in mc){
        data <- group[[mc.type]]
        gammas <- as.numeric(names(data))
        if(length(genes.per.gamma) == 0){
            genes.per.gamma <- sapply(seq_along(gammas), function(x) list())
        }
        for(i in seq_along(gammas)){
            ge <- data[[i]][[1]]$mc_info$mc@mc_fp
            
            
            #find intersection of all genes
            if(length(genes.per.gamma[[i]]) == 0){
                genes.per.gamma[[i]] <- rownames(ge)
            }else{
                genes.per.gamma[[i]] <- intersect(genes.per.gamma[[i]],
                                                  rownames(ge))
            }
        }
    }
    return(genes.per.gamma)
}

# reorder metacells in a better data structure for the Differential expression
create_metacell <- function(genes, groups, mc, mc.type, single_data){
    MC <- sapply(seq_len(20), function(x) list())
    MC.mcClass <- sapply(seq_len(20), function(x) list())
    
    # Rearrange data
    for(group in groups){
        data <- mc[[group]][[mc.type]]
        gammas <- as.numeric(names(data))
        for(i in seq_along(gammas)){
            keep.genes <- genes[[i]]
            ge <- data[[i]][[1]]$mc_info$mc@mc_fp[keep.genes, ]
            e_gc <- data[[i]][[1]]$mc_info$mc@e_gc[keep.genes, ]
            if(length(MC[[i]]) ==  0){
                MC[[i]]$ge <- ge
                MC[[i]]$sample <- rep(group, ncol(ge))
                MC[[i]]$size <- table(data[[i]][[1]]$mc_info$mc@mc)
                MC[[i]]$membership <- data[[i]][[1]]$mc_info$mc@mc
                MC.mcClass[[i]] <- data[[i]][[1]]$mc_info$mc
            }else{
                MC[[i]]$ge <- cbind(MC[[i]]$ge,
                                    ge)
                MC[[i]]$sample <- c(MC[[i]]$sample,
                                    rep(group, ncol(ge)))
                MC[[i]]$size <- c(MC[[i]]$size, 
                                  table(data[[i]][[1]]$mc_info$mc@mc))
                names(MC[[i]]$size) <- seq_along(MC[[i]]$size)
                MC[[i]]$membership <- c(MC[[i]]$membership,
                                        data[[i]][[1]]$mc_info$mc@mc + max(MC[[i]]$membership))
                
                MC.mcClass[[i]]@mc <- c(MC.mcClass[[i]]@mc,
                                        
                                        data[[i]][[1]]$mc_info$mc@mc + max(MC.mcClass[[i]]@mc))
            }
        }
    }
    MC <- MC[sapply(MC, function(x) length(x) > 0)]
    N.sc <- ncol(single_data)
    # Rename with gammas
    for(i in seq_along(MC)){
        ge <- MC[[i]]$ge
        new.colnames <- paste(MC[[i]]$sample, seq_len(ncol(ge)), sep = '_')
        colnames(MC[[i]]$ge) <- new.colnames
        MC[[i]]$gamma <- round(N.sc / ncol(ge))
        
        used.cells <- names(MC.mcClass[[i]]@mc)
        used.cells.idx <- which(used.cells %in% colnames(single_data@assays$RNA@counts))
        used.cells <- used.cells[used.cells.idx]
        MC.mcClass[[i]]@mc <- MC.mcClass[[i]]@mc[used.cells]
        MC.mcClass[[i]] <- mc_compute_fp(mc = MC.mcClass[[i]],
                                         us = single_data@assays$RNA@counts[, used.cells])
        
        MC[[i]]$mc_fp_updated <- MC.mcClass[[i]] 
    }
    
    names(MC) <- sapply(MC, function(x) x$gamma)
    return(MC)
}


# ---------------------------------------------------------
# Meta parameters
# ---------------------------------------------------------
config <- config::get(file = args[1])

filename <- config$filename
data_folder <- file.path("data", config$intermediaryDataFile)
results_folder <- file.path("data", config$resultsFile)
dir.create(results_folder, showWarnings = F, recursive = T)

split.by <- config$splitBy
data_folder <- file.path("data", config$intermediaryDataFile)


# ---------------------------------------------------------
# Data loadings
# ---------------------------------------------------------
sc_data <- readRDS(file = file.path(data_folder, "singleCellData.rds"))
groups <- unique(sc_data[[split.by]][[1]])

set.seed(0)
# ---------------------------------------------------------
# Metacell creation
# ---------------------------------------------------------
SC.mc <- list()
for(group in groups){
  cells.id <- which(sc_data[[split.by]] == as.character(group))
  filename <- paste('superCells', split.by, sep = '_')
  SC.list <- compute_supercells(
    sc = sc_data[, cells.id],
    ToComputeSC = T,
    data.folder = data_folder,
    filename = filename,
    gamma.seq = c(1, 2, 10, 20, 50, 100),
    n.var.genes = 1000,
    k.knn = 5,
    n.pc = 10,
    approx.N = 1000,
    fast.pca = TRUE,
    genes.use = NULL, 
    genes.exclude = NULL,
    seed.seq = 0,
    split.by = split.by
  )
  
  
  res <- compute_supercells_metacells_with_min_mc_size(
      sc.counts = sc_data@assays$RNA@counts[, cells.id],
      SC.list = SC.list,
      min_mc_size_seq =  c(1, 10, 20, 30, 50),
      proj.name = 'metacell',
      ToComputeSC = T, 
      mc.k.knn = 100,
      T_vm_def = 0.08,
      MC.folder = file.path(data_folder, 'MC'), 
      MC_gene_settings = c('Metacell_default', 'Metacell_SC_like')
  )
  SC.mc[[as.character(group)]] <- res
}


# ---------------------------------------------------------
# Metacell cleaning and reshaping
# ---------------------------------------------------------
for(mc.type in c('Metacell_default', 'Metacell_SC-like')){
    sample_vs_condition <- paste('mc', split.by, sep = '_')
    mc_folder <- file.path(results_folder, sample_vs_condition)
    N.sc <- ncol(sc_data)
    
    # Default
    genes <- find_common_genes(mc = SC.mc,
                               mc.type = mc.type)
    
    MC <- create_metacell(genes = genes,
                          groups = groups,
                          mc = SC.mc,
                          mc.type = mc.type,
                          single_data = sc_data)
    if(mc.type == 'Metacell_default'){
        filename <- 'mc_default.rds'
    }else{
        filename <- 'mc_SC_like.rds'
    }
    if(!dir.exists(mc_folder)){
      dir.create(mc_folder, recursive = T)
    }
    saveRDS(MC, file.path(mc_folder, filename))
}
