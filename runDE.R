# 5 January 2022

# Run all differential expression analysis steps and saves results
# Everything is run according to the config file provided

# ---------------------------------------------------------
# Header
# ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
args <- 'configs/hagai_mouse_lps_config.yml'
if (length(args) == 0){
    stop('You must provide a configuration file', call. = FALSE)
}

# SHOULD BE CHANGED ACCORDINGLY TO LOCATIONS OF R LIBRARIES
.libPaths("C:/Users/miche/OneDrive/Documents/R/win-library/4.1")


# ---------------------------------------------------------
# Libraries and dependencies
# ---------------------------------------------------------
library(Seurat)
library(dplyr)
library(Matrix)
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
library(weights)
library(zoo)
library(SuperCellBM)
library(SuperCell)

source('src/utility.R')
source('src/supercells.R')
source('src/analysis.R')
source('src/compute_de.R')

# ---------------------------------------------------------
# Meta parameters
# ---------------------------------------------------------
config <- config::get(file = args[1])

filename <- config$filename
data_folder <- file.path("data", config$intermediaryDataFile)
results_folder <- file.path("data", config$resultsFile)
dir.create(results_folder, showWarnings = F, recursive = T)

# Split either by condition or sample
split.by <- config$splitBy

# Algorithms to run
algos <- config$algos

# Which DE we should compute
computeSingle <- config$DE$computeSingle
computePseudo <- config$DE$computePseudo
computeBulk <- config$DE$computeBulk
computeSuper <- config$DE$computeSuper
computeMeta <- config$DE$computeMeta
computeMetaSC <- config$DE$computeMetaSC
computeRandom <- config$DE$computeRandom
computeSubSampling <- config$DE$computeSubSampling

# Gammas for supercells
gammas <- config$gammas

# Should we compute Supercell or load available files
force_compute <- as.logical(config$compute_supercell)

set.seed(0)


# ---------------------------------------------------------
# Data loadings
# ---------------------------------------------------------
if(!dir.exists(data_folder)){
    stop(sprintf("Cannot load data from folder %s, does not exist", data_folder))
}

single_data <- readRDS(file = file.path(data_folder, "singleCellData.rds"))
pseudobulk_data <- readRDS(file = file.path(data_folder, "pseudoBulkData.rds"))
bulk_data <- readRDS(file = file.path(data_folder, "bulkData.rds"))


# ---------------------------------------------------------
# DE bulk
# ---------------------------------------------------------
if(computeBulk){
    message('Computing Bulk DE')
    for(algo in algos){
        DE <- compute_DE_bulk(data = bulk_data, 
                        labels = bulk_data$label,
                        algo = algo)
        saveMarkers(markers = DE, 
                    algo = algo,
                    split.by = NULL,
                    base.path = results_folder,
                    kind = 'bulk')
    }
    message('Done computing Bulk DE')
}


# ---------------------------------------------------------
# DE pseudobulk
# ---------------------------------------------------------
if(computePseudo){
    message('Computing Pseudo-bulk DE')
    for(algo in algos){
        DE <- compute_DE_bulk(data = pseudobulk_data, 
                              labels = pseudobulk_data$label,
                              algo = algo)
        saveMarkers(markers = DE, 
                    algo = algo,
                    split.by = NULL,
                    base.path = results_folder,
                    kind = 'pseudo')
    }
    message('Done computing Pseudo-bulk DE')
}


# ---------------------------------------------------------
# DE supercells
# ---------------------------------------------------------
# TODO: change this to have each algo ran separately to save them and not lose
# everything if crash happens
if(computeSuper){
    message('Computing SuperCell DE')
    memory.limit(size=56000)
    data <- single_data
    DEs <- list()
    for(gamma in gammas){
        super <- createSuperCellsBM(data = data, 
                                  gamma = gamma,
                                  data_folder = data_folder,
                                  split.by = split.by,
                                  arithmetic = T,
                                  SC.type = 'Exact',
                                  force_compute = force_compute)
        for(algo in algos){
            DE <- compute_supercell_DE(super, algo)
            DEs[[algo]][[as.character(gamma)]] <- DE
        }
    }
    for(algo in algos){
        saveMarkers(markers = DEs[[algo]], 
                    algo = algo,
                    split.by = split.by,
                    base.path = results_folder,
                    kind = 'super')
    }
    message('Done computing SuperCell DE')
}


# ---------------------------------------------------------
# DE single cells
# ---------------------------------------------------------
if(computeSingle){
    message('Computing Single cells DE')
    for(algo in algos){
        DE <- compute_DE_single(single_data, algo)
        saveMarkers(markers = DE, 
                    algo = algo,
                    split.by = split.by,
                    base.path = results_folder,
                    kind = 'single')
    }
    message('Done computing Single cells DE')
}


# ---------------------------------------------------------
# Metacell Footprint
# ---------------------------------------------------------
# TODO: like supercells
if(computeMeta){
    message('Computing MetaCells DE')
    mc.type <- paste('mc', split.by, sep = '_')
    mc <- readRDS(file.path(data_folder, mc.type, 'mc_default.rds'))
    mc_gammas <- names(mc)
    DEs <- list()
    for(mc_gamma in mc_gammas){
        for(algo in algos){
            DE <- compute_DE_meta(mc[[mc_gamma]], algo)
            DEs[[algo]][[mc_gamma]] <- DE
        }
    }
    for(algo in algos){
        saveMarkers(markers = DEs[[algo]], 
                    algo = algo,
                    split.by = split.by,
                    base.path = results_folder,
                    kind = 'meta')
    }
    message('Done computing MetaCells DE')
}


# ---------------------------------------------------------
# Metacell SuperCell-like
# ---------------------------------------------------------
# TODO: like supercells
if(computeMetaSC){
    message('Computing MetaCells SC-like DE')
    mc.type <- paste('mc', split.by, sep = '_')
    mc <- readRDS(file.path(data_folder, mc.type, 'mc_SC_like.rds'))
    mc_gammas <- names(mc)
    DEs <- list()
    for(mc_gamma in mc_gammas){
        for(algo in algos){
            DE <- compute_DE_metasc(data = mc[[mc_gamma]],
                                    single_data = single_data,
                                    algo = algo)
            DEs[[algo]][[mc_gamma]] <- DE
        }
    }
    for(algo in algos){
        saveMarkers(markers = DEs[[algo]], 
                    algo = algo,
                    split.by = split.by,
                    base.path = results_folder,
                    kind = 'metasc')
    }
    message('Done computing MetaCells SC-like DE')
}


# ---------------------------------------------------------
# Random grouping
# ---------------------------------------------------------
# TODO: like supercells
if(computeRandom){
    message('Computing Random grouping DE')
    DEs <- list()
    for(gamma in gammas){
        super <- createSuperCellsBM(data = single_data, 
                                  gamma = gamma,
                                  data_folder = data_folder,
                                  arithmetic = T,
                                  split.by = split.by,
                                  SC.type = 'Random',
                                  force_compute = force_compute)
        for(algo in algos){
            DE <- compute_supercell_DE(super, algo)
            DEs[[algo]][[as.character(gamma)]] <- DE
        }
        
        
    }
    for(algo in algos){
        saveMarkers(markers = DEs[[algo]], 
                    algo = algo,
                    split.by = split.by,
                    base.path = results_folder,
                    kind = 'random')
    }
    message('Done computing Random grouping DE')
}


# ---------------------------------------------------------
# Subsampling
# ---------------------------------------------------------
# TODO: like supercells
if(computeSubSampling){
    message('Computing Subsampling DE')
    DEs <- list()
    data <- single_data
    for(gamma in gammas){
        filename_no_extension <- paste('superCells', gamma, split.by, sep = '_')
        filename <- paste0(filename_no_extension, '.Rds')
        ToComputeSC <- force_compute | !file.exists(data_folder, 'SC', filename)
        SC.list <- compute_supercells(
            sc = single_data,
            ToComputeSC = ToComputeSC,
            data.folder = data_folder,
            filename = filename_no_extension,
            gamma.seq = c(gamma),
            n.var.genes = 1000,
            k.knn = 5,
            n.pc = 10,
            approx.N = 1000,
            fast.pca = TRUE,
            genes.use = NULL, 
            genes.exclude = NULL,
            seed.seq = c(0)
        )
        cells.idx <- SC.list$Subsampling[[as.character(gamma)]][[1]]$cells.use.idx
        
        data$keep <- rep(F, ncol(single_data))
        data$keep[cells.idx] <- T
        keep.data <- subset(data, subset = keep == T)
        for(algo in algos){
            DE <- compute_DE_single(keep.data, algo)
            DEs[[algo]][[as.character(gamma)]] <- DE
        }
    }
    for(algo in algos){
        saveMarkers(markers = DEs[[algo]], 
                    algo = algo,
                    split.by = split.by,
                    base.path = results_folder,
                    kind = 'subsampling')
    }
    message('Done computing Subsampling DE')
}
