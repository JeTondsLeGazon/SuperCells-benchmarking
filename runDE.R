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
library(limma)
library(ggpubr)
library(ggrepel)
library(stringr)
library(tidyseurat)
library(ggExtra)
library(weights)
library(zoo)
library(SuperCellBM)

source('src/utility.R')
source('src/supercells.R')
source('src/analysis.R')
source('src/DEbulk.R')

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

set.seed(0)


# ---------------------------------------------------------
# Data loadings
# ---------------------------------------------------------
if(!dir.exists(data_folder)){
    stop(sprintf("Cannot load data from folder %s, does not exist", data_folder))
}

single_data <- readRDS(file = file.path(data_folder, "singleCellClusteredNormalized.rds"))
pseudobulk_data <- readRDS(file = file.path(data_folder, "pseudoBulkNormalized.rds"))
bulk_data <- readRDS(file = file.path(data_folder, "bulkFilteredNormalized.rds"))


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
        super <- superCellWrapper(data = data, 
                                  gamma = gamma, 
                                  split.by = split.by,
                                  arithmetic = T,
                                  SC.type = 'Exact')
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
# TODO: include seurat computation as algorithm
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
# Metacell via SuperCell t-test
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
# Random grouping t-test
# ---------------------------------------------------------
if(computeRandomGroup){
    message('Computing Random grouping t-test')
    data = sc_clustered_data
    DEs <- list()
    for(gamma in gammas){
        super <- superCellWrapper(data = data, 
                                  gamma = gamma, 
                                  arithmetic = T,
                                  split.by = split.by,
                                  SC.type = 'Random',
                                  bm = T)

        rdm_markers <- supercell_FindMarkers(ge = super$GE,
                                       supercell_size = super$supercell_size,
                                       clusters = super$cell_line,
                                       ident.1 = 'treat',
                                       ident.2 = 'ctrl',
                                       logfc.threshold = 0,
                                       only.pos = F,
                                       do.bootstrapping = F,
                                       test.use = stat.test)
        
        DEs[[as.character(gamma)]] <- arrangeDE(rdm_markers)
    }
    saveMarkers(markers = DEs, 
                algo = 't-test',
                split.by = split.by,
                base.path = results_folder,
                kind = 'random')
    message('Done computing Random grouping t-test')
}


# ---------------------------------------------------------
# Random grouping DESeq2
# ---------------------------------------------------------
if(computeRandomGroupDes){
    memory.limit(size=56000)
    message('Computing Random grouping DESeq2')
    DEs <- list()
    data = sc_filtered_data
    for(gamma in gammas){
        super <- superCellWrapper(data = data, 
                                  gamma = gamma, 
                                  arithmetic = F,
                                  split.by = split.by,
                                  SC.type = 'Random',
                                  bm = T,
                                  norm = F)
        ge <- floor(sweep(super$GE, 2, super$supercell_size, '*'))
        dds <- DESeqDataSetFromMatrix(ge,
                                      colData = data.frame(design = super$cell_line), 
                                      design = ~ design)
        
        dds_wald <- DESeq(dds, test = 'Wald', minReplicatesForReplace = Inf)
        results_wald <- results(dds_wald)
        
        DEs[[as.character(gamma)]] <- arrangeDE(results_wald, 
                                                   oldNameLog = 'log2FoldChange',
                                                   oldNameP = 'padj')
    }
    saveMarkers(markers = DEs, 
                algo = 'DESeq2',
                split.by = split.by,
                base.path = results_folder,
                kind = 'random')
    message('Done computing Random grouping DESeq2')
}



# ---------------------------------------------------------
# Random grouping EdgeR
# ---------------------------------------------------------
if(computeRandomGroupEdge){
    memory.limit(size=56000)
    message('Computing Random grouping EdgeR')
    DEs <- list()
    data = sc_filtered_data
    for(gamma in gammas){
        super <- superCellWrapper(data = data, 
                                  gamma = gamma, 
                                  arithmetic = F,
                                  split.by = split.by,
                                  SC.type = 'Random',
                                  bm = T,
                                  norm = F)
        ge <- floor(sweep(super$GE, 2, super$supercell_size, '*'))

        edge <- DGEList(counts = ge, group = super$cell_line)
        edge <- calcNormFactors(edge)
        model <- model.matrix(~super$cell_line)
        edge <- estimateDisp(edge, model)
        
        # Quasi likelihood test
        fit <- glmQLFit(edge, model)
        DEs[[as.character(gamma)]] <- arrangeDE(glmQLFTest(fit, coef= 2)$table,
                                                   oldNameP = 'PValue')
    }
    saveMarkers(markers = DEs, 
                algo = 'EdgeR',
                split.by = split.by,
                base.path = results_folder,
                kind = 'random')
    message('Done computing Random grouping EdgeR')
}


# ---------------------------------------------------------
# Subsampling
# ---------------------------------------------------------
if(computeSubSampling){
    message('Computing Subsampling t-test')
    data = sc_filtered_data
    DEs <- list()
    for(gamma in gammas){
        SC.list <- compute_supercells(
            sc = data,
            ToComputeSC = T,
            data.folder = 'data/',
            filename = paste0('superCells', gamma),
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
        data$keep <- rep(F, ncol(data))
        data$keep[cells.idx] <- T
        keep.data <- subset(data, subset = keep == T)
        sub_markers <- find_markers(keep.data, stat.test)
        DEs[[as.character(gamma)]] <- sub_markers
    }
    saveMarkers(markers = DEs, 
                algo = 't-test',
                split.by = split.by,
                base.path = results_folder,
                kind = 'subsampling')
    message('Done computing Subsampling t-test')
}


# ---------------------------------------------------------
# Subsampling DESeq2
# ---------------------------------------------------------
if(computeSubSamplingDes){
    message('Computing Subsampling DESeq2')
    data = sc_filtered_data
    DEs <- list()
    for(gamma in gammas){
        SC.list <- compute_supercells(
            sc = data,
            ToComputeSC = T,
            data.folder = 'data/',
            filename = paste0('superCells', gamma),
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
        ge <- GetAssayData(sc_filtered_data)[, cells.idx]
        labels <- data.frame(label = data$label[cells.idx])
        dds <- DESeqDataSetFromMatrix(ge,
                                      colData = labels, 
                                      design = ~ label)
        
        dds_wald <- DESeq(dds, test = 'Wald', minReplicatesForReplace = Inf)
        results_wald <- results(dds_wald)
        
        DEs[[as.character(gamma)]] <- arrangeDE(results_wald, 
                                                   oldNameLog = 'log2FoldChange',
                                                   oldNameP = 'padj')
    }
    saveMarkers(markers = DEs, 
                algo = 'DESeq2',
                split.by = split.by,
                base.path = results_folder,
                kind = 'subsampling')
    message('Done computing Subsampling DESeq2')
}

# ---------------------------------------------------------
# Subsampling EdgeR
# ---------------------------------------------------------
if(computeSubSamplingEdge){
    message('Computing Subsampling EdgeR')
    data = sc_filtered_data
    DEs <- list()
    for(gamma in gammas){
        SC.list <- compute_supercells(
            sc = data,
            ToComputeSC = T,
            data.folder = 'data/',
            filename = paste0('superCells', gamma),
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
        ge <- GetAssayData(sc_filtered_data)[, cells.idx]
        labels <- data$label[cells.idx]
        
        # EdgeR DE
        edge <- DGEList(counts = ge, group = labels)
        edge <- calcNormFactors(edge)
        model <- model.matrix(~labels)
        edge <- estimateDisp(edge, model)
        
        # Quasi likelihood test
        fit <- glmQLFit(edge, model)
        DEs[[as.character(gamma)]] <- arrangeDE(glmQLFTest(fit, coef= 2)$table,
                                                oldNameP = 'PValue')
    }
    saveMarkers(markers = DEs, 
                algo = 'EdgeR',
                split.by = split.by,
                base.path = results_folder,
                kind = 'subsampling')
    message('Done computing Subsampling EdgeR')
}
