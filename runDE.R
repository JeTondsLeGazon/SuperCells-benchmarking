# 5 January 2022

# Run all differential expression analysis steps and saves results
# Everything is run according to the config file provided

# ---------------------------------------------------------
# Header
# ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

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

source('src/utility.R')
source('src/supercells.R')
source('src/analysis.R')


# ---------------------------------------------------------
# Meta parameters
# ---------------------------------------------------------
config <- config::get(file = args[1])

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
computeSuperEdge <- config$DE$computeSuperEdge
computeMetaDes <- config$DE$computeMetaDes
computeMetaEdge <- config$DE$computeMetaEdge
computeMetaSuper <- config$DE$computeMetaSuper
computeMetaSuperDes <- config$DE$computeMetaSuperDes
computeMetaSuperEdge <- config$DE$computeMetaSuperEdge
computeRandomGroup <- config$DE$computeRandomGroup

gammas <- config$gammas

set.seed(0)


# ---------------------------------------------------------
# Data loadings
# ---------------------------------------------------------
if(!dir.exists(data_folder)){
    stop(sprintf("Cannot load data from folder %s, does not exist", data_folder))
}

sc_clustered_data <- readRDS(file = file.path(data_folder, "singleCellClusteredNormalized.rds"))
sc_filtered_data <- readRDS(file = file.path(data_folder, "singleCellFiltered.rds"))
pseudobulk_data <- readRDS(file = file.path(data_folder, "pseudoBulk.rds"))
pseudobulk_norm <- readRDS(file = file.path(data_folder, "pseudoBulkNormalized.rds"))
bulk_filtered_data <- readRDS(file = file.path(data_folder, "bulkFiltered.rds"))
bdata <- readRDS(file = file.path(data_folder, "bulkFilteredNormalized.rds"))


# ---------------------------------------------------------
# DE bulk
# ---------------------------------------------------------
if(computeBulk){
    message('Computing Bulk DE genes')
    bulk_markers <- compute_DE_bulk(bulk_filtered_data)
    volcano_plot(bulk_markers$`edgeR-QLF`) +
        ggtitle('Volcano plot of bulk data from DESeq2 (wald)') +
        theme(plot.title = element_text(hjust = 0.5))
    
    bulk_markers <- lapply(bulk_markers, function(x) x %>% subset(logFC > 0))
    bulk_markers <- lapply(bulk_markers, function(x) if(nrow(x) == 0){x <- NULL}else{x})
    bulk_markers <- bulk_markers[!unlist(lapply(bulk_markers, is.null))]
    saveRDS(bulk_markers, file.path(results_folder, "bulkMarkers.rds"))
    message('Done computing Bulk DE genes')
}



# ---------------------------------------------------------
# DE bulk manual
# ---------------------------------------------------------
if(computeBulkManual){
    message('Computing Bulk DE genes manually')
    manual_bulk_markers <- arrangeDE(find_markers_bulk(bulk_filtered_data, stat.test))
    saveRDS(manual_bulk_markers, file.path(results_folder, "bulkMarkersManual.rds"))
    message('Done computing Bulk DE genes manually')
    
}


# ---------------------------------------------------------
# DE pseudobulk DESeq2
# ---------------------------------------------------------
if(computePseudo){
    message('Computing Pseudo-bulk DE genes')
    pseudo_markers <- compute_DE_bulk(pseudobulk_data)
    volcano_plot(pseudo_markers$`edgeR-QLF`) +
        ggtitle('Volcano plot of pseudo bulk data from DESeq2 (wald)') +
        theme(plot.title = element_text(hjust = 0.5))
    
    pseudo_markers <- lapply(pseudo_markers, function(x) subset(x, logFC > 0))
    pseudo_markers <- lapply(pseudo_markers, function(x) if(nrow(x) == 0){x <- NULL}else{x})
    pseudo_markers <- pseudo_markers[!unlist(lapply(pseudo_markers, is.null))]
    saveRDS(pseudo_markers, file.path(results_folder, "pseudoMarkers.rds"))
    message('Done computing Pseudo-bulk DE genes')
}


# ---------------------------------------------------------
# DE pseudobulk manual
# ---------------------------------------------------------
if(computePseudoManual){
    message('Computing Pseudo-bulk DE genes manually')
    pseudo_markers_manual <- find_markers_bulk(pseudobulk_norm, stat.test) %>%
        arrange(adj.p.value, 1 / (abs(logFC) + 1)) %>%
        mutate(gene = row.names(.)) %>%
        subset(logFC > 0)
    saveRDS(pseudo_markers_manual, file.path(results_folder, "pseudoMarkersManual.rds"))
    message('Done computing Pseudo-bulk DE genes manually')
}


# ---------------------------------------------------------
# DE supercells (weighted and unweighted)
# ---------------------------------------------------------
if(computeSuper){
    message('Computing SuperCell DE genes')
    memory.limit(size=56000)
    super_markers <- superCells_DEs(sc_clustered_data, gammas, 5,
                                    weighted = F,
                                    test.use = stat.test)
    
    volcano_plot(super_markers$`1`, logfc.thres = 0.5) +
        ggtitle('Volcano plot of SuperCells at level gamma = 5') +
        theme(plot.title = element_text(hjust = 0.5))
    
    super_markers <- lapply(super_markers, function(x) subset(x, logFC > 0))
    saveRDS(super_markers, file.path(results_folder, "superMarkers.rds"))
    message('Done computing SuperCell DE genes')
    
    message('Computing SuperCell DE genes with weighted t-test')
    super_markers_weighted <- superCells_DEs(sc_clustered_data, gammas, 5,
                                    weighted = T,
                                    test.use = stat.test)
    
    super_markers_weighted <- lapply(super_markers_weighted, 
                                     function(x) subset(x, logFC > 0))
    saveRDS(super_markers_weighted, file.path(results_folder, "superMarkersWeighted.rds"))
    message('Done computing SuperCell DE genes with weighted t-test')
}


# ---------------------------------------------------------
# DE supercells DESeq2
# ---------------------------------------------------------
if(computeSuperDes){
    message('Computing SuperCell DE genes with DESeq2')
    super_markers_des <- list()
    
    for(gamma in gammas){
        message('\tGamma = ', gamma)
        # supercells creation with geometric mean as using counts here
        super <- createSuperCells(sc_filtered_data, gamma, arithmetic = FALSE)
        ge <- floor(sweep(super$GE, 2, super$supercell_size, '*'))
        dds <- DESeqDataSetFromMatrix(ge,
                                      colData = super$sc.cell.annotation., 
                                      design = ~ sc.cell.annotation.)
        
        dds_wald <- DESeq(dds, test = 'Wald', minReplicatesForReplace = Inf)
        results_wald <- results(dds_wald)
        
        super_markers_des[[as.character(gamma)]] <- arrangeDE(results_wald, 
                                                              oldNameLog = 'log2FoldChange',
                                                              oldNameP = 'padj')
    }
    saveRDS(super_markers_des, file.path(results_folder, "superMarkersDes.rds"))
    message('Done computing SuperCell DE genes with DESeq2')
}


# ---------------------------------------------------------
# DE supercells EdgeR
# ---------------------------------------------------------
# Saves at each gamma to avoid loss of information in case of crash
if(computeSuperEdge){
    message('Computing SuperCell DE genes with EdgeR')
    for(gamma in gammas){
        message('\tGamma = ', gamma)
        # supercells creation with geometric mean as using counts here
        super <- createSuperCells(sc_filtered_data, gamma, arithmetic = FALSE)
        
        ge <- floor(sweep(super$GE, 2, super$supercell_size, '*'))
        edge <- DGEList(counts = ge, group = super$sc.cell.annotation.)
        edge <- calcNormFactors(edge)
        model <- model.matrix(~design)
        edge <- estimateDisp(edge, model)
        
        # Quasi likelihood test
        fit <- glmQLFit(edge, model)
        qlf <- arrangeDE(glmQLFTest(fit, coef= 2)$table, oldNameP = 'PValue')
        
        saveRDS(qlf, file.path(results_folder, sprintf("superMarkersEdge%s.rds", gamma)))
        
        # Memory space savings to avoid crash ...
        rm('super')
        rm('qlf')
    }
    files <- list.files(results_folder)
    files <- files[grep('superMarkersEdge.*.rds', files)]
    super_markers_edge <- list()
    for(f in files){
        gamma <- str_split(str_split(f, 'Edge', simplify = T)[2], '.rds', simplify = T)[1]
        if(gamma %in% gammas){
            d <- readRDS(file.path(results_folder, f))
            super_markers_edge[[gamma]] <- d
        }
    }
    # reordering
    idx <- sort(as.numeric(names(super_markers_edge)), index.return = T)$ix
    saveRDS(super_markers_edge[idx], file.path(results_folder, "superMarkersEdge.rds"))
    message('Done computing SuperCell DE genes with EdgeR')
}


# ---------------------------------------------------------
# DE single cells (seurat)
# ---------------------------------------------------------
if(computeSingle){
    message('Computing Single cells DE genes with Seurat')
    single_markers <- singleCell_DE(sc_clustered_data, 
                                    var.features = 500,
                                    stat.test)
    volcano_plot(single_markers, logfc.thres = 0.5) +
        ggtitle('Volcano plot of single cells from FindAllMarkers (seurat)') +
        theme(plot.title = element_text(hjust = 0.5))
    
    single_markers <- single_markers %>% subset(logFC > 0)
    saveRDS(single_markers, file.path(results_folder, "singleMarkers.rds"))
    message('Done computing Single cells DE genes with Seurat')
}


# ---------------------------------------------------------
# DE single cells (manual)
# ---------------------------------------------------------
if(computeSingleManual){
    message('Computing Single cells DE genes manually')
    manual_single_markers <- find_markers_bulk(sc_clustered_data, stat.test)
    manual_single_markers <- manual_single_markers %>% 
        arrange(adj.p.value, 1 / (abs(logFC) + 1), T) %>%
        mutate(gene = rownames(.)) %>%
        subset(logFC > 0)
    saveRDS(manual_single_markers, file.path(results_folder, "singleMarkersManual.rds"))
    message('Done computing Single cells DE genes manually')
}


# ---------------------------------------------------------
# Metacell GE t-test
# ---------------------------------------------------------
# should run runMeta.R first
if(computeMeta){
    # Metacell own GE matrix from pipeline, manual t-test
    message('Computing MetaCells DE genes with t-test')
    GEs <- createMCGE(sc_filtered_data, results_folder)
    mc_gammas <- names(GEs)
    DEs <- list()
    for(mc_gamma in mc_gammas){
        # create seurat object
        ge <- GEs[[mc_gamma]]
        labels <- rep('ctrl', ncol(ge))
        labels[grep('treat', colnames(ge))] <- 'treat'
        meta <- data.frame(label = labels, row.names = colnames(ge))
        mc_data <- CreateSeuratObject(counts = ge, meta.data = meta)
        Idents(mc_data) <- 'label'
        
        # Run DE
        metacell_markers_manual <- find_markers_bulk(mc_data, stat.test)
        metacell_markers_manual <- arrangeDE(metacell_markers_manual)
        DEs[[as.character(mc_gamma)]] <- metacell_markers_manual
    }
    saveRDS(DEs, file.path(results_folder, "metaGEMarkersManual.rds"))
    message('Done computing MetaCells DE genes with t-test')
}


# ---------------------------------------------------------
# Metacell GE DESeq2
# ---------------------------------------------------------
if(computeMetaDes){
    # Metacell own GE matrix from pipeline, DESeq2
    message('Computing MetaCells DE genes with DESeq2')
    GEs <- createMCGE(sc_filtered_data, results_folder)
    metaInfo <- createMCMembership(sc_filtered_data, results_folder)
    mc_gammas <- names(GEs)
    DEs <- list()
    for(mc_gamma in mc_gammas){
        ge <- GEs[[mc_gamma]]
        membership <- metaInfo$membership[[mc_gamma]]
        sizes <- table(membership)
        labels <- rep('ctrl', ncol(ge))
        labels[grep('treat', colnames(ge))] <- 'treat'
        labels <- data.frame(label = labels, row.names = colnames(ge))
        
        # Run DE
        ge <- floor(sweep(ge, 2, sizes, '*'))
        dds <- DESeqDataSetFromMatrix(ge,
                                      colData = labels, 
                                      design = ~ label)
        
        dds_wald <- DESeq(dds, test = 'Wald', minReplicatesForReplace = Inf)
        results_wald <- results(dds_wald)
        
        DEs[[as.character(mc_gamma)]] <- arrangeDE(results_wald, 
                                                oldNameLog = 'log2FoldChange',
                                                              oldNameP = 'padj')
    }
    saveRDS(DEs, file.path(results_folder, "metaGEMarkersDes.rds"))
    message('Done computing MetaCells DE genes with DESeq2')
}


# ---------------------------------------------------------
# Metacell GE EdgeR
# ---------------------------------------------------------
if(computeMetaEdge){
    # Metacell own GE matrix from pipeline, EdgeR
    message('Computing MetaCells DE genes with EdgeR')
    GEs <- createMCGE(sc_filtered_data, results_folder)
    metaInfo <- createMCMembership(sc_filtered_data, results_folder)
    mc_gammas <- names(GEs)
    DEs <- list()
    for(mc_gamma in mc_gammas){
        ge <- GEs[[mc_gamma]]
        membership <- metaInfo$membership[[mc_gamma]]
        sizes <- table(membership)
        labels <- rep('ctrl', ncol(ge))
        labels[grep('treat', colnames(ge))] <- 'treat'

        # Run DE
        ge <- floor(sweep(ge, 2, sizes, '*'))
        edge <- DGEList(counts = ge, group = labels)
        edge <- calcNormFactors(edge)
        model <- model.matrix(~labels)
        edge <- estimateDisp(edge, model)
        
        # Quasi likelihood test
        fit <- glmQLFit(edge, model)
        DEs[[as.character(mc_gamma)]] <- arrangeDE(glmQLFTest(fit, coef= 2)$table,
                                                   oldNameP = 'PValue')
    }
    saveRDS(DEs, file.path(results_folder, "metaGEMarkersEdge.rds"))
    message('Done computing MetaCells DE genes with EdgeR')
}


# ---------------------------------------------------------
# Metacell via SuperCell t-test
# ---------------------------------------------------------
if(computeMetaSuper){
    # Metacell membership fed into superCell pipeline
    message('Computing MetaCells through SuperCell, t-test')
    metaInfo <- createMCMembership(sc_filtered_data, results_folder)
    mc_gammas <- names(metaInfo$membership)
    DEs <- list()
    sc_filtered_data$cell_name <- colnames(sc_filtered_data)
    for(mc_gamma in mc_gammas){
        sc_data <- subset(sc_filtered_data, subset = cell_name %in% names(metaInfo$membership[[mc_gamma]]))
        super <- createSuperCells(sc_data, 
                                  as.numeric(mc_gamma),
                                  arithmetic = FALSE)
        super$membership <- metaInfo$membership[[mc_gamma]]
        tmp <- data.frame(membership = super$membership,
                          cLine = sapply(metaInfo$annotation[[mc_gamma]], function(x) substr(x, 1, nchar(x) - 1)))
        cell_line <- unique(tmp) %>%
            arrange(membership)
        rownames(cell_line) <- cell_line$membership
        super$cell_line <- cell_line$cLine
        super$supercell_size <- table(metaInfo$membership[[mc_gamma]])
        super$ge <- supercell_GE(GetAssayData(sc_data), super$membership)
        DE <- supercell_FindMarkers(ge = super$GE,
               supercell_size = super$supercell_size,
               clusters = super$cell_line,
               ident.1 = 'treat',
               ident.2 = 'ctrl',
               logfc.threshold = 0,
               only.pos = F,
               do.bootstrapping = F,
               test.use = stat.test)
        DEs[[as.character(mc_gamma)]] <- arrangeDE(DE)
    }
    saveRDS(DEs, file.path(results_folder, "metaSuperMarkers.rds"))
    message('Done computing MetaCells through SuperCell, t-test')
}


# ---------------------------------------------------------
# Metacell via SuperCell DESeq2
# ---------------------------------------------------------
if(computeMetaSuperDes){
    # Metacell membership fed into superCell pipeline
    message('Computing MetaCells through SuperCell, DESeq2')
    metaInfo <- createMCMembership(sc_filtered_data, results_folder)
    mc_gammas <- names(metaInfo$membership)
    DEs <- list()
    sc_filtered_data$cell_name <- colnames(sc_filtered_data)
    for(mc_gamma in mc_gammas){
        sc_data <- subset(sc_filtered_data, subset = cell_name %in% names(metaInfo$membership[[mc_gamma]]))
        super <- createSuperCells(sc_data, 
                                  as.numeric(mc_gamma),
                                  arithmetic = FALSE)
        super$membership <- metaInfo$membership[[mc_gamma]]
        tmp <- data.frame(membership = super$membership,
                          cLine = sapply(metaInfo$annotation[[mc_gamma]], function(x) substr(x, 1, nchar(x) - 1)))
        cell_line <- unique(tmp) %>%
            arrange(membership)
        rownames(cell_line) <- cell_line$membership
        super$cell_line <- cell_line$cLine
        super$supercell_size <- table(metaInfo$membership[[mc_gamma]])
        super$ge <- supercell_GE(GetAssayData(sc_data), super$membership)
        
        ge <- floor(sweep(super$ge, 2, super$supercell_size, '*'))
        dds <- DESeqDataSetFromMatrix(ge,
                                      colData = cell_line, 
                                      design = ~ membership)
        
        dds_wald <- DESeq(dds, test = 'Wald', minReplicatesForReplace = Inf)
        results_wald <- results(dds_wald)
        
        DEs[[as.character(mc_gamma)]] <- arrangeDE(results_wald, 
                                                   oldNameLog = 'log2FoldChange',
                                                   oldNameP = 'padj')
    }
    saveRDS(DEs, file.path(results_folder, "metaSuperMarkersDes.rds"))
    message('Done computing MetaCells through SuperCell, DESeq2')
}


# ---------------------------------------------------------
# Metacell via SuperCell EdgeR
# ---------------------------------------------------------
if(computeMetaSuperEdge){
    # Metacell membership fed into superCell pipeline
    message('Computing MetaCells through SuperCell, EdgeR')
    metaInfo <- createMCMembership(sc_filtered_data, results_folder)
    mc_gammas <- names(metaInfo$membership)
    DEs <- list()
    sc_filtered_data$cell_name <- colnames(sc_filtered_data)
    for(mc_gamma in mc_gammas){
        sc_data <- subset(sc_filtered_data, subset = cell_name %in% names(metaInfo$membership[[mc_gamma]]))
        super <- createSuperCells(sc_data, 
                                  as.numeric(mc_gamma),
                                  arithmetic = FALSE)
        super$membership <- metaInfo$membership[[mc_gamma]]
        tmp <- data.frame(membership = super$membership,
                          cLine = sapply(metaInfo$annotation[[mc_gamma]], function(x) substr(x, 1, nchar(x) - 1)))
        cell_line <- unique(tmp) %>%
            arrange(membership)
        rownames(cell_line) <- cell_line$membership
        super$cell_line <- cell_line$cLine
        super$supercell_size <- table(metaInfo$membership[[mc_gamma]])
        super$ge <- supercell_GE(GetAssayData(sc_data), super$membership)
        
        ge <- floor(sweep(super$ge, 2, super$supercell_size, '*'))
        edge <- DGEList(counts = ge, group = super$cell_line)
        edge <- calcNormFactors(edge)
        model <- model.matrix(~super$cell_line)
        edge <- estimateDisp(edge, model)
        
        # Quasi likelihood test
        fit <- glmQLFit(edge, model)
        DEs[[as.character(mc_gamma)]] <- arrangeDE(glmQLFTest(fit, coef= 2)$table,
                                                   oldNameP = 'PValue')
    }
    saveRDS(DEs, file.path(results_folder, "metaSuperMarkersEdge.rds"))
    message('Done computing MetaCells through SuperCell, EdgeR')
}


# ---------------------------------------------------------
# Random grouping t-test
# ---------------------------------------------------------
if(computeRandomGroup){
    message('Computing Random grouping t-test')
    ge <- GetAssayData(sc_filtered_data)
    samples <- unique(sc_filtered_data$sample)
    DEs <- list()
    for(gamma in gammas){
        tmp <- lapply(samples, function(x) randomGrouping(ge[, which(sc_filtered_data$sample == x)], gamma))
        randomGrp <- tmp[[1]]
        labels <- rep(samples[1], ncol(tmp[[1]]))
        for(i in seq(2:length(tmp))){
            randomGrp <- cbind(randomGrp, tmp[[i]])
            labels <- c(labels, rep(samples[i], ncol(tmp[[i]])))
        }
        colnames(randomGrp) <- as.character(seq(1:ncol(randomGrp)))
        rdmData <- CreateSeuratObject(randomGrp, 
                                      meta.data = data.frame(label = labels, rownames = colnames(randomGrp)))
        Idents(rdmData) <- 'label'
        DE <- find_markers_bulk(rdmData, stat.test)
        DEs[[as.character(gamma)]] <- arrangeDE(DE)
    }
    saveRDS(DEs, file.path(results_folder, 'randomGrouping.rds'))
    message('Done computing Random grouping t-test')
}


# ---------------------------------------------------------
# Random grouping DESeq2
# ---------------------------------------------------------
if(computeRandomGroup){
    memory.limit(size=56000)
    message('Computing Random grouping DESeq2')
    ge <- GetAssayData(sc_filtered_data)
    samples <- unique(sc_filtered_data$sample)
    DEs <- list()
    for(gamma in gammas){
        tmp <- lapply(samples, function(x) randomGrouping(ge[, which(sc_filtered_data$sample == x)], 
                                                          gamma,
                                                          operation = 'sum'))
        randomGrp <- tmp[[1]]
        labels <- rep(samples[1], ncol(tmp[[1]]))
        for(i in seq(2:length(tmp))){
            randomGrp <- cbind(randomGrp, tmp[[i]])
            labels <- c(labels, rep(samples[i], ncol(tmp[[i]])))
        }
        colnames(randomGrp) <- as.character(seq(1:ncol(randomGrp)))
        labels <- sapply(labels, function(x) substr(x, 1, nchar(x) - 1))
        dds <- DESeqDataSetFromMatrix(randomGrp,
                                      colData = data.frame(labels = labels), 
                                      design = ~ labels)
        
        dds_wald <- DESeq(dds, test = 'Wald', minReplicatesForReplace = Inf)
        results_wald <- results(dds_wald)
        
        DEs[[as.character(gamma)]] <- arrangeDE(results_wald, 
                                                oldNameLog = 'log2FoldChange',
                                                oldNameP = 'padj')
    }
    saveRDS(DEs, file.path(results_folder, 'randomGroupingDes.rds'))
    message('Done computing Random grouping DESeq2')
}



# ---------------------------------------------------------
# Random grouping EdgeR
# ---------------------------------------------------------
if(computeRandomGroup){
    message('Computing Random grouping EdgeR')
    ge <- GetAssayData(sc_filtered_data)
    samples <- unique(sc_filtered_data$sample)
    DEs <- list()
    for(gamma in gammas){
        tmp <- lapply(samples, function(x) randomGrouping(ge[, which(sc_filtered_data$sample == x)], 
                                                          gamma,
                                                          operation = 'sum'))
        randomGrp <- tmp[[1]]
        labels <- rep(samples[1], ncol(tmp[[1]]))
        for(i in seq(2:length(tmp))){
            randomGrp <- cbind(randomGrp, tmp[[i]])
            labels <- c(labels, rep(samples[i], ncol(tmp[[i]])))
        }
        colnames(randomGrp) <- as.character(seq(1:ncol(randomGrp)))
        labels <- sapply(labels, function(x) substr(x, 1, nchar(x) - 1))
        edge <- DGEList(counts = randomGrp, group = labels)
        edge <- calcNormFactors(edge)
        model <- model.matrix(~labels)
        edge <- estimateDisp(edge, model)
        
        # Quasi likelihood test
        fit <- glmQLFit(edge, model)
        DEs[[as.character(mc_gamma)]] <- arrangeDE(glmQLFTest(fit, coef= 2)$table,
                                                   oldNameP = 'PValue')
    }
    saveRDS(DEs, file.path(results_folder, 'randomGroupingEdge.rds'))
    message('Done computing Random grouping EdgeR')
}
