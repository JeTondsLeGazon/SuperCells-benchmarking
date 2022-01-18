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
    bulk_markers <- compute_DE_bulk(bulk_filtered_data)
    volcano_plot(bulk_markers$`edgeR-QLF`) +
        ggtitle('Volcano plot of bulk data from DESeq2 (wald)') +
        theme(plot.title = element_text(hjust = 0.5))
    
    bulk_markers <- lapply(bulk_markers, function(x) x %>% subset(logFC > 0))
    bulk_markers <- lapply(bulk_markers, function(x) if(nrow(x) == 0){x <- NULL}else{x})
    bulk_markers <- bulk_markers[!unlist(lapply(bulk_markers, is.null))]
    saveRDS(bulk_markers, file.path(results_folder, "bulkMarkers.rds"))
}



# ---------------------------------------------------------
# DE bulk manual
# ---------------------------------------------------------
if(computeBulkManual){
    manual_bulk_markers <- find_markers_bulk(bulk_filtered_data, stat.test) %>%
        arrange(adj.p.value, 1 / (abs(logFC) + 1)) %>%
        mutate(gene = row.names(.)) %>%
        subset(logFC > 0)
    saveRDS(manual_bulk_markers, file.path(results_folder, "bulkMarkersManual.rds"))
}


# ---------------------------------------------------------
# DE pseudobulk DESeq2
# ---------------------------------------------------------
if(computePseudo){
    pseudo_markers <- compute_DE_bulk(pseudobulk_data)
    volcano_plot(pseudo_markers$`edgeR-QLF`) +
        ggtitle('Volcano plot of pseudo bulk data from DESeq2 (wald)') +
        theme(plot.title = element_text(hjust = 0.5))
    
    pseudo_markers <- lapply(pseudo_markers, function(x) subset(x, logFC > 0))
    pseudo_markers <- lapply(pseudo_markers, function(x) if(nrow(x) == 0){x <- NULL}else{x})
    pseudo_markers <- pseudo_markers[!unlist(lapply(pseudo_markers, is.null))]
    saveRDS(pseudo_markers, file.path(results_folder, "pseudoMarkers.rds"))
}

# ---------------------------------------------------------
# DE pseudobulk manual
# ---------------------------------------------------------
if(computePseudoManual){
    pseudo_markers_manual <- find_markers_bulk(pseudobulk_norm, stat.test) %>%
        arrange(adj.p.value, 1 / (abs(logFC) + 1)) %>%
        mutate(gene = row.names(.)) %>%
        subset(logFC > 0)
    saveRDS(pseudo_markers_manual, file.path(results_folder, "pseudoMarkersManual.rds"))
}

# ---------------------------------------------------------
# DE supercells
# ---------------------------------------------------------
if(computeSuper){
    memory.limit(size=56000)
    super_markers <- superCells_DEs(sc_clustered_data, gammas, 5,
                                    weighted = F,
                                    test.use = stat.test)
    
    volcano_plot(super_markers$`1`, logfc.thres = 0.5) +
        ggtitle('Volcano plot of SuperCells at level gamma = 5') +
        theme(plot.title = element_text(hjust = 0.5))
    
    super_markers <- lapply(super_markers, function(x) subset(x, logFC > 0))
    saveRDS(super_markers, file.path(results_folder, "superMarkers.rds"))
    
    super_markers_weighted <- superCells_DEs(sc_clustered_data, gammas, 5,
                                    weighted = T,
                                    test.use = stat.test)
    
    super_markers_weighted <- lapply(super_markers_weighted, 
                                     function(x) subset(x, logFC > 0))
    saveRDS(super_markers_weighted, file.path(results_folder, "superMarkersWeighted.rds"))
}


# ---------------------------------------------------------
# DE supercells DESeq2
# ---------------------------------------------------------
if(computeSuperDes){
    super_markers_des <- list()
    
    for(gamma in gammas[c(-1,-2)]){
        super <-  SCimplify(GetAssayData(sc_filtered_data),
                            cell.annotation = sc_filtered_data$sample,
                            k.knn = 5,
                            gamma = gamma,
                            n.var.genes = 1000,
                            directed = FALSE
        )
        
        super$cell_line <- supercell_assign(clusters = sc_filtered_data$sample,
                                            supercell_membership = super$membership,
                                            method = "jaccard")
        
        super$GE <- supercell_GE(GetAssayData(sc_filtered_data), super$membership)
        colData <- rep('ctrl', ncol(super$GE))
        colData[grep('treat', super$cell_line)] <- 'treat'
        super$design <- data.frame(colData)
        colnames(super$design) <- 'design'
        dds <- DESeqDataSetFromMatrix(floor(sweep(super$GE, 2, super$supercell_size, '*')),
                                      colData = super$design, 
                                      design = ~ design)
        
        dds_wald <- DESeq(dds, test = 'Wald', minReplicatesForReplace = Inf)
        
        
        
        results_wald <- results(dds_wald)
        
        super_markers_des[[as.character(gamma)]] <- as.data.frame(results_wald) %>%
            dplyr::rename(logFC = log2FoldChange, adj.p.value = padj) %>% 
            mutate(gene = rownames(.))
    }
    
    
    super_markers_des <- lapply(super_markers_des, function(x) x %>% 
                                    arrange(adj.p.value, 1/(abs(logFC) + 1)) %>%
                                    subset(logFC > 0))
    saveRDS(super_markers_des, file.path(results_folder, "superMarkersDes.rds"))
}


# ---------------------------------------------------------
# DE supercells EdgeR
# ---------------------------------------------------------
if(computeSuperEdge){
    for(gamma in gammas[c(-1,-2)]){
        message(sprintf('Running for gamma = %s', gamma))
        super <-  SCimplify(GetAssayData(sc_filtered_data),
                            cell.annotation = sc_filtered_data$sample,
                            k.knn = 5,
                            gamma = gamma,
                            n.var.genes = 1000,
                            directed = FALSE
        )
        
        super$cell_line <- supercell_assign(clusters = sc_filtered_data$sample,
                                            supercell_membership = super$membership,
                                            method = "jaccard")
        
        super$GE <- supercell_GE(GetAssayData(sc_filtered_data), super$membership)
        
        # remove last char from sample to have conditions
        design <- sapply(super$cell_line, function(x) substr(x, 1, nchar(x) - 1))
        ge <- floor(sweep(super$GE, 2, super$supercell_size, '*'))
        edge <- DGEList(counts = ge, group = design)
        edge <- calcNormFactors(edge)
        model <- model.matrix(~design)
        edge <- estimateDisp(edge, model)
        
        # Quasi likelihood test
        fit <- glmQLFit(edge, model)
        qlf <- glmQLFTest(fit, coef= 2)$table %>%
            data.frame() %>%
            dplyr::rename(adj.p.value = PValue) %>% 
            mutate(gene = rownames(.))  %>%
            arrange(adj.p.value, 1/(abs(logFC) + 1)) %>%
            subset(logFC > 0)
        message(sprintf('Done for gamma = %s', gamma))
        saveRDS(qlf, file.path(results_folder, sprintf("superMarkersEdge%s.rds", gamma)))
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
}


# ---------------------------------------------------------
# DE single cells
# ---------------------------------------------------------
if(computeSingle){
    single_markers <- singleCell_DE(sc_clustered_data, 
                                    var.features = 500,
                                    stat.test)
    volcano_plot(single_markers, logfc.thres = 0.5) +
        ggtitle('Volcano plot of single cells from FindAllMarkers (seurat)') +
        theme(plot.title = element_text(hjust = 0.5))
    
    single_markers <- single_markers %>% subset(logFC > 0)
    saveRDS(single_markers, file.path(results_folder, "singleMarkers.rds"))
}

# ---------------------------------------------------------
# DE single cells (manual)
# ---------------------------------------------------------
if(computeSingleManual){
    manual_single_markers <- find_markers_bulk(sc_clustered_data, stat.test)
    manual_single_markers <- manual_single_markers %>% 
        arrange(adj.p.value, 1 / (abs(logFC) + 1), T) %>%
        mutate(gene = rownames(.)) %>%
        subset(logFC > 0)
    saveRDS(manual_single_markers, file.path(results_folder, "singleMarkersManual.rds"))
    
}


# ---------------------------------------------------------
# Metacell (unix or mac systems only)
# ---------------------------------------------------------
# should run metacells composition first

if(F){
    stop <- F
    myMcFile <- file.path(results_folder, 'mcComposition.rds')
    if(!file.exists(myMcFile)){
        message('Could not find MetaCell cell composition, cannot run DE analysis on MetaCells')
        stop <- T
    }
    if(stop){break}
    mcComposition <- readRDS(file = myMcFile)
    mcComposition <- mcComposition[!unlist(lapply(mcComposition, is.null))]
    for(sample in names(mcComposition)){
        mcComposition[[sample]] <- mcComposition[[sample]][!unlist(lapply(mcComposition[[sample]], is.null))]
    }
    
    nbGammas <- length(mcComposition[[1]])
    DEs <- list()
    DEs_des <- list()
    DEs_edge <- list()
    for(gam in c(1, nbGammas)){
        mc_membership <- c()
        mc_sizes <- c()
        curMax <- 0
        
        for(tmp in mcComposition){
            mc_membership <- c(mc_membership, tmp[[gam]] + curMax)
            mc_sizes <- c(mc_sizes, as.numeric(table(tmp[[gam]])))
            curMax <- curMax + max(tmp[[gam]] )
        }
        mydata <- sc_clustered_data
        usedCells <- match(names(mc_membership), names(mydata$sample))
        message(sprintf('Dropped %s cells compared to standard supercells for mc to work', abs(length(usedCells) - ncol(mydata))))
        super <-  SCimplify(GetAssayData(mydata)[ ,usedCells],
                            cell.annotation = mydata$sample[usedCells],
                            k.knn = 5,
                            gamma = 50,
                            n.var.genes = 1000,
                            directed = FALSE
        )
        
        super$membership <- mc_membership
        super$cell_line <- supercell_assign(clusters = mydata$label[usedCells],
                                            supercell_membership = super$membership,
                                            method = "jaccard")
        super$supercell_size <- mc_sizes
        super$GE_t <- log(supercell_GE(expm1(GetAssayData(mydata)[ ,usedCells]), super$membership) + 1)
        super$GE_des <- supercell_GE(GetAssayData(sc_filtered_data)[ ,usedCells], super$membership)
        
        DE <- supercell_FindMarkers(ge = super$GE_t,
                     supercell_size = super$supercell_size,
                     clusters = super$cell_line,
                     ident.1 = 'treat',
                     ident.2 = 'ctrl',
                     logfc.threshold = 0,
                     only.pos = F,
                     do.bootstrapping = F,
                     test.use = 't') %>%
            mutate(gene = rownames(.)) %>%
            arrange(adj.p.value, 1 / (abs(logFC) + 1), T) %>%
            subset(!duplicated(.))
        rownames(DE) <- DE$gene
        DEs[[as.character(mean(super$supercell_size))]] <- DE
        
        # ----------- DESeq2 ---------------------
        super$design <- data.frame(super$cell_line)
        colnames(super$design) <- 'design'
        m <- matrix(as.numeric(floor(sweep(as.matrix(super$GE_des), 2, super$supercell_size, '*'))), ncol = ncol(super$GE_des), byrow = T) + 1
        rownames(m) <- rownames(super$GE_des)
        dds <- DESeqDataSetFromMatrix(m,
                                      colData = super$design, 
                                      design = ~ design)
        #dds <- estimateSizeFactors(dds, type = 'iterate')
        dds_wald <- DESeq(dds, test = 'Wald', minReplicatesForReplace = Inf)
        results_wald <- results(dds_wald)
        
        DEs_des[[as.character(mean(super$supercell_size))]] <- as.data.frame(results_wald) %>%
            dplyr::rename(logFC = log2FoldChange, adj.p.value = padj) %>% 
            mutate(gene = rownames(.)) %>%
            arrange(adj.p.value, 1/(abs(logFC) + 1)) %>%
            subset(logFC > 0)
        message('Done with DESeq2')
        
    }
    saveRDS(DEs, file.path(results_folder, "metaCellMarkers.rds"))
    saveRDS(DEs_des, file.path(results_folder, "metaCellMarkersDes.rds"))
}

if(computeMeta){
    stop <- F
    myMcFile <- file.path(results_folder, 'mcComposition.rds')
    if(!file.exists(myMcFile)){
        message('Could not find MetaCell cell composition, cannot run DE analysis on MetaCells')
        stop <- T
    }
    if(stop){break}
    mcComposition <- readRDS(file = myMcFile)
    mcComposition <- mcComposition[!unlist(lapply(mcComposition, is.null))]
    for(sample in names(mcComposition)){
        mcComposition[[sample]] <- mcComposition[[sample]][!unlist(lapply(mcComposition[[sample]], is.null))]
    }
    
    nbGammas <- length(mcComposition[[1]])
    DEs <- list()
    DEs_des <- list()
    DEs_edge <- list()
    for(gam in c(1, nbGammas)){
        mc_membership <- c()
        mc_sizes <- c()
        curMax <- 0
        
        for(tmp in mcComposition){
            mc_membership <- c(mc_membership, tmp[[gam]] + curMax)
            mc_sizes <- c(mc_sizes, as.numeric(table(tmp[[gam]])))
            curMax <- curMax + max(tmp[[gam]] )
        }
        mydata <- sc_clustered_data
        usedCells <- match(names(mc_membership), names(mydata$sample))
        message(sprintf('Dropped %s cells compared to standard supercells for mc to work', abs(length(usedCells) - ncol(mydata))))
        super <-  SCimplify(GetAssayData(mydata)[ ,usedCells],
                            cell.annotation = mydata$sample[usedCells],
                            k.knn = 5,
                            gamma = 50,
                            n.var.genes = 1000,
                            directed = FALSE
        )
        
        super$membership <- mc_membership
        super$cell_line <- supercell_assign(clusters = mydata$label[usedCells],
                                            supercell_membership = super$membership,
                                            method = "jaccard")
        super$supercell_size <- mc_sizes
        super$GE_des <- supercell_GE(GetAssayData(sc_filtered_data)[ ,usedCells], super$membership)
        
        # ----------- edgeR ---------------------
        m <- matrix(as.numeric(floor(sweep(as.matrix(super$GE_des), 2, super$supercell_size, '*'))), ncol = ncol(super$GE_des), byrow = T) + 1
        rownames(m) <- rownames(super$GE_des)
        # edgeR approach
        design <- sapply(super$cell_line, function(x) substr(x, 1, nchar(x) - 1))
        edge <- DGEList(counts = m, group = design)
        edge <- calcNormFactors(edge)
        model <- model.matrix(~design)
        edge <- estimateDisp(edge, model)
        
        # Quasi likelihood test
        fit <- glmQLFit(edge, model)
        qlf <- glmQLFTest(fit,coef= 2)$table %>%
            data.frame() %>%
            dplyr::rename(adj.p.value = PValue) %>% 
            mutate(gene = rownames(.))  %>%
            arrange(adj.p.value, 1/(abs(logFC) + 1)) %>%
            subset(logFC > 0)
        DEs_edge[[as.character(mean(super$supercell_size))]] <- qlf
        message('Done with EdgeR')
        
    }
    saveRDS(DEs_edge, file.path(results_folder, "metaCellMarkersEdge.rds"))
}

