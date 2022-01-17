# 5 January 2022

# Run all processing steps for a set of data and saves all processed data
# Everything is run according to the config file provided

# ---------------------------------------------------------
# Header
# ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0){
    stop('You must provide a configuration file', call. = FALSE)
}

# SHOULD BE CHANGED ACCORDINGLY TO LOCATIONS OF R LIBRARIES IF NEEDED
.libPaths("C:/Users/miche/OneDrive/Documents/R/win-library/4.1")


# ---------------------------------------------------------
# Libraries and dependencies
# ---------------------------------------------------------
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(ggrepel)
library(ggExtra)
library(data.table)
library(tidyr)
library(stringr)
library(config)

source('src/utility.R')
source('src/processing.R')

# ---------------------------------------------------------
# Meta parameters
# ---------------------------------------------------------
config <- config::get(file = args[1])

filename <- config$filename
data_folder <- file.path("data", config$intermediaryDataFile)

ctrl_vs_treat <- list(ctrl = config$ctrl_vs_treat$ctrl,
                      treat = config$ctrl_vs_treat$treat)

filtering_param <- list(max.doublet.percentile = config$filteringParam$doubletMaxPercentile,
                        min.gene.per.cell = config$filteringParam$minGenePerCell,
                        min.count.per.cell = config$filteringParam$minCountPerCell,
                        min.count.per.genes = config$filteringParam$minCountPerGene,
                        max.ribo.percent = config$filteringParam$maxRiboPercent,
                        max.mito.percent = config$filteringParam$maxMitoPercent,
                        max.hb.percent = config$filteringParam$maxHbPercent)

normMethod <- config$normMethod

centers <- config$centers

compute_cluter <- config$computeCluster

# ---------------------------------------------------------
# Data loading
# ---------------------------------------------------------
scpath <- 'data/sc_rnaseq/rds'
bulkpath <- 'data/bulk_rnaseq/rds'

sc_data <- readRDS(file.path(scpath, filename))
bulk_raw <- readRDS(file.path(bulkpath, filename))

#Label modification to ensure homogeneity between datasets
sc_data$label <- factor(sc_data$label, 
                        labels = names(ctrl_vs_treat), 
                        levels = ctrl_vs_treat)


rownames(bulk_raw$meta) <- bulk_raw$meta$sample
bulk_data <- CreateSeuratObject(bulk_raw$assay, meta.data = bulk_raw$meta)

bulk_data$label <- factor(tolower(bulk_raw$meta$label), 
                          labels = names(ctrl_vs_treat), 
                          levels = ctrl_vs_treat)

bulk_data$sample <- createSample(bulk_data)
sc_data$sample <- createSample(sc_data)

Idents(bulk_data) <- 'label'
Idents(sc_data) <- 'label'
set.seed(0)


# ---------------------------------------------------------
# QC and filtering
# ---------------------------------------------------------
sc_data <- singleCell_qc(sc_data)
sc_filtered_data <- singleCell_filtering(sc_data, filtering_param)
bulk_filtered_data <- bulk_qc_and_filtering(bulk_data)


# ---------------------------------------------------------
# Intersection (to avoid dropout problems)
# ---------------------------------------------------------
# should intersect common genes between sc and bulk before DEA to ensure they
# are comparable
use_genes <- intersect(rownames(bulk_filtered_data), rownames(sc_filtered_data))
sc_filtered_data <- sc_filtered_data[use_genes, ]
bulk_filtered_data <- bulk_filtered_data[use_genes, ]


# ---------------------------------------------------------
# Normalization
# ---------------------------------------------------------
# Only on single cell data as bulk data should be used with raw counts
sc_normalized_data <- NormalizeObject(sc_filtered_data, method = normMethod)


# ---------------------------------------------------------
# Pseudo bulk creation
# ---------------------------------------------------------
pseudobulk_data <- create_pseudobulk(sc_filtered_data)
pseudobulk_norm <- NormalizeObject(pseudobulk_data, method = normMethod)


# ---------------------------------------------------------
# Clustering of single cells
# ---------------------------------------------------------
# Subgroups may appear in the single cell data, which should be clustered to
# compare these groups between each other
sc_clustered_data <- sub_cluster(sc_normalized_data)
if(compute_cluter){
    sc_clustered_data <- reIdent(sc_clustered_data, centers)
}else{
    # group according to samples and conditions
    Idents(sc_clustered_data) <- 'sample'
}


# bulk data plot
bdata <- NormalizeData(bulk_filtered_data) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 500)

plot1 <- VariableFeaturePlot(bdata)
plot2 <- LabelPoints(plot = plot1, 
                     points = head(VariableFeatures(bdata), 10),
                     repel = T,
                     max.overlaps = Inf)
print(plot2)


# ---------------------------------------------------------
# Data saves
# ---------------------------------------------------------
dir.create(data_folder, showWarnings = F, recursive = T)
saveRDS(sc_clustered_data, file = file.path(data_folder, "singleCellClusteredNormalized.rds"))
saveRDS(sc_filtered_data, file = file.path(data_folder, "singleCellFiltered.rds"))
saveRDS(pseudobulk_data, file = file.path(data_folder, "pseudoBulk.rds"))
saveRDS(pseudobulk_norm, file = file.path(data_folder, "pseudoBulkNormalized.rds"))
saveRDS(bulk_filtered_data, file = file.path(data_folder, "bulkFiltered.rds"))
saveRDS(bdata, file = file.path(data_folder, "bulkFilteredNormalized.rds"))
