# 5th November 2021

# Script to perform DE analysis on cell_lines from Tian et al., 2019 using 
# different techniques such as scRNA or pseudobulk to compare to the results of 
# SuperCell DE analysiss

# TODO: use logr to create logfile


library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)
library(data.table)

source('utility.R')
source('supercells.R')
source('analysis.R')
source('processing.R')

# data loading
load('sincell_with_class_5cl.RData')

# Convert to seurat for DE analysis and pseudobulk build
meta <- data.frame(cell_type = sce_sc_10x_5cl_qc@colData$cell_line,
                   row.names = rownames(sce_sc_10x_5cl_qc@colData),
                   doublet = sce_sc_10x_5cl_qc$demuxlet_cls)

sc <- CreateSeuratObject(counts = sce_sc_10x_5cl_qc@assays$data$counts, 
                         meta.data = meta)

# ---------------------------------------------------------
# QC
# ---------------------------------------------------------

filtered_data <- qc_and_filtering(sc)

# ---------------------------------------------------------
# Clustering
# ---------------------------------------------------------

# #VizDimLoadings(sc, dims = 1:2, reduction = "pca")
#print(DimPlot(sc, reduction = "pca", group.by = 'cell_type'))
#print(DimPlot(sc, reduction = 'tsne', group.by = 'cell_type'))


# ---------------------------------------------------------
# DE
# ---------------------------------------------------------
Idents(filtered_data) <- 'cell_type'

markers <- compute_DE_sc(filtered_data, c('t', 'wilcox', 'negbinom', 'bimod'), 1)

# ---------------------------------------------------------
# Comparison
# ---------------------------------------------------------
# AUC computation
gammas <- c(1, 1.5, 2, 5, 10, 20, 50)
topn = 50
df_results <- run_all_aucc(filtered_data, gammas, markers, 5, topn)
auc_results <- df_results$auc
plot_auc_results(auc_results, topn)

# TPR computation
topn = 20
tpr_results <- compute_tpr(markers1 = markers, markers2 = df_results$markers, size = topn)
tpr_results <- tpr_results %>% 
    rbindlist() %>%
    t() %>%
    data.frame() %>%
    mutate(gamma = gammas)
colnames(tpr_results)[1:(length(tpr_results)-1)] <- names(markers)
plot_tpr_results(tpr_results, topn)


# Visualization
cl <- 'A549'
# show variation of genes between clusters for single cell matrix
#display_significant_genes(sc, cl, super_markers, markers)

# show variation of genes between clusters for super cell matrix
#display_significant_genes(supercells, cl, super_markers, markers, sample = ncol(supercells))
