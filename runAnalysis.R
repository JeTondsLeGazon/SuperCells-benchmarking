# 5 January 2022

# Run all analysis scrips on specified data
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

source('utility.R')
source('analysis.R')
source('supercells.R')

# ---------------------------------------------------------
# Meta parameters
# ---------------------------------------------------------
config <- config::get(file = args[1])

filename <- config$filename
data_folder <- file.path("data", config$intermediaryDataFile)
results_folder <- file.path("data", config$resultsFile)

stat.test <- config$statTest
weighted <- config$weightedSuperCells

gammas <- config$gammas


# ---------------------------------------------------------
# Data loadings
# ---------------------------------------------------------
if(!dir.exists(data_folder)){
    stop(sprintf("Cannot load data from folder %s, does not exist", data_folder))
}

if(!dir.exists(results_folder)){
    stop(sprintf("Cannot load data from folder %s, does not exist", results_folder))
}

sc_clustered_data <- readRDS(file = file.path(data_folder, "singleCellClusteredNormalized.rds"))
sc_filtered_data <- readRDS(file = file.path(data_folder, "singleCellFiltered.rds"))
pseudobulk_data <- readRDS(file = file.path(data_folder, "pseudoBulk.rds"))
pseudobulk_norm <- readRDS(file = file.path(data_folder, "pseudoBulkNormalized.rds"))
bulk_filtered_data <- readRDS(file = file.path(data_folder, "bulkFiltered.rds"))
bdata <- readRDS(file = file.path(data_folder, "bulkFilteredNormalized.rds"))

bulk_markers <- readRDS(file.path(results_folder, "bulkMarkers.rds"))
bulk_markers_manual <- readRDS(file.path(results_folder, "bulkMarkersManual.rds"))
pseudo_markers <- readRDS(file.path(results_folder, "pseudoMarkers.rds"))
pseudo_markers_manual <- readRDS(file.path(results_folder, "pseudoMarkersManual.rds"))
super_markers <- readRDS(file.path(results_folder, "superMarkers.rds"))
super_markers_weighted <- readRDS(file.path(results_folder, "superMarkersWeighted.rds"))
super_markers_des <- readRDS(file.path(results_folder, "superMarkersDes.rds"))
single_markers <- readRDS(file.path(results_folder, "singleMarkers.rds"))
mc_markers <- readRDS(file.path(results_folder, 'metaCellMarkers.rds'))
mc_markers_edge <- readRDS(file.path(results_folder, 'metaCellMarkersEdge.rds'))
mc_markers_des <- readRDS(file.path(results_folder, 'metaCellMarkersDes.rds'))

# ---------------------------------------------------------
# Comparison
# ---------------------------------------------------------
plot_results(mc_markers, 
             list('single cells (t-test)' = single_markers, 
                  'bulk (t-test)' = pseudo_markers_manual, 
                  'bulk (DESeq2)' = bulk_markers$`DESeq2-Wald`,
                  'pseudobulk (DESeq2)' = pseudo_markers$`DESeq2-Wald`,
                  'pseudobulk (t-test)' = pseudo_markers_manual),
             super.type = 'DESeq2',
             score.type = 'tpr')

# ---------------------------------------------------------
# LogFC - LogFC graphs
# ---------------------------------------------------------

p1 <- LogFcLogFcPlot(single_markers, super_markers$`1`) + 
    ylab('LogFC Super Cells gamma = 1') +
    xlab('LogFC single cells (seurat)')
p2 <- LogFcLogFcPlot(single_markers, super_markers$`2`) +
    ylab('LogFC Super Cells gamma = 2') +
    xlab('LogFC single cells (seurat)')
p3 <- LogFcLogFcPlot(single_markers, super_markers$`5`) + 
    ylab('LogFC Super Cells gamma = 5') +
    xlab('LogFC single cells (seurat)')
p4 <- LogFcLogFcPlot(single_markers, super_markers$`10`) + 
    ylab('LogFC Super Cells gamma = 10') +
    xlab('LogFC single cells (seurat)')

fig <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
annotate_figure(fig, top = text_grob('LogFC vs LogFC graph for SuperCells vs single cells', 
                                     face = 'bold', size = 14))


# ---------------------------------------------------------
# Rank top n genes
# ---------------------------------------------------------
concerned_genes <- single_markers$gene[1:100]
rank_plot(concerned_genes, bulk_markers$`DESeq2-Wald`, super_markers)


# ---------------------------------------------------------
# Weighted vs unweighted
# ---------------------------------------------------------
plot(NULL, 
     xlab = 'Unweighted -log10(p values)', 
     ylab = 'Weighted -log10(p values)',
     main = 'P values of t-test vs weighted t-test at different gammas',
     xlim = c(0, 10), 
     ylim = c(0, 10))
colors <- c('red', 'blue', 'green', 'black')
legend_names <- c()
for(i in 1:4){
    unweighted <- super_markers[[gammas[i]]]$adj.p.value
    weighted <- super_markers_weighted[[gammas[i]]]$adj.p.value
    
    row_sub <- which(unweighted != 0 & weighted != 0)
    unweighted <- unweighted[row_sub]
    weighted <- weighted[row_sub]
    
    lines(-log10(unweighted), -log10(weighted), 
          col = colors[i], lty = 1, type = 'l', lwd = 2)
    legend_names <- c(legend_names, sprintf('Gamma = %s', gammas[i]))
}
legend('topright', legend = legend_names, col = colors, lty = 1)

# ---------------------------------------------------------
# Gene selection analysis
# ---------------------------------------------------------
# Selection of genes found as DE by bulk and pseudo bulk -> gt genes that should
# also be found at supercell level
# TODO: update this deprecated functionnality
# 
# df_seurat <- list(single = sc_clustered_data,
#                   pseudo = pseudobulk_norm,
#                   bulk = bdata)
# df_super <- list(super5 = superCells_GE(sc_clustered_data, 5),
#                  super50 = superCells_GE(sc_clustered_data, 50),
#                  super200 = superCells_GE(sc_clustered_data, 200))
# 
# gt.genes <- intersect(bulk_markers$`DESeq2-Wald`$gene, 
#                       pseudo_markers$`DESeq2-Wald`$gene)
# not_in_single <- gt.genes[!(gt.genes %in% single_markers$gene)]
# use.genes <- sample(not_in_single, 10)
# use.genes <- c('Dyrk2')
# # Same but with ggplot format
# all.values <- c()
# all.cells.levels <- c()
# all.grp.levels <- c()
# for(sub in names(df_seurat)){
#     for(cond in c('treat|LPS4', 'ctrl|UNST')){
#         values <- GetAssayData(df_seurat[[sub]])[use.genes[1], grep(cond, Idents(df_seurat[[sub]]))]
#         cell.levels <- rep(sub, length(values))
#         grp.levels <- rep(strsplit(cond, '|', fixed = T)[[1]][1], length(values))
#         all.values <- c(all.values, values)
#         all.cells.levels <- c(all.cells.levels, cell.levels)
#         all.grp.levels <- c(all.grp.levels, grp.levels)
#     }
# }
# 
# for(sub in names(df_super)){
#     for(cond in c('treat|LPS4', 'ctrl|UNST')){
#         values <- df_super[[sub]]$GE[match(use.genes[1], rownames(sc_clustered_data)), grep(cond, df_super[[sub]]$SC.cell.annotation.)]
#         cell.levels <- rep(sub, length(values))
#         grp.levels <- rep(strsplit(cond, '|', fixed = T)[[1]][1], length(values))
#         all.values <- c(all.values, values)
#         all.cells.levels <- c(all.cells.levels, cell.levels)
#         all.grp.levels <- c(all.grp.levels, grp.levels)
#     }
# }
# df <- data.frame(values = all.values, grp = all.grp.levels, size = all.cells.levels)
# txt_df <- data.frame(size = unique(df$size), 
#                      stats = matrix(c(single_markers[use.genes[1], c('logFC', 'adj.p.value')],
#                                       pseudo_markers$`DESeq2-Wald`[use.genes[1], c('logFC', 'adj.p.value')],
#                                       bulk_markers$`DESeq2-Wald`[use.genes[1], c('logFC', 'adj.p.value')],
#                                       super_markers$`5`[use.genes[1], c('logFC', 'adj.p.value')],
#                                       super_markers$`50`[use.genes[1], c('logFC', 'adj.p.value')],
#                                       super_markers$`200`[use.genes[1], c('logFC', 'adj.p.value')]), byrow = T, ncol = 2),
#                      vals = rep(5, length(unique(df$size))))
# level_order <- c('single', 'super5', 'super50', 'super200', 'pseudo', 'bulk')
# ggplot(data = df, aes(x = factor(size, level = level_order), y = values, fill = grp)) +
#     geom_boxplot() +
#     xlab('') +
#     scale_y_log10() +
#     ylab('log of logcounts') +
#     ggtitle(use.genes[1]) +
#     theme_classic() +
#     annotate('text', x = factor(txt_df$size, level = level_order), y = txt_df$vals, label = paste0('logFC: ', round(unlist(txt_df$stats.1), 2), '\n', 'p value: ', format(txt_df$stats.2, scientific = T))) +
#     theme(axis.text=element_text(size=12),
#           axis.title=element_text(size=14,face="bold"))

# ---------------------------------------------------------
# Fraction
# ---------------------------------------------------------
frac <- fractionGenes(list(single = single_markers, 
                           super10 = super_markers$`5`, 
                           super100 = super_markers$`100`, 
                           bulk = bulk_markers$`DESeq2-Wald`))

fracPercent <- sapply(frac, function(x) table(x)/length(x) * 100)
fracPercent <- lapply(seq_along(fracPercent), function(i) setNames(data.frame(fracPercent[i]), 
                                                                   c('n', names(fracPercent[i]))))

mydf <- fracPercent[[1]]
for(i in 2:length(fracPercent))
{
    mydf <- merge(mydf, fracPercent[[i]], by = 'n', all = T)
}
rownames(mydf) <- mydf$n
mydf$n <- NULL
mydf <- mydf %>%
    as.matrix() %>%
    melt()
ggplot(data = mydf, aes(x = Var2, y = value, fill = Var1)) + 
    geom_bar(stat = 'identity', position = 'stack') +
    xlab('') +
    ylab('Percentage') +
    ggtitle('Repartition of DE genes found among different sources')