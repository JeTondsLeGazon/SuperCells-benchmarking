# 5 January 2022

# Run all analysis scrips on specified data
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
library(SuperCell)
library(SuperCellBM)

source('src/utility.R')
source('src/analysis.R')
source('src/supercells.R')

# ---------------------------------------------------------
# Meta parameters
# ---------------------------------------------------------
config <- config::get(file = args[1])

filename <- config$filename
data_folder <- file.path("data", config$intermediaryDataFile)
results_folder <- file.path("data", config$resultsFile)
markers_folder <- paste('markers', config$splitBy, sep = '_')

stat.test <- config$statTest
weighted <- config$weightedSuperCells

gammas <- config$gammas

algos <- config$algos
split.by <- config$splitBy

set.seed(0)


# ---------------------------------------------------------
# Data loadings
# ---------------------------------------------------------
if(!dir.exists(data_folder)){
    stop(sprintf("Cannot load data from folder %s, does not exist", data_folder))
}

if(!dir.exists(results_folder)){
    stop(sprintf("Cannot load data from folder %s, does not exist", results_folder))
}

# Data
sc_data <- readRDS(file = file.path(data_folder, "singleCellData.rds"))
pseudobulk_data <- readRDS(file = file.path(data_folder, "pseudoBulkData.rds"))
bulk_data <- readRDS(file = file.path(data_folder, "bulkData.rds"))

# DE markers
markers.type <- c('bulk', 'super', 'meta', 'metasc', 'random', 'subsampling')
markers <- list()
for(marker in markers.type){
    markers[[marker]] <- load_markers(marker, algos, split.by, results_folder)
}


# ---------------------------------------------------------
# Benchmarking plot to compare different techniques
# ---------------------------------------------------------
# Benchmarking plots to compare supercells, metacells, random grouping and subsampling
# Aims to follow plots given in SuperCell manuscript, written by Mariia
for(algo in algos){
    stat.method <- switch(algo, 'DESeq2' = 'des', 'EdgeR' = 'edge', 't-test' = 't')
    to_compare <- list()
    for(marker in setdiff(markers.type, 'bulk')){
        type_algo <- paste(marker, stat.method, sep = '_')
        to_compare[[type_algo]] <- markers[[marker]][[algo]]
    }

    for(score.method in c('auc', 'tpr', 'tpr_100')){
        plot_results_BM(to_compare, 
                        markers$bulk[[algo]],
                        GT.type = algo,
                        score.type = score.method)
    }
}


# ---------------------------------------------------------
# Volcano Plots
# ---------------------------------------------------------
# Create some volcano plots to observe results for different markers
fig <- volcano_plot(markers$bulk[[algos[1]]])
title <- sprintf('Volcano Plot of bulk for %s', algos[1])
fig <- annotate_figure(fig, top = text_grob(title, face = 'bold', color = 'red', 
                                     size = 14))
print(fig)

fig <- volcano_plot(markers$super[[algos[1]]][[as.character(gammas[1])]])
title <- sprintf('Volcano Plot of SuperCell for %s at gamma = %s', algos[1], gammas[1])
fig <- annotate_figure(fig, top = text_grob(title, face = 'bold', color = 'red', 
                                            size = 14))
print(fig)


# ---------------------------------------------------------
# LogFC - LogFC graphs (with t-test)
# ---------------------------------------------------------
# This figure aims to show the evolution of logFC at different gammas when compared to 
# the single-cell level (gamma = 1). It helped implement the arithmetic mean instead of
# of the geometric one in the computation of the ge of supercells to have constant
# LogFCs whatever the gamma
if('t-test' %in% algos){
    super_markers <- markers$super[[1]]
    p1 <- LogFcLogFcPlot(super_markers$`1`, super_markers$`1`) + 
        ylab('LogFC Super Cells gamma = 1') +
        xlab('LogFC single cells (seurat)')
    p2 <- LogFcLogFcPlot(super_markers$`1`, super_markers$`2`) +
        ylab('LogFC Super Cells gamma = 2') +
        xlab('LogFC single cells (seurat)')
    p3 <- LogFcLogFcPlot(super_markers$`1`, super_markers$`5`) + 
        ylab('LogFC Super Cells gamma = 5') +
        xlab('LogFC single cells (seurat)')
    p4 <- LogFcLogFcPlot(super_markers$`1`, super_markers$`10`) + 
        ylab('LogFC Super Cells gamma = 10') +
        xlab('LogFC single cells (seurat)')
    
    fig <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
    fig <- annotate_figure(fig, top = text_grob('LogFC vs LogFC graph for SuperCells vs single cells', 
                                         face = 'bold', color = 'red', size = 14))
    print(fig)
}


# This figure here shows what happens when we use the geometric mean for the computation
# of the ge of supercells and show the evolution of the logFC
geom_markers <- list()
for(gamma in gammas[c(1,2,3,4)]){
    filename_no_extension <- paste('superCells', gamma, split.by, sep = '_')
    filename <- paste0(filename_no_extension, '.Rds')
    ToComputeSC <- !file.exists(file.path(data_folder, 'SC', filename))
    SC.list <- compute_supercells(
        sc = sc_data,
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
        seed.seq = c(0),
        split.by = split.by
    )
    super <- SC.list$Exact[[as.character(gamma)]][[1]]
    super$label <- supercell_assign(clusters = sc_data$label,
                                    supercell_membership = super$membership,
                                    method = "jaccard")
    
    # Geometric average as we compute average on log normalized counts
    super$GE <- supercell_GE(sc_data@assays$RNA@data, super$membership)
    DE <- compute_supercell_DE(super, 't-test')
    geom_markers[[as.character(gamma)]] <- DE
}
p1 <- LogFcLogFcPlot(geom_markers$`1`, geom_markers$`1`) + 
    ylab('LogFC Super Cells gamma = 1') +
    xlab('LogFC single cells (seurat)')
p2 <- LogFcLogFcPlot(geom_markers$`1`, geom_markers$`2`) +
    ylab('LogFC Super Cells gamma = 2') +
    xlab('LogFC single cells (seurat)')
p3 <- LogFcLogFcPlot(geom_markers$`1`, geom_markers$`5`) + 
    ylab('LogFC Super Cells gamma = 5') +
    xlab('LogFC single cells (seurat)')
p4 <- LogFcLogFcPlot(geom_markers$`1`, geom_markers$`10`) + 
    ylab('LogFC Super Cells gamma = 10') +
    xlab('LogFC single cells (seurat)')

fig <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
fig <- annotate_figure(fig, top = text_grob('LogFC vs LogFC graph for SuperCells (geometric mean) vs single cells', 
                                     face = 'bold', color = 'red', size = 14))
print(fig)


# ---------------------------------------------------------
# LogFC comparison
# ---------------------------------------------------------
if(length(intersect(c('t-test', 'DESeq2', 'EdgeR'), algos)) == 3){
    pairs <- combn(algos, 2)
    order.genes <- markers$bulk$`t-test`$gene
    for(i in 1:ncol(pairs)){
        plot(markers$bulk[[pairs[1, i]]][order.genes, 'logFC'], 
             markers$bulk[[pairs[2, i]]][order.genes, 'logFC'], 
             xlab = pairs[1, i],
             ylab = pairs[2, i],
             main = 'LogFC comparison between statistical packages')
        abline(a = 0, b = 1, col = 'red')
    }
}


# ---------------------------------------------------------
# Rank top n genes
# ---------------------------------------------------------
# We compare the ranking of the n first genes between two methods to show how
# the ranking is affected with different techniques/gammas
# Note that no genes may appear at a certain gamma, as the selected genes (concerned genes)
# are those of supercell at gamma = 1

# Supercells vs Bulk
N <- 100
for(algo in algos){
    super_markers <- markers$super[[algo]]
    concerned_genes <- super_markers$`1`$gene[1:N]
    fig <- rank_plot(concerned_genes, markers$bulk[[algo]], super_markers)
    title <- sprintf('Comparison of Bulk and Supercell, calculated with %s',
                     algo)
    fig <- annotate_figure(fig, top = text_grob(title, color = "red", 
                                         face = "bold", size = 14))
    print(fig)
}

# Supercells vs Metacells at gamma ~70
for(algo in algos){
    super_markers <- markers$super[[algo]]
    concerned_genes <- super_markers$`1`$gene[1:N]
    fig <- rank_plot(concerned_genes, markers$metasc[[algo]][[2]], super_markers)
    title <- sprintf('Comparison of Metacell (gamma = 32) and Supercell, calculated with %s',
                     algo)
    fig <- annotate_figure(fig, top = text_grob(title, color = "red", 
                                                face = "bold", size = 14))
    print(fig)
}

# ---------------------------------------------------------
# Weighted vs unweighted t-test
# ---------------------------------------------------------
# Compare weithed vs unweighted t-test for p-value computation of supercell markers
# At first, weighted t-test was used with supercell size as weitghts. However, the 
# weighted t-test lowers all p values by a huge margin, resulting in some cases in all
# genes being DE. Unweighted t-test was therefore selected for the analysis
if('t-test' %in% algos){
    # Random genes that do not have p-value = 0 for both weighted and unweighted
    selected.genes <-  c('Vps37a', 'Maip1', 'Mitf', 'Calu', 'Xrn1', 'Dnajc1', 'Ppfibp2')
    points_weighted <- c()
    points_unweighted <- c()
    for(gamma in gammas[c(1,2,3,4)]){
        super <- createSuperCellsBM(sc_data, 
                                    gamma = gamma,
                                    data_folder = data_folder,
                                    arithmetic = TRUE,
                                    split.by = split.by,
                                    SC.type = 'Exact',
                                    force_compute = as.logical(config$compute_supercell))
        idx <- super$label == 'treat'
        idy <- super$label == 'ctrl'
        weighted <- apply(super$GE[selected.genes, ], 1, 
                     function(row) weights::wtd.t.test(x = row[idx],
                                                       y = row[idy],
                                                       weight = super$supercell_size[idx],
                                                       weighty = super$supercell_size[idy])$coefficients[['p.value']])
        unweighted <- unlist(markers$super$`t-test`[[as.character(gamma)]][selected.genes, 'p.value'])
        weighted[weighted == 0] <- 1
        points_weighted <- c(points_weighted, weighted)
        points_unweighted <- c(points_unweighted, unweighted)
    }
    N = length(selected.genes)
    K = length(points_weighted) / N
    plot(NULL, xlim = c(0,300), ylim = c(0,300),
         xlab = 'Unweighted p values',
         ylab = 'Weighted p values',
         main = 'Weighted vs unweighted p values')
    colors <- c('blue', 'red', 'green', 'black')
    legends <- c()
    for(k in seq_len(K)){
        points(-log10(points_unweighted[((k-1)*N+1):(k*N)]), 
             -log10(points_weighted[((k-1)*N+1):(k*N)]),
             pch = 16,
             col = colors[k])
        legends <- c(legends, sprintf('Gamma = %s', gammas[k]))
    }
    legend(x="topleft", 
           legend=legends,
           col=colors, lwd=1, lty=c(1,2), 
           pch=16) 
}


# ---------------------------------------------------------
# Fraction
# ---------------------------------------------------------
# Computes the fraction of genes associated with different cluster. It can help
# understand how the fraction of DE genes may vary with different gammas, ie if
# a higher gamma will produce DE genes closely related to Bulk or not.
frac <- fractionGenes(list(single = markers$super$`t-test`$`1`, 
                           super10 = markers$super$`t-test`$`10`, 
                           bulk = markers$bulk$`t-test`))

fracPercent <- sapply(frac, function(x) table(x)/length(x) * 100)
fracPercent <- melt(fracPercent)
ggplot(data = fracPercent, aes(x = Var2, y = value, fill = Var1)) + 
    geom_bar(stat = 'identity', position = 'stack') +
    xlab('') +
    ylab('Percentage') +
    ggtitle('Repartition of DE genes found among different sources')