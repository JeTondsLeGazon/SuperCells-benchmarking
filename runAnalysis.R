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
library(limma)
library(ggpubr)
library(ggrepel)
library(stringr)
library(tidyseurat)
library(ggExtra)
library(weights)
library(zoo)

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
sc_data <- readRDS(file = file.path(data_folder, "singleCellClusteredNormalized.rds"))
pseudobulk_data <- readRDS(file = file.path(data_folder, "pseudoBulkNormalized.rds"))
bulk_data <- readRDS(file = file.path(data_folder, "bulkFilteredNormalized.rds"))

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

    for(score.method in c('auc', 'tpr', 'match')){
        plot_results_BM(to_compare, 
                        markers$bulk[[algo]],
                        GT.type = algo,
                        score.type = score.method)
    }
}


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
    annotate_figure(fig, top = text_grob('LogFC vs LogFC graph for SuperCells vs single cells', 
                                         face = 'bold', color = 'red', size = 14))
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
    selected.genes <- markers$bulk$`t-test`$gene[1:10]
    points <- data.frame(x = c(), y = c(), gamma = c())
    for(gamma in gammas[c(1,2,3,4)]){
        super <- createSuperCellsBM(sc_data, 
                                    gamma = gamma,
                                    results_folder = results_folder,
                                    arithmetic = TRUE,
                                    split.by = split.by,
                                    SC.type = 'Exact',
                                    force_compute = as.logical(config$compute_supercell))
        idx <- super$cell_line == 'treat'
        idy <- super$cell_line == 'ctrl'
        weighted <- apply(super$GE[selected.genes, ], 1, 
                     function(row) weights::wtd.t.test(x = row[idx],
                                                       y = row[idy],
                                                       weight = super$supercell_size[idx],
                                                       weighty = super$supercell_size[idy])$coefficients[['p.value']])
        unweighted <- markers$super$`t-test`[[gamma]]$p.value[selected.genes]
        points$x <- c(points$x, weighted)
        points$y <- c(points$y, unweighted)
        points$gamma <- c(points$gamma, rep(gamma, length(selected.genes)))
    }
}


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