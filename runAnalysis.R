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
bulk_markers <- load_markers('bulk', algos, split.by, results_folder)
super_markers <- load_markers('super', algos, split.by, results_folder)
mc_markers <- load_markers('meta', algos, split.by, results_folder)
mc_sc_markers <- load_markers('metasc', algos, split.by, results_folder)
random_markers <- load_markers('random', algos, split.by, results_folder)
sub_markers <- load_markers('subsampling', algos, split.by, results_folder)


# ---------------------------------------------------------
# Comparison
# ---------------------------------------------------------
su_mc <- list('super_t' = super_markers,
              'meta_t' = mc_markers,
              'metasc_t' = mc_sc_markers, 
              'random_t' = random_markers,
              'sub_t' = sub_sampling)

comp <- list('Bulk (t-test)' = bulk_markers)
for(score.type in c('auc', 'tpr', 'match')){
    plot_results_flex(su_mc, comp, score.type = score.type)
}

su_mc <- list('super_des' = super_markers_des,
              'meta_des' = mc_markers_des,
              'metasc_des' = mc_sc_markers_des, 
              'random_des' = random_markers_des,
              'sub_des' = sub_sampling_des)

comp <- list('Bulk (DESeq2)' = bulk_markers_des)
for(score.type in c('auc', 'tpr', 'match')){
    plot_results_flex(su_mc, comp, score.type = score.type)
}
su_mc <- list('super_edge' = super_markers_edge,
              'meta_edge' = mc_markers_edge,
              'metasc_edge' = mc_sc_markers_edge, 
              'random_edge' = random_markers_edge,
              'sub_edge' = sub_sampling_edge)

comp <- list('Bulk (EdgeR)' = bulk_markers_edge)
for(score.type in c('auc', 'tpr', 'match')){
    plot_results_flex(su_mc, comp, score.type = score.type)
}


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