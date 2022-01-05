# 11th November 2021

# Performs comparison between DEA from different sources, mainly SuperCells and
# ground truth (bulk DNA) from various dataset
if(!('SuperCells-benchmarking' %in% unlist(strsplit(getwd(), '/')))){
    setwd(paste0(getwd(), '/SuperCells-benchmarking'))
}
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)
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
source('supercells.R')
source('analysis.R')
source('processing.R')


config <- config::get(file = 'hagai_mouse_lps_config.yml')

# ---------------------------------------------------------
# Meta parameters
# ---------------------------------------------------------
filename <- config$filename
data_folder <- file.path("data", config$intermediaryDataFile)
results_folder <- file.path("data", config$resultsFile)

ctrl_vs_treat <- list(ctrl = config$ctrl_vs_treat$ctrl,
                      treat = config$ctrl_vs_treat$treat)

stat.test <- config$statTest
weighted <- config$weightedSuperCells
resetSingle <- F
resetSuper <- F
resetSuperDes <- T

filtering_param <- list(max.doublet.percentile = config$filteringParam$doubletMaxPercentile,
                        min.gene.per.cell = config$filteringParam$minGenePerCell,
                        min.count.per.cell = config$filteringParam$minCountPerCell,
                        min.count.per.genes = config$filteringParam$minCountPerGene,
                        max.ribo.percent = config$filteringParam$maxRiboPercent,
                        max.mito.percent = config$filteringParam$maxMitoPercent,
                        max.hb.percent = config$filteringParam$maxHbPercent)

normMethod <- config$normMethod
#filtering_param_others <- list(max.doublet.percentile = 0.95,
#                               min.gene.per.cell = 300,
#                               min.count.per.cell = 500,
#                               min.count.per.genes = 400,
#                               max.ribo.percent = 40,
#                               max.mito.percent = 20,
#                               max.hb.percent = 5)

centers <- config$centers

compute_cluter <- config$computeCluster

gammas <- config$gammas

# ---------------------------------------------------------
# Data loading
# ---------------------------------------------------------
#available: Hagai2018_mouse-lps, Hagai2018_mouse-pic, Hagai2018_pig-lps, 
# Hagai2018_rabbit-lps, Hagai2018_rat-lps, Hagai2018_rat-pic, Angelidis2019_pneumo,
# Angelidis2019_alvmac, CanoGamez2020_naive-iTreg, CanoGamez2020_memory-Th17

scpath <- '../sc_rnaseq/rds'
bulkpath <- '../bulk_rnaseq/rds'

sc_data <- readRDS(file.path(scpath, filename))
bulk_raw <- readRDS(file.path(bulkpath, filename))

#Label modification to ensure homogeneity between datasets
sc_data$label <- factor(sc_data$label, 
                        labels = names(ctrl_vs_treat), 
                        levels = ctrl_vs_treat)


rownames(bulk_raw$meta) <- bulk_raw$meta$sample
bulk_data <- CreateSeuratObject(bulk_raw$assay, meta.data = bulk_raw$meta)

# if CanoGamez, troubles in the meta ...
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
# Effects of normalization:
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
sc_data <- singleCell_qc(sc_data)
sc_filtered_data <- singleCell_filtering(sc_data, 
                                         filtering_param)
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
# ------------------------------------------------------------------------------
# CanoGamez2020_memory-Th17: matrix(c(-10, 10, -2, 2), ncol = 2)
# mouse-lps : matrix(c(5, -5, 5, -5, 5, 0, -10, -10), ncol = 2)
# pig-lps: matrix(c(-5, 5, 0, 0), ncol = 2)
# rat-lps: matrix(c(-5, -5, 7, 7,   -5, 5, 5, -5), ncol = 2)
# rabbit-lps: matrix(c(7, -5, 7, -5, -5, -5, 7, 7), ncol = 2)
if(compute_cluter){
    sc_clustered_data <- reIdent(sc_clustered_data, centers)
}else{
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


# ---------------------------------------------------------
# Data loadings
# ---------------------------------------------------------
sc_clustered_data <- readRDS(file = file.path(data_folder, "singleCellClusteredNormalized.rds"))
sc_filtered_data <- readRDS(file = file.path(data_folder, "singleCellFiltered.rds"))
pseudobulk_data <- readRDS(file = file.path(data_folder, "pseudoBulk.rds"))
pseudobulk_norm <- readRDS(file = file.path(data_folder, "pseudoBulkNormalized.rds"))
bulk_filtered_data <- readRDS(file = file.path(data_folder, "bulkFiltered.rds"))
bdata <- readRDS(file = file.path(data_folder, "bulkFilteredNormalized.rds"))

# ---------------------------------------------------------
# DE bulk
# ---------------------------------------------------------

bulk_markers <- compute_DE_bulk(bulk_filtered_data)
volcano_plot(bulk_markers$`edgeR-QLF`) +
    ggtitle('Volcano plot of bulk data from DESeq2 (wald)') +
    theme(plot.title = element_text(hjust = 0.5))

bulk_markers <- lapply(bulk_markers, function(x) x %>% subset(logFC > 0))
bulk_markers <- lapply(bulk_markers, function(x) if(nrow(x) == 0){x <- NULL}else{x})
bulk_markers <- bulk_markers[!unlist(lapply(bulk_markers, is.null))]


# ---------------------------------------------------------
# DE bulk manual
# ---------------------------------------------------------

manual_bulk_markers <- find_markers_bulk(bulk_filtered_data, stat.test) %>%
    arrange(adj.p.value, 1 / (abs(logFC) + 1)) %>%
    mutate(gene = row.names(.)) %>%
    subset(logFC > 0)


# ---------------------------------------------------------
# DE pseudobulk DESeq2
# ---------------------------------------------------------
pseudo_markers2 <- compute_DE_bulk(pseudobulk_data)
volcano_plot(pseudo_markers2$`edgeR-QLF`) +
    ggtitle('Volcano plot of pseudo bulk data from DESeq2 (wald)') +
    theme(plot.title = element_text(hjust = 0.5))

pseudo_markers2 <- lapply(pseudo_markers2, function(x) subset(x, logFC > 0))
pseudo_markers2 <- lapply(pseudo_markers2, function(x) if(nrow(x) == 0){x <- NULL}else{x})
pseudo_markers2 <- pseudo_markers2[!unlist(lapply(pseudo_markers2, is.null))]


# ---------------------------------------------------------
# DE pseudobulk manual
# ---------------------------------------------------------
pseudo_markers_manual <- find_markers_bulk(pseudobulk_norm, stat.test) %>%
    arrange(adj.p.value, 1 / (abs(logFC) + 1)) %>%
    mutate(gene = row.names(.)) %>%
    subset(logFC > 0)


# ---------------------------------------------------------
# DE supercells
# ---------------------------------------------------------
memory.limit(size=56000)
super_markers <- superCells_DEs(sc_clustered_data, gammas, 5, 
                                resetData = resetSuper,
                                weighted = weighted,
                                test.use = stat.test)

volcano_plot(super_markers$`1`, logfc.thres = 0.5) +
    ggtitle('Volcano plot of SuperCells at level gamma = 5') +
    theme(plot.title = element_text(hjust = 0.5))

super_markers <- lapply(super_markers, function(x) subset(x, logFC > 0))

# ---------------------------------------------------------
# DE supercells DESeq2
# ---------------------------------------------------------
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


# ---------------------------------------------------------
# DE single cells
# ---------------------------------------------------------
# Single cell markers

single_markers <- singleCell_DE(sc_clustered_data, var.features = 500, 
                                resetData = resetSingle, stat.test)
volcano_plot(single_markers, logfc.thres = 0.5) +
    ggtitle('Volcano plot of single cells from FindAllMarkers (seurat)') +
    theme(plot.title = element_text(hjust = 0.5))

single_markers <- single_markers %>% subset(logFC > 0)

# ---------------------------------------------------------
# DE single cells (manual)
# ---------------------------------------------------------
manual_single_markers <- find_markers_bulk(sc_clustered_data, stat.test)
manual_single_markers <- manual_single_markers %>% 
    arrange(adj.p.value, 1 / (abs(logFC) + 1), T) %>%
    mutate(gene = rownames(.)) %>%
    subset(logFC > 0)


# ---------------------------------------------------------
# saving
# ---------------------------------------------------------
dir.create(results_folder, showWarnings = F, recursive = T)
saveRDS(single_markers, file.path(results_folder, "singleMarkers.rds"))
saveRDS(super_markers, file.path(results_folder, "superMarkers.rds"))
saveRDS(super_markers_des, file.path(results_folder, "superMarkersDes.rds"))
saveRDS(pseudo_markers2, file.path(results_folder, "pseudoMarkers.rds"))
saveRDS(pseudo_markers_manual, file.path(results_folder, "pseudoMarkersManual.rds"))
saveRDS(bulk_markers, file.path(results_folder, "bulkMarkers.rds"))
saveRDS(manual_bulk_markers, file.path(results_folder, "bulkMarkersManual.rds"))


# ---------------------------------------------------------
# Comparison
# ---------------------------------------------------------
score_results <- compute_score(super_markers, manual_bulk_markers, 'aucc') %>%
    melt %>%
    mutate(gammas = rep(gammas[seq_along(super_markers)], 3))
plot_score_results(score_results)
score_results <- compute_score(super_markers, manual_bulk_markers, 'tpr')

eff <- 'DESeq2'

plot_matches(super_markers_des, 
             list(single = single_markers, 
               bulk_man = manual_bulk_markers, 
               bulk_des = bulk_markers$`DESeq2-Wald`,
               pseudo_des = pseudo_markers2$`DESeq2-Wald`,
               pseudo_man = pseudo_markers_manual),
             c(sprintf('SuperCells (%s-test) vs single cells (%s-test)', eff, stat.test), 
               sprintf('SuperCells (%s-test) vs Bulk (%s-test)', eff, stat.test),
               sprintf('SuperCells (%s-test) vs Bulk (DESeq2)', eff),
               sprintf('SuperCells (%s-test) vs pseudo-bulk (DESeq2)', eff),
               sprintf('SuperCells (%s-test) vs pseudo-bulk (%s test)', eff, stat.test)),
             'match')

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
# Statistics comparison at different gammas
# ---------------------------------------------------------
compare_statistics(super_markers)


# ---------------------------------------------------------
# Weighted vs unweighted (should uncomment lines in super_DE)
# ---------------------------------------------------------
l <- nrow(super_markers$`1`)
df <- data.frame(p1 = super_markers$`1`$adj.p.value[1:l],
                 wp1 = super_markers$`1`$weighted_p[1:l],
                 p2 = super_markers$`2`$adj.p.value[1:l],
                 wp2 = super_markers$`2`$weighted_p[1:l],
                 p5 = super_markers$`5`$adj.p.value[1:l],
                 wp5 = super_markers$`5`$weighted_p[1:l],
                 p10 = super_markers$`10`$adj.p.value[1:l],
                 wp10 = super_markers$`10`$weighted_p[1:l])

plot(-log10(df$p1), -log10(df$wp1), col = 'red', lty = 1, type = 'l', lwd = 2,
     xlab = 'Unweighted -log10(p values)', ylab = 'Weighted -log10(p values)',
     main = 'P values of t-test vs weighted t-test at different gammas')
lines(-log10(df$p2),  -log10(df$wp2), col = 'blue', lty = 1, lwd = 2)
lines(-log10(df$p5),  -log10(df$wp5), col = 'green', lty = 1, lwd = 2)
lines(-log10(df$p10),  -log10(df$wp10), col = 'black', lty = 1, lwd = 2)
legend('topright', c('Gamma = 1', 'Gamma = 2', 'Gamma = 5', 'Gamma = 10'))

# ---------------------------------------------------------
# Gene selection analysis
# ---------------------------------------------------------
# Selection of genes found as DE by bulk and pseudo bulk -> gt genes that should
# also be found at supercell level


df_seurat <- list(single = sc_clustered_data,
          pseudo = pseudobulk_norm,
          bulk = bdata)
df_super <- list(super5 = superCells_GE(sc_clustered_data, 5),
                 super50 = superCells_GE(sc_clustered_data, 50),
                 super200 = superCells_GE(sc_clustered_data, 200))#,
                 #super1000 = superCells_GE(sc_clustered_data, 1000))


gt.genes <- intersect(bulk_markers$`DESeq2-Wald`$gene, 
                      pseudo_markers2$`DESeq2-Wald`$gene)
not_in_single <- gt.genes[!(gt.genes %in% single_markers$gene)]
use.genes <- sample(not_in_single, 10)
use.genes <- c('Dyrk2')
# Same but with ggplot format
all.values <- c()
all.cells.levels <- c()
all.grp.levels <- c()
for(sub in names(df_seurat)){
    for(cond in c('treat|LPS4', 'ctrl|UNST')){
        values <- GetAssayData(df_seurat[[sub]])[use.genes[1], grep(cond, Idents(df_seurat[[sub]]))]
        cell.levels <- rep(sub, length(values))
        grp.levels <- rep(strsplit(cond, '|', fixed = T)[[1]][1], length(values))
        all.values <- c(all.values, values)
        all.cells.levels <- c(all.cells.levels, cell.levels)
        all.grp.levels <- c(all.grp.levels, grp.levels)
    }
}

for(sub in names(df_super)){
    for(cond in c('treat|LPS4', 'ctrl|UNST')){
        values <- df_super[[sub]]$GE[match(use.genes[1], rownames(sc_clustered_data)), grep(cond, df_super[[sub]]$SC.cell.annotation.)]
        cell.levels <- rep(sub, length(values))
        grp.levels <- rep(strsplit(cond, '|', fixed = T)[[1]][1], length(values))
        all.values <- c(all.values, values)
        all.cells.levels <- c(all.cells.levels, cell.levels)
        all.grp.levels <- c(all.grp.levels, grp.levels)
    }
}
df <- data.frame(values = all.values, grp = all.grp.levels, size = all.cells.levels)
txt_df <- data.frame(size = unique(df$size), 
                     stats = matrix(c(single_markers[use.genes[1], c('logFC', 'adj.p.value')],
                                               pseudo_markers2$`DESeq2-Wald`[use.genes[1], c('logFC', 'adj.p.value')],
                                               bulk_markers$`DESeq2-Wald`[use.genes[1], c('logFC', 'adj.p.value')],
                                               super_markers$`5`[use.genes[1], c('logFC', 'adj.p.value')],
                                               super_markers$`50`[use.genes[1], c('logFC', 'adj.p.value')],
                                               super_markers$`200`[use.genes[1], c('logFC', 'adj.p.value')]), byrow = T, ncol = 2),
                     vals = rep(5, length(unique(df$size))))
level_order <- c('single', 'super5', 'super50', 'super200', 'pseudo', 'bulk')
ggplot(data = df, aes(x = factor(size, level = level_order), y = values, fill = grp)) +
    geom_boxplot() +
    xlab('') +
    scale_y_log10() +
    ylab('log of logcounts') +
    ggtitle(use.genes[1]) +
    theme_classic() +
    annotate('text', x = factor(txt_df$size, level = level_order), y = txt_df$vals, label = paste0('logFC: ', round(unlist(txt_df$stats.1), 2), '\n', 'p value: ', format(txt_df$stats.2, scientific = T))) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))

# ---------------------------------------------------------
# Fraction
# ---------------------------------------------------------
frac <- fractionGenes(list(single = single_markers, 
                   super5 = super_markers$`5`, 
                   super10 = super_markers$`10`, 
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