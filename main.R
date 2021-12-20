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


source('utility.R')
source('supercells.R')
source('analysis.R')
source('processing.R')


# ---------------------------------------------------------
# Meta parameters
# ---------------------------------------------------------
filename <- 'Hagai2018_mouse-lps.rds'

stat.test <- 'wilcox'
weighted <- F
resetSingle <- F  # saved with wilcox
resetSuper <- T


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

rownames(bulk_raw$meta) <- bulk_raw$meta$sample
bulk_data <- CreateSeuratObject(bulk_raw$assay, meta.data = bulk_raw$meta)

# if CanoGamez, troubles in the meta ...
Idents(bulk_data) <- bulk_raw$meta$label
bulk_data$label <- Idents(bulk_data)

Idents(sc_data) <- 'label'
set.seed(0)
# ---------------------------------------------------------
# QC and filtering
# ---------------------------------------------------------
# Effects of normalization:
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
sc_data <- singleCell_qc(sc_data)
sc_filtered_data <- singleCell_filtering(sc_data, 
                                         max.ribo.percent = 55)#,
                                         #min.gene.per.cell = 300,
                                         #min.count.per.cell = 400)
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
sc_normalized_data <- NormalizeObject(sc_filtered_data, method = 'else')


# ---------------------------------------------------------
# Pseudo bulk creation
# ---------------------------------------------------------
pseudobulk_data <- create_pseudobulk(rename_sample(sc_filtered_data, 
                                                   c('treat1', 'ctrl1', 'treat2', 'ctrl2', 'treat3', 'ctrl3')))
pseudobulk_norm <- NormalizeObject(pseudobulk_data, method = 'else')


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
centers <- matrix(c(5, -5, 5, -5, 5, 0, -10, -10), ncol = 2)
sc_clustered_data <- reIdent(sc_clustered_data, 
                             initial_centers = NULL, 
                             labels = c('treat', 'ctrl'),
                             replicate = T)


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
# DE bulk
# ---------------------------------------------------------

bulk_markers <- compute_DE_bulk(bulk_filtered_data)
volcano_plot(bulk_markers$`edgeR-QLF`) +
    ggtitle('Volcano plot of bulk data from DESeq2 (wald)') +
    theme(plot.title = element_text(hjust = 0.5))

bulk_markers <- lapply(bulk_markers, function(x) subset(x, adj.p.value < 0.05 & logFC > 0))
bulk_markers <- lapply(bulk_markers, function(x) if(nrow(x) == 0){x <- NULL}else{x})
bulk_markers <- bulk_markers[!unlist(lapply(bulk_markers, is.null))]


# ---------------------------------------------------------
# DE bulk manual
# ---------------------------------------------------------

manual_bulk_markers <- find_markers_bulk(bulk_filtered_data, stat.test) %>%
    arrange(adj.p.value, 1 / (abs(logFC) + 1)) %>%
    mutate(gene = row.names(.)) %>%
    subset(adj.p.value < 0.05 & logFC > 0) 


# ---------------------------------------------------------
# DE pseudobulk Seurat
# ---------------------------------------------------------
pseudo_markers <- pseudobulk_data %>% 
                    FindVariableFeatures(nfeatures = 500) %>%
                    FindMarkers(ident.1 = 'treat',
                                ident.2 = 'ctrl',
                                only.pos = T, 
                                logfc.threshold = 0, 
                                test.use = stat.test) %>%
    mutate(gene = rownames(.)) %>%
    dplyr::rename(adj.p.value = p_val_adj, logFC = avg_log2FC) %>%
    arrange(adj.p.value, 1 / (abs(logFC) + 1), T) %>%
    subset(logFC > 0 & adj.p.value < 0.05)


# ---------------------------------------------------------
# DE pseudobulk DESeq2
# ---------------------------------------------------------
pseudo_markers2 <- compute_DE_bulk(pseudobulk_data)
volcano_plot(pseudo_markers2$`edgeR-QLF`) +
    ggtitle('Volcano plot of pseudo bulk data from DESeq2 (wald)') +
    theme(plot.title = element_text(hjust = 0.5))

pseudo_markers2 <- lapply(pseudo_markers2, function(x) subset(x, adj.p.value < 0.05 & logFC > 0))
pseudo_markers2 <- lapply(pseudo_markers2, function(x) if(nrow(x) == 0){x <- NULL}else{x})
pseudo_markers2 <- pseudo_markers2[!unlist(lapply(pseudo_markers2, is.null))]


# ---------------------------------------------------------
# DE pseudobulk manual
# ---------------------------------------------------------
half.id <- ncol(pseudobulk_norm) / 2
if(stat.test == 'wilcox'){
    hyp.test <- wilcox.test
}else{
    hyp.test <- t.test
}
pseudo_markers_manual <- lapply(1:nrow(pseudobulk_norm), 
                  function(i) unlist(hyp.test(x = GetAssayData(pseudobulk_norm)[i, (half.id + 1) : ncol(pseudobulk_norm)],
                              y = GetAssayData(pseudobulk_norm)[i, 1: half.id])))
if(stat.test == 'wilcox'){
    pseudo_markers_manual <- data.frame(matrix(unlist(pseudo_markers_manual), ncol = 6, byrow = T)[, c(1, 2)], 
                                        row.names = rownames(pseudobulk_norm))
    colnames(pseudo_markers_manual) <- c('w', 'p.value')
}else{
    pseudo_markers_manual <- data.frame(matrix(unlist(pseudo_markers_manual), ncol = 12, byrow = T)[, c(1, 2, 3, 9)], 
                                        row.names = rownames(pseudobulk_norm))
    colnames(pseudo_markers_manual) <- c('t', 'df', 'p.value', 'std.err')
}

pseudo_markers_manual$adj.p.value <- p.adjust(pseudo_markers_manual$p.value, method = 'BH')
pseudo_markers_manual$logFC <- log(rowSums(expm1(GetAssayData(pseudobulk_norm)[ , (half.id + 1) : ncol(pseudobulk_norm)])) / rowSums(expm1(GetAssayData(pseudobulk_norm)[ , 1 : half.id])))
pseudo_markers_manual <- pseudo_markers_manual %>%
    mutate(gene = rownames(.)) %>%
    arrange(adj.p.value, 1 / (abs(logFC) + 1)) %>%
    subset(logFC > 0 & adj.p.value < 0.05)


# ---------------------------------------------------------
# DE supercells
# ---------------------------------------------------------
gammax <- ncol(sc_clustered_data) / length(unique(Idents(sc_clustered_data))) / 3
if(gammax > 1000){
    gammas <- c(1, 2, 5, 10, 50, 100, 200, 1000)
}else{
    gammas <- c(1, 2, 5, 10, 50, 100, 200)
}
memory.limit(size=56000)
super_markers <- superCells_DEs(sc_clustered_data, gammas, 5, 
                                resetData = resetSuper,
                                weighted = weighted,
                                test.use = stat.test)

volcano_plot(super_markers$`1`, logfc.thres = 0.5) +
    ggtitle('Volcano plot of SuperCells at level gamma = 5') +
    theme(plot.title = element_text(hjust = 0.5))

super_markers <- lapply(super_markers ,function(x) x %>% 
                            subset(adj.p.value < 0.05 & logFC > 0))


# ---------------------------------------------------------
# DE single cells
# ---------------------------------------------------------
# Single cell markers

single_markers <- singleCell_DE(sc_clustered_data, var.features = 500, 
                                resetData = resetSingle, stat.test)
volcano_plot(single_markers, logfc.thres = 0.5) +
    ggtitle('Volcano plot of single cells from FindAllMarkers (seurat)') +
    theme(plot.title = element_text(hjust = 0.5))
single_markers <- single_markers %>%
    subset(adj.p.value < 0.05 & logFC > 0)


# ---------------------------------------------------------
# Comparison
# ---------------------------------------------------------
score_results <- compute_score(super_markers, manual_bulk_markers, 'aucc') %>%
    melt %>%
    mutate(gammas = rep(gammas[seq_along(super_markers)], 3))
plot_score_results(score_results)
score_results <- compute_score(super_markers, manual_bulk_markers, 'tpr')


plot_matches(super_markers, 
             list(single = single_markers, 
               bulk_man = manual_bulk_markers, 
               bulk_des = bulk_markers$`DESeq2-Wald`,
               pseudo_des = pseudo_markers2$`DESeq2-Wald`,
               pseudo_man = pseudo_markers_manual ),
             c(sprintf('SuperCells (%s-test) vs single cells (%s-test seurat)', stat.test, stat.test), 
               sprintf('SuperCells (%s-test) vs Bulk (manual %s-test)', stat.test, stat.test),
               sprintf('SuperCells (%s-test) vs Bulk (DESeq2)', stat.test),
               sprintf('SuperCells (%s-test) vs pseudo-bulk (DESeq2)', stat.test),
               sprintf('SuperCells (%s-test) vs pseudo-bulk (manual %s test)', stat.test, stat.test)))

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

           