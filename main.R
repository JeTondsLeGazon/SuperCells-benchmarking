# 11th November 2021

# Performs comparison between DEA from different sources, mainly SuperCells and
# ground truth (bulk DNA) from various dataset
setwd(paste0(getwd(), '/SuperCells-benchmarking'))

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
library(d)


source('utility.R')
source('supercells.R')
source('analysis.R')
source('processing.R')
#source('redefined_functions.R')

# Redefine wtd.t.test to get smaller p-values
#tmpfun <- get("wtd.t.test", envir = asNamespace("weights"))
#environment(custom.wtd.t.test) <- environment(tmpfun)
#attributes(custom.wtd.t.test) <- attributes(tmpfun)  # don't know if this is really needed
#assignInNamespace("wtd.t.test", custom.wtd.t.test, ns = "weights")

# ---------------------------------------------------------
# Data loading
# ---------------------------------------------------------
#available: Hagai2018_mouse-lps, Hagai2018_mouse-pic, Hagai2018_pig-lps, 
# Hagai2018_rabbit-lps, Hagai2018_rat-lps, Hagai2018_rat-pic, Angelidis2019_pneumo,
# Angelidis2019_alvmac, CanoGamez2020_naive-iTreg, CanoGamez2020_memory-Th17

scpath <- '../sc_rnaseq/rds'
bulkpath <- '../bulk_rnaseq/rds'
filename <- 'Hagai2018_mouse-lps.rds'

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
rename_sample <- function(data, new){
    old <- unique(data$sample)
    data$sample <- new[match(data$sample, old)]
    data
}


create_pseudobulk <- function(data){
    counts <- as.matrix(GetAssayData(data))
    sample.id <- sort(unique(data$sample))
    gene.names <- rownames(data)
    bulk <- lapply(sample.id, 
                   function(id) rowSums(counts[, which(data$sample == id)]))
    df <- data.frame(matrix(unlist(bulk), byrow = T, ncol = length(sample.id)), 
               row.names = gene.names)
    colnames(df) <- sample.id
    meta <- data.frame(row.names = sample.id, 
                       label = c('ctrl', 'treat')[as.numeric(sapply(sample.id, function(x) grepl('treat', x))) + 1],
                       replicate = rep(c(1, 2, 3), 2))
    sc <- CreateSeuratObject(df, meta.data = meta)
    Idents(sc) <- 'label'
    return(sc)
}

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
                             initial_centers = centers, 
                             labels = c('treat_grp1', 'ctrl_grp1'))


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

manual_bulk_markers <- find_markers_bulk(bulk_filtered_data) %>%
    arrange(adj.p.value, 1 / (abs(logFC) + 1)) %>%
    mutate(gene = row.names(.)) %>%
    subset(adj.p.value < 0.05 & logFC > 0) 


# ---------------------------------------------------------
# DE pseudobulk
# ---------------------------------------------------------
pseudo_markers <- pseudobulk_norm %>% 
                    FindVariableFeatures(nfeatures = 500) %>%
                    FindMarkers(ident.1 = 'ctrl',
                                ident.2 = 'treat',
                                only.pos = T, 
                                logfc.threshold = 0, 
                                test.use = 't') %>%
    mutate(gene = rownames(.)) %>%
    dplyr::rename(ajd.p.value = p_val_adj, logFC = avg_log2FC) %>%
    arrange(ajd.p.value, 1 / (abs(logFC) + 1), T)


# ---------------------------------------------------------
# DE pseudobulk
# ---------------------------------------------------------

# ---------------------------------------------------------
# DE supercells
# ---------------------------------------------------------

gammas <- c(1, 2, 5, 10, 50, 100, 200)
memory.limit(size=56000)

super_markers <- superCells_DEs(sc_clustered_data, gammas, 5)

volcano_plot(super_markers$`1`, logfc.thres = 0.5) +
    ggtitle('Volcano plot of SuperCells at level gamma = 5') +
    theme(plot.title = element_text(hjust = 0.5))

super_markers <- lapply(super_markers ,function(x) x %>% 
                            subset(adj.p.value < 0.05 & logFC > 0))


# ---------------------------------------------------------
# DE single cells
# ---------------------------------------------------------
# Single cell markers

single_markers <- singleCell_DE(sc_clustered_data, var.features = 500)
volcano_plot(single_markers, logfc.thres = 0.5) +
    ggtitle('Volcano plot of single cells from FindAllMarkers (seurat)') +
    theme(plot.title = element_text(hjust = 0.5))
single_markers <- single_markers %>%
    subset(adj.p.value < 0.05 & logFC > 0)


# ---------------------------------------------------------
# Comparison
# ---------------------------------------------------------
score_results <- compute_score(super_markers, manual_bulk_markers, which.score) %>%
    melt %>%
    mutate(gammas = rep(gammas, 3))
plot_score_results(score_results)

gt <- manual_bulk_markers
gt2 <- bulk_markers$`DESeq2-Wald`
match_scores <- data.frame(
    super_vs_single = unlist(lapply(super_markers, function(x) gene_match(x$gene, single_markers$gene))),
    super_vs_bulk = unlist(lapply(super_markers, function(x) gene_match(x$gene, gt$gene))),
    single_vs_bulk = c(gene_match(single_markers$gene, gt$gene), rep(NA, length(gammas)-1)),
    bulk_vs_bulk = c(gene_match(bulk_markers$`DESeq2-Wald`$gene, manual_bulk_markers$gene), rep(NA, length(gammas)-1)),
    super_vs_bulk2 = unlist(lapply(super_markers, function(x) gene_match(x$gene, gt2$gene))),
    
    gammas = gammas)

plot(gammas, match_scores$super_vs_single, type = 'b', 
     ylim = c(0, 1), col = 'blue', pch = 19, log = 'x')
lines(gammas, match_scores$super_vs_bulk, col = 'red')
points(gammas, match_scores$super_vs_bulk, col = 'red', pch = 19)

lines(gammas, match_scores$super_vs_bulk2, col = 'orange')
points(gammas, match_scores$super_vs_bulk2, col = 'orange', pch = 19)
points(gammas, match_scores$single_vs_bulk, col = 'green', pch=19)
points(gammas, match_scores$bulk_vs_bulk, col = 'black', pch=19)
legend('topright', legend = c('SuperCells (t-test) vs single cells (t-test)', 
                     'SuperCells (t-test) vs Bulk (t-test)',
                     'SuperCells (t-test) vs Bulk (DESeq2)',
                     'Single cells (t-test) vs Bulk (t-test)',
                     'Bulk (DESeq2) vs Bulk (t-test)'),
       col = c('blue', 'red', 'orange', 'green', 'black'), pch = 19)
grid()

# ---------------------------------------------------------
# Stats comparison
# ---------------------------------------------------------

log.thres <- 1
log.thres2 <- 1
gt_coverage_wrapper(subset(single_markers, logFC > log.thres), 
                    lapply(super_markers, function(x) subset(x, logFC > log.thres2)))

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

# Show log FC are different between bulk and single cells -> hence difference
p1 <- LogFcLogFcPlot(single_markers, bulk_markers$`edgeR-QLF`) + 
    xlab('LogFC single cells (seurat)') +
    ylab('LogFC bulk (DESeq2 Wald)')
p2 <- LogFcLogFcPlot(bulk_markers$`edgeR-QLF`, bulk_markers$`DESeq2-Wald`) +
    xlab('LogFC bulk (edgeR QLF)') +
    ylab('LogFC bulk (DESeq2 Wald)')
p3 <- LogFcLogFcPlot(bulk_markers$`DESeq2-Wald`, super_markers$`1`) + 
    xlab('LogFC bulk (DESeq2 Wald)') +
    ylab('LogFC Super Cells gamma = 1')
p4 <- LogFcLogFcPlot(bulk_markers$`DESeq2-Wald`, super_markers$`5`) + 
    xlab('LogFC bulk (DESeq2 Wald)') +
    ylab('LogFC Super Cells gamma = 5')

fig <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
annotate_figure(fig, top = text_grob('LogFC vs LogFC graph for Bulk and single cells', 
                                     face = 'bold', size = 14))

# ---------------------------------------------------------
# Analysis pvalue and logFC
# ---------------------------------------------------------
selected_genes <- super_markers$`1`$gene[1:1000]
logFCs <- data.frame(lapply(super_markers, function(x) x[selected_genes, 'logFC']), row.names = selected_genes, check.rows = F)
colnames(logFCs) <- gammas
logFCs$diff1_50 <- logFCs[1] - logFCs[5]

#Observe difference between gammas in terms of logFCs -> select those
target_genes <- rownames(logFCs)[logFCs$diff1_50/logFCs[1] > 0.5]

tmp <- logFCs[target_genes, -6] %>% group_by(gene = rownames(.)) %>% melt()
ggplot(data = tmp, aes(x = variable, y = value, group = gene, color = gene)) + 
    geom_line() + 
    geom_point()

# Observe counts and normalized counts of genes with drops in LogFCS
df <- as.matrix(GetAssayData(sc_filtered_data)[target_genes, ]) %>% melt() %>% data.frame() %>% mutate(grp = Idents(sc_filtered_data)[Var2])
p1 <- ggplot(data = df, aes(x = value, y = Var1, color = grp)) + 
    geom_boxplot()
df <- as.matrix(GetAssayData(sc_filtered_data)[selected_genes[1:3], ]) %>% melt() %>% data.frame() %>% mutate(grp = Idents(sc_filtered_data)[Var2])
p2 <- ggplot(data = df, aes(x = value, y = Var1, color = grp)) + 
    geom_boxplot()
df <- as.matrix(GetAssayData(sc_normalized_data)[target_genes, ]) %>% melt() %>% data.frame() %>% mutate(grp = Idents(sc_filtered_data)[Var2])
p3 <- ggplot(data = df, aes(x = value, y = Var1, color = grp)) + 
    geom_boxplot()
df <- as.matrix(GetAssayData(sc_normalized_data)[selected_genes[1:3], ]) %>% melt() %>% data.frame() %>% mutate(grp = Idents(sc_filtered_data)[Var2])
p4 <- ggplot(data = df, aes(x = value, y = Var1, color = grp)) + 
    geom_boxplot()


fig <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
annotate_figure(fig, top = text_grob('Counts and normalized gene expression for dropped LogFC genes and top 3 genes', 
                                     face = 'bold', size = 14))

# ---------------------------------------------------------
# Rank top n genes
# ---------------------------------------------------------
concerned_genes <- single_markers$gene[1:100]
rank_plot(concerned_genes, bulk_markers$`DESeq2-Wald`, super_markers)

common.genes <- intersect(super_markers$`1`$gene, super_markers$`50`$gene)
df <- data.frame(p1 = super_markers$`1`[common.genes, 'adj.p.value'],
                 p50 = super_markers$`50`[common.genes, 'adj.p.value'],
                 std1 = super_markers$`1`[common.genes, 'std.err'],
                 std50 = super_markers$`50`[common.genes, 'std.err'],
                 df1 = super_markers$`1`[common.genes, 'df'],
                 df50 = super_markers$`50`[common.genes, 'df'],
                 t1 = super_markers$`1`[common.genes, 't.value'],
                 t50 = super_markers$`50`[common.genes, 't.value'],
                 top100 = common.genes %in% concerned_genes,
                 m1.1 = super_markers$`1`[common.genes, 'w.mean.1'],
                 m1.50 = super_markers$`50`[common.genes, 'w.mean.1'],
                 m2.1 = super_markers$`1`[common.genes, 'w.mean.2'],
                 m2.50 = super_markers$`50`[common.genes, 'w.mean.2'])
p1 <- ggplot(data = df, aes(x = -log10(p1), y = -log10(p50), color = top100)) +
        geom_point() +
        geom_smooth() +
        ylab('-log10 p gamma = 50') +
        xlab('-log10 p gamma = 1')
p2 <- ggplot(data = df, aes(x = std1, y = std50, color = top100)) +
    geom_point() +
    geom_smooth() +
    ylab('std err gamma = 50') +
    xlab('std err gamma = 1')
p3 <- ggplot(data = df, aes(x = as.numeric(df1), y = as.numeric(df50), color = top100)) +
    geom_point() +
    geom_smooth() +
    ylab('df gamma = 50') +
    xlab('df gamma = 1')
p4 <- ggplot(data = df, aes(x = as.numeric(t1), y = as.numeric(t50), color = top100)) +
    geom_point() +
    geom_smooth() +
    ylab('t gamma = 50') +
    xlab('t gamma = 1')
p5 <- ggplot(data = df, aes(x = as.numeric(m1.1), y = as.numeric(m1.50), color = top100)) +
    geom_point() +
    geom_smooth() +
    ylab('mean 1 gamma = 50') +
    xlab('mean 1 gamma = 1')
p6 <- ggplot(data = df, aes(x = as.numeric(m2.1), y = as.numeric(m2.50), color = top100)) +
    geom_point() +
    geom_smooth() +
    ylab('mean 2 gamma = 50') +
    xlab('mean 2 gamma = 1')
ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)

# ---------------------------------------------------------
# Weighted vs unweighted
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
