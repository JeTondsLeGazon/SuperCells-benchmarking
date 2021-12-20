# 8th November 2021

# Utility functions for the supercells DE analysis benchmark

library(testit)


# Computes the Area under the Concordance curve between two sets of genes
aucc <- function(set1,  # first set of DE genes, ordered by p-values
                 set2,  # second set of DE genes, ordered by p-values
                 k)      # number of top interactions to perform
{
    # Check variables
    assert('K should be a positive non-null integer', k >= 0 && k%%1 == 0)
    assert('Size of compared sets of DE genes should be greater or equal than k', 
           length(set1) >= k)
    assert('Size of compared sets of DE genes should be greater or equal than k', 
           length(set2) >= k)
    
    res <- vector()
    for (i in 1:k){
        # here the order is not considered important
        res <- append(res, as.integer(sum(!is.na(match(set1[1:i], set2[1:i])))))
    }
    max_s <- k * (k+1) / 2
    return(sum(res)/max_s)
}


# Computes true positive rate (TPR) between ground truth and a candidate
tpr <- function(gt, other){
    sum(!is.na(match(gt, other))) / length(other)
}



plot_score_results <- function(scores){
   df <- data.frame(scores)
   labels <- sapply(unique(df$L1), function(x) paste0('Top n = ', x))
   ggplot(data = df, aes(x = as.factor(gammas), 
                         y = value, 
                         label = as.factor(gammas),
                         fill = log(gammas))) +
       geom_bar(stat = 'identity', alpha = 0.5) +
       facet_wrap(~as.numeric(L1), labeller = labeller(labels)) +
       xlab('Graining level') +
       ylab('AUCC score') +
       ylim(c(0, 1)) +
       ggtitle(sprintf('AUCC score between superCell and Ground truth for top %s, %s and %s Genes', unique(df$L1)[1], unique(df$L1)[2], unique(df$L1)[3])) +
       theme(legend.position = "none",
             axis.title.y = element_text(size=14, face="bold"),
             axis.title.x = element_text(size=14, face="bold"),
             plot.title = element_text(hjust = 0.5, size = 15, face = 'bold'))
}



gt_coverage <- function(de_gt, de_other, sig_level = 0.05){
    s_gt <- de_gt[de_gt$adj.p.value < sig_level, 'gene']
    s_other <- de_other[de_other$adj.p.value < sig_level, 'gene']
    sum(!is.na(match(s_gt, s_other))) / length(s_gt)
}

gt_coverage_wrapper <- function(de_gt, des_other){
    
    sig_levels <- c(0.05)
    coverages <- lapply(sig_levels, function(y) lapply(des_other, function(x) gt_coverage(de_gt, x, y)))
    sizes_gts <- lapply(sig_levels, function(x) nrow(de_gt[de_gt$adj.p.value < x, ]))
    sizes_other <- lapply(sig_levels, function (y) lapply(des_other, function(x) nrow(x[x$adj.p.value < y, ])))
    
    data <- data.frame(gammas = as.numeric(rep(names(des_other), length(sig_levels))), 
                       sig_level = as.factor(sort(rep(sig_levels, length(names(des_other))), decreasing = T)), 
                       values = unlist(coverages) * 100,
                       Ngt = sort(rep(unlist(sizes_gts), length(names(des_other))), decreasing = T),
                       Nsuper = unlist(sizes_other))
    
    ggplot(data = data, aes(x = gammas, y = values)) +
        geom_line() +
        geom_point(size = 2) + 
        geom_text(aes(label = paste0('N gt = ', Ngt, '\n', 'N super = ', Nsuper)), 
                  nudge_y = -1, size = 3) +
        scale_x_continuous(trans = 'log2') +
        xlab('Graining level') + 
        ylab('Percentage of GT DE genes recovery') +
        ggtitle('Percentage of DE GT genes recovery from superCells at different 
                graining levels and different significance level') +
        theme(plot.title = element_text(hjust = 0.5))
}

scaling_factor <- function(seuratdata, method = 'manual'){
    if(method == 'manual'){
        logcounts <- log(as.matrix(GetAssayData(seuratdata)) + 1)
        geom_avg <- rowMeans(logcounts)
        
        scaling_log_factor <- apply(logcounts - geom_avg, 2, median)
        return(exp(scaling_log_factor))
    }else{
        dds <- DESeqDataSetFromMatrix(GetAssayData(seuratdata), 
                                      colData = seuratdata@meta.data, 
                                      design = ~ label)
        dds <- estimateSizeFactors(dds)
        return(sizeFactors(dds))
    }
}


# Computes match between two sets, proportionally to the size of the first
compute_match <- function(set1, set2){
    sum(!is.na(match(set1, set2))) / length(set2)
}


# Volcano plot to show DE genes according to p value and logFC
volcano_plot <- function(stats, 
                         logfc.thres = 1, 
                         thres.annotated.fc = NULL, 
                         thres.annotated.p = NULL,
                         top.genes = 15){
    stats <- stats[!is.na(stats$adj.p.value), ]
    
    # In case of p value to 0, put them at the minimum
    if(min(stats$adj.p.value) == 0){
        stats$adj.p.value[stats$adj.p.value == 0] <- min(stats$adj.p.value[stats$adj.p.value != 0]) / 10
    }
    

    
     p <- ggplot(data = stats, aes(x = logFC, y = -log10(adj.p.value), color = adj.p.value < 0.05 & abs(logFC) > logfc.thres, label = gene)) +
            geom_point() +
            geom_vline(xintercept = logfc.thres, linetype = 'dashed') +
            geom_vline(xintercept = -logfc.thres, linetype = 'dashed') +
            geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
            theme(legend.position = 'none')
     
     # Use top.genes instead of thresholds
    if (is.null(thres.annotated.fc) & is.null(thres.annotated.p)){
        sub.stats <- stats %>%
            subset(logFC > 0) %>%
            slice_head(n = top.genes)
        thres.annotated.fc <- min(sub.stats$logFC)
        thres.annotated.p <- max(sub.stats$adj.p.value)
    }

    p + geom_text_repel(aes(label=ifelse(logFC >= thres.annotated.fc & adj.p.value <= thres.annotated.p & logFC > logfc.thres, as.character(gene),'')), 
                            box.padding = 0.5, color = 'black', max.overlaps = Inf)
}


# Shows rank of given genes between a set of markers and superCells markers
rank_plot <- function(concerned_genes, test_markers, super_markers){

    rank1 <- match(concerned_genes, test_markers$gene)
    rank2 <- match(concerned_genes, super_markers$`1`$gene)
    rank3 <- match(concerned_genes, super_markers$`2`$gene)
    rank4 <- match(concerned_genes, super_markers$`5`$gene)
    rank5 <- match(concerned_genes, super_markers$`10`$gene)
    
    subplot <- function(x, y, g){
        qplot(x, 
              y, 
              xlab = 'Single-cell (Seurat)', 
              ylab = 'Super-cells') +
            geom_abline(intercept = 0, slope = 1, color = 'red', size = 1) +
            ggtitle(paste0('Top 100 DE genes comparison at gamma = ', g)) +
            xlim(c(0, 100)) +
            ylim(c(0, 100))
    }
    p1 <- subplot(rank1, rank2, 1)
    p2 <- subplot(rank1, rank3, 2)
    p3 <- subplot(rank1, rank4, 5)
    p4 <- subplot(rank1, rank5, 10)
    ggarrange(p1, p2, p3, p4)
}


# Check percentage of matching top n DE genes between supercells and other DEs
plot_matches <- function(supercell_res, others, legends = NULL){
    gammas <- as.numeric(names(supercell_res))
    plot(NULL, ylim=c(0,1), xlim=c(min(gammas), max(gammas)), 
         ylab="Match scores", xlab="Gammas", log = 'x')
    chr <- c(8, 15, 16, 17, 18, 4, 3)
    for (i in seq_along(others)){
        matching <- unlist(lapply(supercell_res, function(x) gene_match(x$gene, others[[i]]$gene)))
        points(gammas, matching, col = i, pch = chr[i], cex = 1.5)
        lines(gammas, matching, col = i, lwd = 1.5)
    }
    legend('topright', 
           legend = legends, col = seq_along(others), pch = chr[seq_along(others)])
    grid()

    scores <- c()
    names <- c()
    for (i in seq_along(others)){
        for(j in seq_len(length(others) - i)){
            matching <- gene_match(others[[i]]$gene, others[[j + i]]$gene)
            scores <- c(scores, matching)
            names <- c(names, paste0(names(others)[i], '_vs_', names(others)[j + i]))
        }
    }
    p <- barplot(scores, las=2, ylim = c(0,1))
    text(x = p, y = scores + 0.17, labels = names, srt = 90)
    print(p)
}


# Compare statistics from t-test at different graining levels
compare_statistics <- function(super_markers){
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
}
