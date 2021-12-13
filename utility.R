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

compute_aucc <- function(set1,
                         set2,  # second set of DE genes, ordered by p-values
                         k = 200)      # number of top interactions to perform
{
    results <- aucc(set1, set2, k = k)
    
    # s <- length(results[[2]])
    # plot(1:s, 
    #      results[[2]], 
    #      xlab = 'top n significant genes',
    #      ylab = 'Concordance score between sets', 
    #      main = sprintf('Concordance Curve (score AUCC: %s)', round(results[[1]], 3)), 
    #      type = 'l',
    #      lwd = 2,
    #      col = 'red',
    #      pch = 16,
    #      xlim = c(0, s),
    #      ylim = c(0, s))
    # lines(1:s, 1:s, lwd = 2)
    # legend(1, 100, legend=c("Concordance curve", "Optimal"), col=c("red", "black"), lty = 1, cex=0.8, lwd = 2)
    return(results[1])
    }


# Performs AUCC computing per cluster
compute_aucc_per_cluster <- function(DE1,
                                     DE2,
                                     clusters,
                                     k = 200)
{
    all_res <- list()
    for(cl in clusters){
        set1 <- subset(DE1, cluster == cl) %>%
            arrange(adj.p.value, 1 / (logFC + 1)) %>%
            select(gene) %>%
            unlist()

        set2 <- subset(DE2, cluster == cl) %>%
            arrange(adj.p.value, 1 / (logFC + 1)) %>%
            select(gene) %>%
            unlist()
        all_res[cl] <- list(aucc(set1, set2, k))
        s <- length(all_res[[cl]][[2]])
        plot(1:s, 
             all_res[[cl]][[2]], 
             xlab = 'top n significant genes',
             ylab = 'Concordance score between sets', 
             main = sprintf('Concordance Curve cluster %s (score AUCC: %s)', cl, round(all_res[[cl]][[1]], 3)), 
             type = 'l',
             lwd = 2,
             col = 'red',
             pch = 16,
             xlim = c(0, s),
             ylim = c(0, s))
        lines(1:s, 1:s, lwd = 2)
        legend(1, 100, legend=c("Concordance curve", "Optimal"), col=c("red", "black"), lty = 1, cex=0.8, lwd = 2)
    }
    all_res
}

# Save all figures from tmp_dir that are visualized in Rstudio into selected 
# folder
# TODO: create high level function that keeps in memory already copies files
savefig <- function(target_dir_path,  # directory in which figures will be saved
                    files_names = NULL)  # names of figures in appearing order
{
    # Check that target_dir_path exists
    if (!dir.exists(target_dir_path)){
        dir.create(target_dir_path)
    }
    
    plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
    plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
    
    
    file.copy(from = plots.png.paths, to = target_dir_path)
    
    # Check we have enough names for all files (more names will not be used)
    if (!is.null(files_names)){
        # get names of figures only
        splt_names <- unlist(strsplit(plots.png.paths, '/'))
        fig_names <- splt_names[grep(".png$", splt_names)]
        fig_names_clean <- fig_names[!fig_names %in% c('empty.png')]
        
        l <- length(fig_names_clean)
        if (length(files_names) < length(fig_names_clean)){
            message('Size of files_names should be greater or equal to the number of
                    files in tmp_dir')
        }else{
            # rename
            file.rename(file.path(target_dir_path, fig_names_clean), 
                        file.path(target_dir_path, files_names[1:l]))
        }
    }
    invisible(0)
}


# Cleans all .png from a folder
cleanfig <- function(target_dir_path)  # directory that contains the figures
{
    if (dir.exists(target_dir_path)){
        file.remove(list.files(target_dir_path, pattern = '.png$'))
    }else{
        message(sprintf('Could not find directory %s', target_dir_path))
    }
    invisible(0)
}


# Display significant genes from two groups to compare
# TODO: check with gather() for tidyverse of melt() from reshape
display_significant_genes <- function(seurat_obj,  # seurat object with data
                                      cluster_,  # cluster to check vs others
                                      super_markers,  # markers from super cell
                                      single_markers,  # markers from single cell
                                      topn = 10,  # number of top genes to compare
                                      samples = 500)  # subsampling level
{
    sub_sc <- NormalizeData(seurat_obj[, sample(colnames(seurat_obj), samples)])
    
    grouping <- factor(sub_sc$cell_type == cluster_, labels = c(cluster_, 'other'), levels = c(TRUE, FALSE))
    markers1 <- subset(super_markers, cluster == cluster_)$gene[1:topn]
    markers2 <- subset(single_markers, cluster == cluster_)$gene[1:topn]
    tot_markers <- c(markers1, markers2)
    
    ge <- as.matrix(GetAssayData(sub_sc)[tot_markers, ])
    to_plot <- data.frame(cell = NULL, gene = NULL, count = NULL)
    for(j in 1:ncol(ge)){
        to_plot <- rbind(to_plot, data.frame(cell = rep(colnames(ge)[j], nrow(ge)), 
                                             gene = rownames(ge), 
                                             count = ge[, j],
                                             cluster = rep(grouping[colnames(ge)[j]], nrow(ge)),
                                             super = tot_markers %in% markers1,
                                             single = tot_markers %in% markers2))
    }
    
    t <- to_plot[to_plot$super == T, ]
    p1 <- ggplot(data = t, aes(x = gene, y = count, color = cluster)) + 
        geom_boxplot() +
        scale_y_log10() +
        ggtitle(sprintf("Top %s Significant DE Genes for Super cells", topn)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ylab('Log10 normalized count')
    
    t <- to_plot[to_plot$single == T, ]
    p2 <- ggplot(data = t, aes(x = gene, y = count, color = cluster)) + 
        geom_boxplot() +
        scale_y_log10() +
        ggtitle(sprintf("Top %s Significant DE Genes for Single cells", topn)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        ylab('Log10 normalized count')
    
    print(p1 + p2)
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


plot_expression <- function(top_genes, seurat_data, title_append, sampling = F){
    p <- GetAssayData(seurat_data)[top_genes, ] %>%
        as.matrix() %>%
        melt() %>%
        mutate(ctrl = Idents(seurat_data)[Var2])
    if (sampling) {p <- sample(p, 100)}
    
    gplot <- ggplot(data = p, aes(x = Var1, y = value, color = ctrl)) +
        geom_point(size = 5) +
        scale_y_log10() +
        ylab('Log10count of genes') +
        xlab('Top DE genes according to GT') +
        ggtitle(sprintf('Comparison of top 10 DE genes from %s between groups using %s 
            gene expression matrix', title_append[1], title_append[2])) +
        theme(plot.title = element_text(hjust = 0.5))
    return(gplot)
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