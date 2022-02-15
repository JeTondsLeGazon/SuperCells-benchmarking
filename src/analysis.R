# 11th November 2021

# Script that contains all the analysis functions to compare methods



# Compute concordance score between super cells DEA at different graining levels and another DEA
compute_score <- function(DEAs,  # list containing the results of DEA from supercells
                          GT,  # Ground truth for comparison with supercells
                          which.score = 'aucc')  # type of score to compute
{
    if (which.score == 'aucc'){
        max.k <- min(sapply(DEAs, nrow))
        if(max.k > 100){
            ks <- c(10, 50, 100)
        }else{
            ks <- c(round(max.k / 4, 0), round(max.k / 2, 0), max.k)
        }
            scores <- lapply(ks, function(k) sapply(DEAs, function(x) aucc(x$gene, GT$gene, k)))
            names(scores) <- ks
    }else if(which.score == 'tpr'){
        scores <- sapply(DEAs, function(x) tpr(x$gene, GT$gene))
    }
    return(scores)
}


# Computation of logFC and p-values for ge matrix (already normalized)
# Uses t-test as statistical test
find_markers <- function(ge, labels){
    
    message(sprintf('Computing t-test on matrix of dimension %s x %s', dim(ge)[1], dim(ge)[2]))
    treat_grp <- grep('treat', labels)
    ctrl_grp <- grep('ctrl', labels)
    
   
    #Log FC computation, computed on normalized counts instead of ge matrix
    tmp <- expm1(ge) %>% 
        data.frame()
    grp1 <- rowMeans(tmp[, treat_grp])
    grp2 <- rowMeans(tmp[, ctrl_grp])
    logFCs <- log2((grp1 + 0.001) / (grp2 + 0.001))  # must do like this otherwise may / 0
    
    pvals <- apply(ge, 1, function(x) t.test(x[treat_grp], x[ctrl_grp])$p.value)
    padj <- p.adjust(pvals, 'BH', nrow(ge))
    DE <- data.frame(row.names = rownames(ge),
                    p.value = pvals,
                    logFC = logFCs, 
                    adj.p.value = padj)
    return(arrangeDE(DE))
}


# Shows a LogFC x LogFC plots to compare values from two different methods
LogFcLogFcPlot <- function(stats1, stats2, title = ''){
    common.genes <- intersect(stats1$gene, stats2$gene)
    df <- data.frame(log1 = stats1[common.genes, 'logFC'],
                     log2 = stats2[common.genes, 'logFC'])
    ggscatter(df, x = 'log1', y = 'log2',
              add = 'reg.line', conf.int = T, cor.coef = T, cor.method = 'pearson',
              add.params = list(color = "blue", fill = "lightgray")) +
    geom_abline(color = 'red', size = 1)
}


# Plots for benchmarking purpose following special coloring and character plot
plot_results_BM <- function(super_mc, GT, GT.type, score.type = 'tpr_100'){
    plot(NULL, 
         ylim=c(0,1), 
         xlim=c(1, 100), 
         ylab="Scores", 
         xlab="Gammas", 
         log = 'x')
    
    # score to use
    score_func <- switch(score.type, 'tpr_100' = tpr_100, 'auc' = auc, 'tpr' = tpr)

    chr.used <- c()
    colors.used <- c()
    legends <- c()
    for(i in seq_along(super_mc)){
        matching <- unlist(lapply(super_mc[[i]], function(x) score_func(GT, x)))
        gammas <- names(matching)
        if(grepl('metasc', names(super_mc)[i])){
            tag <- 'MetaCells Super'
            chr.used <- c(chr.used, 8)
            colors.used <- c(colors.used, 'gold')
        }else if(grepl('meta', names(super_mc)[i])){
            tag <- 'MetaCells'
            chr.used <- c(chr.used, 4)
            colors.used <- c(colors.used, 'darkgreen')
        }else if(grepl('super', names(super_mc)[i])){
            tag <- 'SuperCells'
            chr.used <- c(chr.used, 1)
            colors.used <- c(colors.used, 'brown2')
        }else if(grepl('sub', names(super_mc)[i])){
            tag <- 'Subsampling'
            chr.used <- c(chr.used, 2)
            colors.used <- c(colors.used, 'black')
        }else{
            tag <- 'Random grouping'
            chr.used <- c(chr.used, 3)
            colors.used <- c(colors.used, 'grey')
        }
    
    
        if(grepl('t$', names(super_mc[i]))){
            testTag <- 't-test'
        }else if(grepl('twt', names(super_mc[i]))){
            testTag <- 'weighted t-test'
        }else if(grepl('[Dd][Ee][Ss]', names(super_mc[i]))){
            testTag <- 'DESeq2'
        }else{
            testTag <- 'EdgeR'
        }
        points(gammas, matching, col = colors.used[length(colors.used)], pch = chr.used[length(chr.used)], cex = 1.5)
        lines(gammas, matching, col = colors.used[length(colors.used)] , lwd = 1.5)
        super_mc_legend <- paste0(tag, ' (', testTag, ')')
        legends <- c(legends, 
                     paste(super_mc_legend, sprintf('Bulk (%s)', GT.type), sep = ' vs '))
    }
    if(score.type == 'tpr_100'){
        title <- sprintf('True positive rate among the top 100 DE genes against Ground truth (Bulk)', tag)
    }else if (score.type == 'auc'){
        title <- sprintf('AUROC of DE genes against Ground truth (Bulk)', tag)
    }else if(score.type == 'tpr'){
        title <- sprintf('True positive rate of DE genes against Ground truth (Bulk)', tag)
    }
    if(score.type == 'tpr' | score.type == 'auc'){
        legend.pos <- 'bottomleft'
    }else{
        legend.pos <- 'topright'
    }
    legend(legend.pos, 
           legend = legends, 
           col = colors.used, 
           pch = chr.used)
    title(title)
    grid()
}


# Calculate fraction of Genes belonging to different groups. Used to create figure
# in runAnalysis
fractionGenes <- function(DEs){
    lapply(DEs, function(DE) sapply(DE$gene, function(x) fractionGene(x, DEs)))
}

fractionGene <- function(gene, DEs){
    res <- sapply(DEs, function(x) gene %in% x$gene)
    paste(sort(names(DEs)[res]), collapse = '-')
}




# Shows rank of given genes between a set of markers and superCells like structure markers
rank_plot <- function(concerned_genes, test_markers, super_markers){
    
    rank1 <- match(concerned_genes, test_markers$gene)
    rank2 <- match(concerned_genes, super_markers$`1`$gene)
    rank3 <- match(concerned_genes, super_markers$`5`$gene)
    rank4 <- match(concerned_genes, super_markers$`50`$gene)
    rank5 <- match(concerned_genes, super_markers$`1000`$gene)
    
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
    p2 <- subplot(rank1, rank3, 5)
    p3 <- subplot(rank1, rank4, 50)
    p4 <- subplot(rank1, rank5, 1000)
    ggarrange(p1, p2, p3, p4)
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