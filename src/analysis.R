# 11th November 2021

# Script that contains all the analysis functions to compare methods



# Computes the match between the topn genes between two sets, normalized [0, 1]
gene_match <- function(set1, set2, topn = 100){
    
    if(topn > min(length(set1), length(set2))){
        topn <- min(length(set1), length(set2))
    }
    
    sum(!is.na(match(set1[1:topn], set2[1:topn]))) / topn
}


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
    logFCs <- log1p(grp1) - log1p(grp2)  # must do like this otherwise may / 0
    
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


# plots different super/meta cells vs others
plot_results_flex <- function(super_mc, others, score.type = 'match'){
    plot(NULL, 
         ylim=c(0,1), 
         xlim=c(1, 100), 
         ylab="Scores", 
         xlab="Gammas", 
         log = 'x')
    
    # score to use
    score_func <- switch(score.type, 'match' = tpr, 'auc' = auc, 'tpr' = tpr)
    if(score.type == 'match'){
        super_mc <- lapply(super_mc, function(x) lapply(x, function(xx) xx[1:100, ]))
        others <- lapply(others, function(x) x[1:100, ])
    }
    
    # Case when we want to compare each element of super_mc to a single in others
    if(length(others) == 1){
        others <- rep(others, length(super_mc), simpifly = F)
    }

    chr.used <- c()
    colors.used <- c()
    legends <- c()
    for (i in seq_along(others)){
        matching <- unlist(lapply(super_mc[[i]], function(x) score_func(others[[i]], x)))
        gammas <- as.numeric(names(super_mc[[i]]))
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
                     paste(super_mc_legend, names(others)[i], sep = ' vs '))
    }
    if(score.type == 'match'){
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