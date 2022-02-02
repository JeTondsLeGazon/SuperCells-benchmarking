# 11th November 2021

# Script that contains all the analysis functions to compare methods


# Performs DE analysis with different methods passed as argument on a dataset
compute_DE_sc <- function(data,  # seurat object containing the normalized data
                       methods,  # list of methods to use
                       fc_thres)  # threshold below which no gene is selected
{
    tests_results <- list()
    for(method in methods){
        tests_results[[method]] <- FindAllMarkers(data, logfc.threshold = fc_thres, test.use = method) %>%
                                    arrange(p_val_adj, 1 / (abs(avg_log2FC) + 1)) %>%
                                    dplyr::rename(adj.p.value = p_val_adj, logFC = avg_log2FC)
    }
    return(tests_results)
}


# Performs DE analysis with different methods passed as argument on a bulk dataset
# Careful !! data should not be normalized at this point!!
compute_DE_bulk <- function(data, meta){
    
    process_table <- function(DE.table, log_fc_name, padj_name){
        
        DE.table %>%
            data.frame() %>%
            dplyr::rename(adj.p.value = padj_name, logFC = log_fc_name) %>%
            arrange(adj.p.value, 1 / (abs(logFC) + 1), decrease = T) %>%
            mutate(gene = row.names(.))
    }
    
    # DESeq2 approach
    dds <- DESeqDataSetFromMatrix(GetAssayData(data), 
                           colData = data@meta.data, 
                           design = ~ label)

    dds_wald <- DESeq(dds, test = 'Wald', minReplicatesForReplace = Inf)
    
    
    group.id <- grep('[Uu][Nn][Ss][Tt]|[Cc][Tt][Rr][Ll]', levels(data))
    
    results_wald <- results(dds_wald, contrast = c('label', levels(data)[c(2, 1)[group.id]], levels(data)[group.id])) %>%
                    process_table(log_fc_name = 'log2FoldChange', padj_name = 'padj')
    
    # edgeR approach
    edge <- DGEList(counts = GetAssayData(data), group = Idents(data))
    edge <- calcNormFactors(edge)
    model <- model.matrix(~Idents(data))
    edge <- estimateDisp(edge, model)
    
    # Quasi likelihood test
    fit <- glmQLFit(edge, model)
    qlf <- glmQLFTest(fit,coef= 2)$table %>%
        process_table(log_fc_name = 'logFC', padj_name = 'PValue')
    if(group.id == 2){
        qlf$logFC <- -(qlf$logFC)  # because does not know how to reverse log ...
    }
    
    
    return(list('DESeq2-Wald' = results_wald,
                'edgeR-QLF' = qlf))
    

}


# Computes differential expression for single-cell data
singleCell_DE <- function(sc_data, var.features, 
                          stat.test = 't', 
                          by.group = F){
    
    file.name <- 'singleCell.RData'

    sc_data <- FindVariableFeatures(sc_data, nfeatures = var.features)
    
    if(by.group){
        # Identify number of groups for paired test
        nb_groups <- sapply(levels(sc_data), 
                            function(x) as.numeric(str_sub(x, -1, -1)))
        
        single_markers <- c()
        for(i in seq_along(max(nb_groups))){
            group_id_treat <- grep(paste0('treat', i), levels(sc_data))
            group_id_ctrl <- grep(paste0('ctrl', i), levels(sc_data))
            markers <- FindMarkers(sc_data,
                                   ident.1 = levels(sc_data)[group_id_treat],
                                   ident.2 = levels(sc_data)[group_id_ctrl],
                                   only.pos = F, 
                                   logfc.threshold = 0, 
                                   test.use = stat.test) %>%
                mutate(gene = rownames(.))
            single_markers <- rbind(single_markers, markers)
        }
    }else{
        single_markers <- FindMarkers(sc_data,
                               ident.1 = levels(sc_clustered_data)[grep('treat', levels(sc_clustered_data))],
                               ident.2 = levels(sc_clustered_data)[grep('ctrl', levels(sc_clustered_data))],
                               only.pos = F, 
                               logfc.threshold = 0, 
                               test.use = stat.test) %>%
            mutate(gene = rownames(.))
    }
    
    single_markers <- single_markers %>%
        dplyr::rename(logFC = avg_log2FC, adj.p.value = p_val_adj) %>%
        arrange(adj.p.value, 1 / (abs(logFC) + 1), T) %>%
        subset(!duplicated(.$gene))
    
    return(single_markers)
}


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


# Computation of logFC and p-values for passed data with seurat or hyp test
# 
find_markers <- function(data, stat.test, seurat = F, norm = T){
    
    if(norm){
        data <- NormalizeData(data)
    }
    
    treat_grp <- grep('treat', Idents(data))
    ctrl_grp <- grep('ctrl', Idents(data))
    
    if(seurat){
        DE <- FindMarkers(data,
                    ident.1 = treat_grp,
                    ident.2 = ctrl_grp,
                    only.pos = F, 
                    logfc.threshold = 0, 
                    test.use = stat.test)
        return(arrangeDE(DE,
                       oldNameLog = 'avg_log2FC',
                       oldNameP = 'p_val_adj'))
    }
    #Log FC computation
    tmp <- expm1(GetAssayData(data)) %>% data.frame()
    grp1 <- rowMeans(tmp[, treat_grp])
    grp2 <- rowMeans(tmp[, ctrl_grp])
    logFCs <- log1p(grp1) - log1p(grp2)  # must do like this otherwise may / 0
    
    # t-test
    if(stat.test == 'wilcox'){
        hyp.test <- wilcox.test
    }else{
        hyp.test <- t.test
    }
    ge <- GetAssayData(data)
    pvals <- apply(ge, 1, function(x) hyp.test(x[treat_grp], x[ctrl_grp])$p.value)
    padj <- p.adjust(pvals, 'BH', nrow(ge))
    r <- data.frame(row.names = rownames(data),
                    p.value = pvals,
                    logFC = logFCs, 
                    adj.p.value = padj)
    return(arrangeDE(r))
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
        }else if(grepl('t', names(super_mc[i]))){
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
        title <- sprintf('True positive rate among the top 100 DE genes against Ground truth(Bulk)', tag)
    }else if (score.type == 'auc'){
        title <- sprintf('AUROC of DE genes against Ground truth(Bulk)', tag)
    }else if(score.type == 'tpr'){
        title <- sprintf('True positive rate of DE genes against Ground truth(Bulk)', tag)
    }
    if(score.type == 'tpr'){
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