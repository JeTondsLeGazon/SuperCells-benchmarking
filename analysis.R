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
    #dds_lrt <- DESeq(dds, test = 'LRT', minReplicatesForReplace = Inf, 
    #                 reduced = ~1)
    dds_wald <- DESeq(dds, test = 'Wald', minReplicatesForReplace = Inf)
    
    #results_lrt <- results(dds_lrt) %>%
    #                process_table(log_fc_name = 'log2FoldChange', padj_name = 'padj')
    
    results_wald <- results(dds_wald, contrast = c('label', levels(data)[1], levels(data)[2])) %>%
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
    qlf$logFC <- -(qlf$logFC)  # because does not know how to reverse log ...
    
    
    return(list('DESeq2-Wald' = results_wald,
                'edgeR-QLF' = qlf))
    

}


# Computes differential expression for single-cell data
singleCell_DE <- function(sc_data, var.features){
    sc_data <- FindVariableFeatures(sc_data, nfeatures = var.features)
    
    # Identify number of groups for paired test
    nb_groups <- sapply(levels(sc_data), 
                        function(x) as.numeric(str_sub(x, -1, -1)))
    
    single_markers <- c()
    for(i in seq_along(max(nb_groups))){
        group_id_treat <- grep(paste0('^treat.+', i, '$'), levels(sc_data))
        group_id_ctrl <- grep(paste0('^ctrl.+', i, '$'), levels(sc_data))
        markers <- FindMarkers(sc_data,
                               ident.1 = levels(sc_data)[group_id_treat],
                               ident.2 = levels(sc_data)[group_id_ctrl],
                               only.pos = F, 
                               logfc.threshold = 0) %>%
            mutate(gene = rownames(.))
        single_markers <- rbind(single_markers, markers)
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
    max.k <- min(sapply(DEAs, nrow))
    if(max.k > 100){
        ks <- c(10, 50, 100)
    }else{
        ks <- c(round(max.k / 4, 0), round(max.k / 2, 0), max.k)
    }
        scores <- lapply(ks, function(k) sapply(DEAs, function(x) aucc(x$gene, GT$gene, k)))
        names(scores) <- ks
    return(scores)
}


# Manual computation of logFC and t-test for bulk data
find_markers_bulk <- function(bulkData){
    
    bulkData <- NormalizeData(bulkData)
    
    nb_samples <- dim(bulkData)[2]
    #Log FC computation
    tmp <- exp(GetAssayData(bulkData)) %>% data.frame()
    grp1 <- rowMeans(tmp[, 1:nb_samples/2])
    grp2 <- rowMeans(tmp[, (nb_samples/2 + 1):nb_samples])
    logFCs <- log1p(grp1 + 1) - log1p(grp2 + 1)
    
    # t-test
    ge <- GetAssayData(bulkData)
    dfs <- lapply(1:nrow(ge), function(i) data.frame(value = ge[i, ], 
                                                     grp = c(rep('lsp', nb_samples/2), rep('ctrl', nb_samples/2))))
    pvals <- sapply(dfs, function(x) t.test(x$value ~ x$grp)$p.value)
    padj <- p.adjust(pvals, 'BH')
    r <- data.frame(row.names = rownames(bulkData), 
                    logFC = logFCs, 
                    adj.p.value = padj)
    return(r)
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