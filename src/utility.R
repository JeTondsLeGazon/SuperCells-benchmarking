# 8th November 2021

# Utility functions for the supercells DE analysis benchmark

library(testit)

# Creates the sample according to the condition (treatment vs control) and the
# replicate if available (mouse/patient 1, mouse/patient 2, ...)
createSample <- function(data){
    if('replicate' %in% names(data[[]])){
        if(nchar(unique(data$replicate)[1]) > 1){
            spt <- str_split(data$replicate, '')
            pos <- sapply(spt, function(x) grep('[0-9]', x))
            samp <- sapply(seq_along(pos), function(i) spt[[i]][pos[i]])
            return(paste0(data$label, samp))
        }else{
            return(paste0(data$label, data$replicate))
        }
        
    }else if('sample' %in% names(data[[]])){
        spt <- str_split(data$sample, '_')
        pos <- sapply(spt, function(x) grep('mouse|patient|subject', x))
        samples <- sapply(seq_along(pos), function(i) spt[[i]][pos[i]])
        
        spt <- str_split(samples, '')
        pos <- sapply(spt, function(x) grep('[0-9]', x))
        samp <- sapply(seq_along(pos), function(i) spt[[i]][pos[i]])
        return(paste0(data$label, samp))
    }else{
        return(0)
    }
}


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


# Computes area under the roc curve from two statistics (gt and other)
auc <- function(gt, other){
    thresholds <- c(1-10, 1-5, 1-3, 0.01, 0.05, 0.1, 0.5, 1)
    pos <- gt %>% subset(adj.p.value < 0.05)
    neg <- gt %>% subset(adj.p.value >= 0.05)
    points <- c(1, 0)
    for(thres in thresholds){
        tp <- sum(!is.na(match(pos$gene, subset(other, other$adj.p.value < thres)$gene)))
        tn <- sum(!is.na(match(neg$gene, subset(other, other$adj.p.value >= thres)$gene)))
        fp <- sum(!is.na(match(neg$gene, subset(other, other$adj.p.value < thres)$gene)))
        fn <- sum(!is.na(match(pos$gene, subset(other, other$adj.p.value >= thres)$gene)))
        sensi <- tp / (tp + fn)
        speci <- tn / (tn + fp)
        points <- c(points, c(speci, sensi))
    }
    points <- data.frame(matrix(points, ncol = 2, byrow = T))
    colnames(points) <- c('specificity', 'sensitivity')
    auc <- sum(diff((1 - points$specificity))*rollmean(points$sensitivity,2))
    return(auc)
}


tpr <- function(gt, other){
    tp <- sum(!is.na(match(subset(gt, gt$adj.p.value < 0.05)$gene, 
                     subset(other, other$adj.p.value < 0.05)$gene)))
    return(tp / nrow(subset(gt, gt$adj.p.value < 0.05)))
    
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


# Compute scores of DE genes correlated between superCells and other markers
plot_results <- function(supercell_res, others, super.type, score.type = 'match'){
    gammas <- as.numeric(names(supercell_res))
    plot(NULL, 
         ylim=c(0,1), 
         xlim=c(min(gammas), max(gammas)), 
         ylab="Scores", 
         xlab="Gammas", 
         log = 'x')
    
    # supercell legend
    super_legend <- sprintf('SuperCells (%s)', super.type)
    
    # score to use
    score_func <- switch(score.type, 'match' = tpr, 'auc' = auc, 'tpr' = tpr)
    if(score.type == 'match'){
        supercell_res <- lapply(supercell_res, function(x) x[1:100, ])
        others <- lapply(others, function(x) x[1:100, ])
    }
    
    if(score.type == 'match'){
        title <- 'True positive rate among the top 100 DE genes between superCell and Ground truth'
    }else if (score.type == 'auc'){
        title <- 'AUROC of DE genes between superCell and Ground truth'
    }else if(score.type == 'tpr'){
        title <- 'True positive rate of DE genes between superCell and Ground truth'
    }
    
    chr <- c(8, 15, 16, 17, 18, 4, 3)
    legends <- c()
    for (i in seq_along(others)){
        matching <- unlist(lapply(supercell_res, function(x) score_func(others[[i]], x)))
        points(gammas, matching, col = i, pch = chr[i], cex = 1.5)
        lines(gammas, matching, col = i, lwd = 1.5)
        legends <- c(legends, paste(super_legend, names(others)[i], sep = ' vs '))
    }
    legend('bottomleft', 
           legend = legends, 
           col = seq_along(others), 
           pch = chr[seq_along(others)])
    title(title)
    grid()
    
    # score for singlets
    if(F){
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


fractionGenes <- function(DEs){
    lapply(DEs, function(DE) sapply(DE$gene, function(x) fractionGene(x, DEs)))
}


fractionGene <- function(gene, DEs){
    res <- sapply(DEs, function(x) gene %in% x$gene)
    paste(sort(names(DEs)[res]), collapse = '-')
}


# Create Gene Expression matrix (GE) of metacells from saved files mcGExxx.rds
createMCGE <- function(sc_data, results_folder){
    if(!any(grepl('mcGE*', list.files(results_folder)))){
        stop("Cannot find metacells gene expression matrices ('mcGExxx.rds')")
    }
    
    mcGEs <- list()
    samples <- unique(sc_data$sample)
    for(sample in samples){
        mcGEs[[sample]] <- readRDS(file.path(results_folder, sprintf("mcGE%s.rds", sample)))
    }
    mcGEs <- lapply(mcGEs, function(x) x[!unlist(lapply(x, is.null))])
    mcGEs <- lapply(seq_along(mcGEs[[1]]), function(i) lapply(mcGEs, function(x) x[[i]]))
    names(mcGEs) <- as.character(c(1, 10, 20, 30, 50))  # TODO: add to configs
    finalGEs <- list()
    for(GE in mcGEs){
        genesToKeep <- Reduce(intersect, sapply(GE, rownames))
        tmp <- lapply(GE, function(x) x[genesToKeep, ])
        resultGE <- data.frame(tmp[[1]], row.names = genesToKeep)
        coln <- sprintf(paste0(substr(samples[1], 1, nchar(samples[1])),'_%s'), seq_len(ncol(tmp[[1]])))
        for(i in 2:length(tmp)){
            resultGE <- cbind(resultGE, tmp[[i]])
            coln <- c(coln, 
                      sprintf(paste0(substr(samples[i], 1, nchar(samples[i])),'_%s'), seq_len(ncol(tmp[[i]]))))
        }
        colnames(resultGE) <- coln
        gamma <- floor(ncol(sc_data) / ncol(resultGE))
        finalGEs[[as.character(gamma)]] <- resultGE
    }
    return(finalGEs)
}


# Create membership of metacells to be fed into supercells from saved files 
# mcCompositionxxx.rds
createMCMembership <- function(sc_data, results_folder){
    if(!any(grepl('mcComposition*', list.files(results_folder)))){
        stop("Cannot find metacells memberships files ('mcCompositionxxx.rds')")
    }
    
    samples <- unique(sc_data$sample)
    mcCompositions <- list()
    for(sample in samples){
        mcCompositions[[sample]] <- readRDS(file.path(results_folder, sprintf("mcComposition%s.rds", sample)))
    }
    mcCompositions <- lapply(mcCompositions, function(x) x[!unlist(lapply(x, is.null))])
    mcCompositions <- lapply(seq_along(mcCompositions[[1]]), function(i) lapply(mcCompositions, function(x) x[[i]]))
    names(mcCompositions) <- as.character(c(1, 10, 20, 30, 50))  # TODO: add to configs
    
    memberships <- list()
    annotations <- list()
    for(mc in mcCompositions){
        curMax <- 0
        membership <- c()
        annotation <- c()
        for(sample in names(mc)){
            membership <- c(membership, mc[[sample]] + curMax)
            annotation <- c(annotation, rep(sample, length(mc[[sample]])))
            curMax <- curMax + max(mc[[sample]])
        }
        gamma <- floor(ncol(sc_data) / length(unique(membership)))
        memberships[[as.character(gamma)]] <- membership
        annotations[[as.character(gamma)]] <- annotation
    }
    return(list(membership = memberships,
                annotation = annotations))
}


# Arrange DE tables by creating a gene columns, arranging in decreasing order
# according to p values and logFC, renaming columns and subseting
arrangeDE <- function(DE, oldNameLog = NULL, oldNameP = NULL, subset_logFC = T){
        DE <- DE %>%
            data.frame() %>%
            dplyr::rename(logFC = oldNameLog, adj.p.value = oldNameP) %>% 
            mutate(gene = row.names(.)) %>%
            arrange(adj.p.value, 1/(abs(logFC) + 1))
    if(subset_logFC){
        subset(DE, logFC > 0)
    }else{
        DE
    }
}


# Create random groups between cells and compute their mean or sum
randomGrouping <- function(data, gamma, operation = 'mean'){
    availableCols <- seq(1:ncol(data))
    grps <- list()
    grpNum <- 1
    if(operation == 'mean'){
        opFunc <- rowMeans
    }else{
        opFunc <- rowSums
    }
    while(length(availableCols) > gamma){
        subSample <- sample(availableCols, gamma)
        availableCols <- setdiff(availableCols, subSample)
        grps[[grpNum]] <- opFunc(as.matrix(data[, subSample]))
        grpNum <- grpNum + 1
    }
    grps[[grpNum]] <- opFunc(as.matrix(data[, availableCols]))
    res <- data.frame(grps)
    colnames(res) <- seq(1:grpNum)
    return(res)
}


# Save markers in appropriate folder with suitable name
saveMarkers <- function(markers,
                        algo,
                        split.by,
                        base.path,
                        kind){
    if((!algo %in% c('DESeq2', 'EdgeR', 't-test')) | 
       (!kind %in% c('super', 'subsampling', 'random', 'bulk', 'meta', 'metasc'))){
        stop('Could not save markers')
    }

    stat.method <- switch(algo, 'DESeq2' = 'des', 'EdgeR' = 'edge', 't-test' = 't')
    filename <- paste0(paste(kind, stat.method, sep = '_'), '.rds')
    if(kind != 'bulk'){
        top.dir <- paste('markers', split.by, sep = '_')
    }else{
        top.dir <- 'GT'
    }
    full.dir <- file.path(base.path, top.dir)
    if(!dir.exists(full.dir)){
        dir.create(full.dir, recursive = T)
    }
    saveRDS(markers, file.path(full.dir, filename))
}