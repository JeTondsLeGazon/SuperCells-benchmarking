# 3 February 2022

# # Contains all differential expression computation functions for 
# bulk / pseudo-bulk, single, metacell, supercell


# Compute DE via DESeq2 on a counts matrix with corresponding labels to columns 
# of counts matrix
computeDESeq2 <- function(counts, labels){
    
    message(sprintf('Computing DESeq2 on matrix of dimension %s', dim(counts)))
    if(!check_one_zero(counts)){  # mandatory for DESeq2 to work
        counts <-  counts + 1
    }
    dds <- DESeqDataSetFromMatrix(counts,
                                  colData = labels, 
                                  design = ~ label)
    dds_wald <- DESeq(dds, test = 'Wald', minReplicatesForReplace = Inf)
    results_wald <- results(dds_wald, contrast = c('label', 'treat', 'ctrl'))
    message('Done computing DESeq2')
    return(arrangeDE(results_wald, 
              oldNameLog = 'log2FoldChange',
              oldNameP = 'padj'))
}


# Compute DE via EdgeR on a counts matrix with corresponding labels to columns 
# of counts matrix
computeEdgeR <- function(counts, labels){
    
    message(sprintf('Computing EdgeR on matrix of dimension %s', dim(counts)))
    edge <- DGEList(counts = counts, group = labels)
    edge <- calcNormFactors(edge)
    model <- model.matrix(~labels)
    edge <- estimateDisp(edge, model)
    
    # Quasi likelihood test
    fit <- glmQLFit(edge, model)
    qlf <- glmQLFTest(fit,coef= 2)$table
    message('Done computing EdgeR')
    return(arrangeDE(qlf, 
                     oldNameLog = 'logFC',
                     oldNameP = 'PValue'))
}


# Wrapper for bulk DE computation using either DESeq2, EdgeR, or t-test
compute_DE_bulk <- function(data, labels, algo){
    
    counts <- data@assays$RNA@counts
    ge <- data@assays$RNA@data
    
    if(algo == 'DESeq2'){
        if(class(labels) != "data.frame"){
            labels <- data.frame(label = labels, row.names = colnames(data))
        }
        DE <- computeDESeq2(counts, labels)
    }else if(algo == 'EdgeR'){
        DE <- computeEdgeR(counts, labels)
    }else if(algo == 't-test'){
        DE <- find_markers(ge, labels)
    }else{
        stop('Cannot compute DE of expression, algorithm passed unknown')
    }
    return(DE)
}


# Wrapper for Single cell DE computation using either DESeq2, EdgeR, or t-test
compute_DE_single <- function(data, algo){
    
    labels <- data$label
    counts <- data@assays$RNA@counts
    ge <- data@assays$RNA@data
    
    if(algo == 'DESeq2'){
        DE <- computeDESeq2(counts, labels)
    }else if(algo == 'EdgeR'){
        DE <- computeEdgeR(counts, labels)
    }else if(algo == 't-test'){
        DE <- find_markers(ge, labels)
    }else{
        stop('Cannot compute DE of expression, algorithm passed unknown')
    }
    return(DE)
}


# Wrapper for Metacell DE computation using either DESeq2, EdgeR, or t-test
compute_DE_meta <- function(data, algo){
    
    ge <- log(data$ge)  # according to metacells vignette
    labels <- unlist(strsplit(data$sample, '[0-9]'))  # in case of split.by = sample
    counts <- floor(sweep(data$ge, 2, data$size, '*'))
    
    if(algo == 'DESeq2'){
        DE <- computeDESeq2(counts, labels)
    }else if(algo == 'EdgeR'){
        DE <- computeEdgeR(counts, labels)
    }else if(algo == 't-test'){
        DE <- find_markers(ge, labels)
    }else{
        stop('Cannot compute DE of expression, algorithm passed unknown')
    }
    return(DE)
}


# Wrapper for Metacell SC-like computation using either DESeq2, EdgeR, or t-test
compute_DE_metasc <- function(data, single_data, algo){
    
    # Metacell drops some cells along the way
    cells.use <- colnames(data)[colnames(data) %in% names(mc[[mc_gamma]]$membership)]
    
    # Input groups of Metacell into SuperCell GE computation
    ge <- supercell_GE(ge = single_data@assays$RNA@data[, cells.use],
                       groups =  data[[mc_gamma]]$membership[cells.use])
    
    # Use normalized counts (exm1(ge)) as input to compute normalized average counts
    counts <- supercell_GE(ge = expm1(single_data@assays$RNA@data[, cells.use]),
                           groups =  data[[mc_gamma]]$membership[cells.use])
    
    if(algo == 'DESeq2'){
        DE <- computeDESeq2(counts, labels)
    }else if(algo == 'EdgeR'){
        DE <- computeEdgeR(counts, labels)
    }else if(algo == 't-test'){
        DE <- find_markers(ge, labels)
    }else{
        stop('Cannot compute DE of expression, algorithm passed unknown')
    }
    return(DE)
}


# Wrapper for SuperCell DE computation using either DESeq2, EdgeR, or t-test
compute_supercell_DE <- function(SC, algo){
    labels <- SC$cell_line
    counts <- floor(sweep(SC$counts, 2, SC$supercell_size, '*'))
    if(algo == 'DESeq2'){
        DE <- computeDESeq2(counts, labels)
    }else if(algo == 'EdgeR'){
        DE <- computeEdgeR(counts, labels)
    }else if(algo == 't-test'){
        DE <- find_markers(SC$GE, labels)
    }else{
        stop('Cannot compute DE of expression, algorithm passed unknown')
    }
    return(DE)
}
