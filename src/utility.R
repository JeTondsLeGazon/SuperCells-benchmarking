# 8th November 2021

# Utility functions for the supercells DE analysis benchmark
 
library(pROC)


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


# Uses function auc from pROC to compute the area under the curve between two
# sets of markers
auc <- function(gt, other, logfc.thresh = 0){
    positives <- subset(gt, adj.p.value < 0.05 & logFC > logfc.thresh)$gene
    labels <- rep(F, nrow(other))
    labels[which(other$gene %in% positives)] <- T
    result <- pROC::auc(labels, other$adj.p.value)
    return(as.numeric(result))
}


# True positive rate (TPR) between ground truth (GT) and another set of markers
# Takes subset of statistically significant markers from ground truth and 
# compares it to the same number of markers in other set
tpr <- function(gt, other, logfc.thresh = 0){
    gt.markers <- subset(gt, adj.p.value < 0.05 & logFC > logfc.thresh)$gene
    N <- length(gt.markers)
    other.markers <- other$gene[1:N]
    tp <- sum(other.markers %in% gt.markers)
    return(tp / N)
}


# True positive rate (TPR) between ground truth and another set of markers for
# the first 100 markers in each set
tpr_100 <- function(gt, other, logfc.thres = NULL){
    N <- 100
    gt.markers <- gt$gene[1:N]
    other.markers <- other$gene[1:N]
    tp <- sum(other.markers %in% gt.markers)
    return(tp / N)
}


# Computes the scaling factor for normalization for a set of data either
# manually or via DESeq2 for normalization over library size and RNA population composition
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


# Arrange DE tables by creating a gene columns, arranging in decreasing order
# according to p values and logFC, renaming columns and subseting
arrangeDE <- function(DE, oldNameLog = NULL, oldNameP = NULL, subset_logFC = F){
        DE <- DE %>%
            data.frame() %>%
            dplyr::rename(logFC = oldNameLog, adj.p.value = oldNameP) %>% 
            mutate(gene = row.names(.)) %>%
            arrange(adj.p.value, 1/(abs(logFC) + 1))
        return(DE)
}


# Save markers in appropriate folder with suitable name
saveMarkers <- function(markers,
                        algo,
                        split.by,
                        base.path,
                        kind){
    cell.grouping <- c('super', 'subsampling', 'random', 'bulk', 'pseudo',
                       'meta', 'metasc', 'single')
    algo.types <- c('DESeq2', 'EdgeR', 't-test')
    if((!algo %in% algo.types) | 
       (!kind %in% cell.grouping)){
        stop('Could not save markers')
    }

    stat.method <- switch(algo, 'DESeq2' = 'des', 'EdgeR' = 'edge', 't-test' = 't')
    filename <- paste0(paste(kind, stat.method, sep = '_'), '.rds')
    if(kind == 'single'){
        top.dir <- 'single'
    }else if(kind != 'bulk'){
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


# Check if all rows contain at least one zero
# Used for DESeq2 as algorithm does not work if all rows contains some zeros
check_one_zero <- function(ge){
    row_contains_zero <- apply(ge, 1, function(row) sum(row == 0) == 0)
    return(sum(row_contains_zero) > 0)
}


# Load markers from different algorithm into a list
load_markers <- function(markers.type, algos, split.by, results.path){
    algo.pairing <- list('DESeq2' = 'des',
                         'EdgeR' = 'edge',
                         't-test' = 't')
    if(sum(algos %in% names(algo.pairing)) != length(algos)){
        stop('Could not find some algorithms in provided list')
    }
    
    markers <- list()
    for(algo in algos){
        stat.method <- algo.pairing[[algo]]
        filename <- paste0(paste(markers.type, stat.method, sep = '_'), '.rds')
        if(markers.type == 'single'){
            top.dir <- 'single'
        }else if(markers.type != 'bulk'){
            top.dir <- paste('markers', split.by, sep = '_')
        }else{
            top.dir <- 'GT'
        }
        full.path <- file.path(results.path, top.dir, filename)
        if(!file.exists(full.path)){
            stop(sprintf('Could not find markers at %s', full.path))
        }
        markers[[algo]] <- readRDS(full.path)
    }
    return(markers)
}
