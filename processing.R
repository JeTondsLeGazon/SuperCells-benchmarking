# 11th November 2021
# Contains all the functions necessary for the processing pipeline

library(scDblFinder)


# Computes scaling factor for RNA composition and library size normalization
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


# Quality control for the single-cell data prior to filtering and downstream
# analysis
singleCell_qc <- function(sc_data){
    
    # Outliers check
    p <- FeatureScatter(sc_data, 
                        feature1 = 'nFeature_RNA', 
                        feature2 = 'nCount_RNA',
                        group.by = 'label') +
        ggtitle('Scatter plot of counts vs number of different genes per cells')
    print(p)
    
    # Check quality with mitochondrial genes
    sc_data <- PercentageFeatureSet(sc_data, "^MT-", col.name = "percent_mito")
    
    # Check quality with ribosomal protein (could also indicate dead cells)
    sc_data <- PercentageFeatureSet(sc_data, "^RP[SL]", col.name = "percent_ribo")
    
    # Check quality with hemoglobin (blood contamination)
    sc_data <- PercentageFeatureSet(sc_data, "^HB[^(P)]", col.name = "percent_hb")
    
    # Counts per cell, should be > 500 to be usable
    p <- tidyseurat::ggplot(data = sc_data, aes(x = nCount_RNA, color = label, fill = label)) +
        geom_density(alpha = 0.2) +
        scale_x_log10() +
        theme_classic() + 
        ylab('Cell density') +
        geom_vline(xintercept = 500, color = 'red', size = 1) +
        ggtitle('Distribution of counts per cell') +
        theme(plot.title = element_text(size = 20, face = "bold"))
    print(p)
    
    # Gene per cell
    p <- tidyseurat::ggplot(data = sc_data, aes(x = nFeature_RNA, color = label, fill = label)) + 
        geom_density(alpha = 0.2) +
        scale_x_log10() +
        theme_classic() + 
        ylab('Cell density') +
        geom_vline(xintercept = 300, color = 'red', size = 1) +
        ggtitle('Distribution of number of different genes per cell') +
        theme(plot.title = element_text(size = 15, face = "bold"))
    print(p)
    
    # Doublet density
    sc_data$doubletScore <- computeDoubletDensity(GetAssayData(sc_data))
    p <- tidyseurat::ggplot(data = sc_data, aes(x = doubletScore, color = label, fill = label)) +
        geom_density(alpha = 0.2) +
        scale_x_log10() +
        theme_classic() +
        ylab('Cell density') +
        geom_vline(xintercept = quantile(sc_data$doubletScore, 0.95), color = 'red', size = 1) +
        ggtitle('Detection of doublet through scDblFinder') +
        theme(plot.title = element_text(size = 20, face = "bold"))
   print(p)
    
    p <- VlnPlot(sc_data, features = c('percent_mito', 'percent_ribo', 'percent_hb'),
            group.by = 'label') +
            #ggtitle('Violin plots of mitochondrial, ribosomal and hemoglobin genes fraction of total genes') +
        theme(plot.title = element_text(size = 10, face = "bold"))
    print(p)
    
    # Transcript capture efficiency
    df <- data.frame(x = log2(rowSums(GetAssayData(sc_data)) + 1),
                     y = rowSums(GetAssayData(sc_data) > 0) / ncol(GetAssayData(sc_data)))
    p <- ggplot(data = df, aes(x = x, y = y)) +
        geom_point(color='darkred') +
        geom_smooth() +
        ylab('Fraction of cells with detected gene') +
        xlab('Total logcounts per gene') +
        ggtitle('Capture efficiency') +
        theme(plot.title = element_text(size = 20, face = "bold")) +
        geom_vline(xintercept = log2(200), color = 'red', size = 1)
    
    p <- ggMarginal(p, margins = 'x', 
                    color="darkred", 
                    fill = 'darkred',
                    alpha = 0.2)
    print(p)
    return (sc_data)
}


# Filtering of single cells according to parameter found in singleCell_qc
singleCell_filtering <- function(sc_data, 
                                 max.doublet.percentile = 0.95,
                                 min.gene.per.cell = 300,
                                 min.count.per.cell = 500,
                                 min.count.per.genes = 200,
                                 max.ribo.percent = 40,
                                 max.mito.percent = 20,
                                 max.hb.percent = 5){
    
    cat('Dimensions prior to filtering: ')
    cat(paste0(dim(sc_data)[1], ' genes x ', dim(sc_data)[2], ' cells\n'))
    genes.use <- rownames(GetAssayData(sc_data))[rowSums(GetAssayData(sc_data)) > min.count.per.genes]
    sc_data <- subset(sc_data, 
                      subset = (percent_mito < max.mito.percent) &
                          (percent_ribo < max.ribo.percent) &
                          (percent_hb < max.hb.percent) &
                          (nFeature_RNA > min.gene.per.cell) &
                          (nCount_RNA > min.count.per.cell) & 
                          (doubletScore < quantile(doubletScore, max.doublet.percentile)),
                      features = genes.use)
    cat('Dimensions after filtering: ')
    cat(paste0(dim(sc_data)[1], ' genes x ', dim(sc_data)[2], ' cells\n'))
    return(sc_data)
}


# Normalization for single cell data based on either seurat normalization or manual
# log-norm from DESeq2
NormalizeObject <- function(data, method = 'DESeq2', scaling.method = 'manual'){
    normalized_data <- data
    if (method == 'DESeq2'){
        norm_factors <- scaling_factor(data, 'else')
        normalized_data[['RNA']] <- CreateAssayObject(as.matrix(log(sweep(GetAssayData(data), 2, norm_factors, '/') + 1)))
        
    }else{
        normalized_data <- NormalizeData(data)
    }
    return(normalized_data)
}



# Quality control and processing of bulk RNA data but without normalization due
# to the use of DESeq2
bulk_qc_and_filtering <- function(dataset){
    
    # Detection based filtering
    transcript_level <- 50
    dropped_genes <- rowSums(GetAssayData(dataset)) < transcript_level
    cat(sprintf('Dropped %s genes after detection based filtering\n', sum(dropped_genes)))
    dataset <- subset(dataset, features = rownames(dataset)[!dropped_genes])
    return(dataset)
}



# Filter cells based on doublet scores
doubletFinder <- function(data,  # seurat object containing the data
                          threshold = 2.5)  # threshold to be considered doublet
{
    scores <- computeDoubletDensity(GetAssayData(data))
    return(colnames(data)[scores > threshold])
    
}


# Clustering into subgroups because of immune cells sub-population
sub_cluster <- function(singleCell_data){
    sc <- singleCell_data %>%
            FindVariableFeatures(selection.method = "vst", nfeatures = 500)
    
    plot1 <- VariableFeaturePlot(sc)
    plot2 <- LabelPoints(plot = plot1, 
                         points = head(VariableFeatures(sc), 20), 
                         repel = TRUE)
    print(plot2)
    
    sc <- sc %>%
        ScaleData() %>%
        RunPCA() %>%
        RunUMAP(dims = 1:10)
    
    p <- DimPlot(sc, reduction = 'umap') +
        ggtitle('Umap dimension reduction for cell idents')
    print(p)
    
    if('replicate' %in% names(sc@meta.data)){
        p <- DimPlot(sc, reduction = 'umap', group.by = 'replicate') +
            ggtitle('Umap dimension reduction between replicates')
        print(p)
    }
    
    batch.id <- grep('batch', names(sc@meta.data))
    if(length(batch.id) > 0){
        p <- DimPlot(sc, reduction = 'umap', group.by = names(sc@meta.data)[batch.id]) +
            ggtitle('Umap dimension reduction between batches')
        print(p)
    }
    return(sc)
}


# Change identification of cells based on clustering / analysis
reIdent <- function(sc, initial_centers = NULL, labels  = NULL){
    if(is.null(initial_centers) & is.null(labels)){
        return(sc)
    }
    if(is.null(initial_centers)){
        new.labels <- labels
        names(new.labels) <- unique(Idents(sc))
        sc <- RenameIdents(sc, new.labels)
        print(DimPlot(sc, reduction = 'umap'))
        return(sc)
    }
    clustering <- kmeans(Embeddings(sc, reduction = 'umap'), 
                         centers = initial_centers, 
                         iter.max = 10, nstart = 1)
    sc$label <- clustering$cluster
    
    sorted_sizes <- clustering$withinss
    new.idents <- c()
    for(i in seq_along(labels)){
        new.idents <- c(new.idents, 
                        as.character(which(clustering$withinss == sorted_sizes[i])))
    }
    new.idents.labels <- labels

    names(new.idents.labels) <- new.idents
    sc$label <- new.idents.labels[sc$label]
    sc <- subset(sc, subset = label %in% labels)
    Idents(sc) <- 'label'
    print(DimPlot(sc, reduction = 'umap'))
    return(sc)
}


# Rename samples according to new vector, careful as no check is done
rename_sample <- function(data, new){
    old <- unique(data$sample)
    data$sample <- new[match(data$sample, old)]
    data
}


# Create pseudobulk from seurat data, expect samples to have treat in their name
create_pseudobulk <- function(data){
    counts <- as.matrix(GetAssayData(data))
    sample.id <- sort(unique(data$sample))
    gene.names <- rownames(data)
    bulk <- lapply(sample.id, 
                   function(id) rowSums(counts[, which(data$sample == id)]))
    df <- data.frame(matrix(unlist(bulk), byrow = F, ncol = length(sample.id)), 
                     row.names = gene.names)
    colnames(df) <- sample.id
    meta <- data.frame(row.names = sample.id, 
                       label = c('ctrl', 'treat')[as.numeric(sapply(sample.id, function(x) grepl('treat', x))) + 1],
                       replicate = rep(c(seq_len(length(sample.id)/2)), 2))
    sc <- CreateSeuratObject(df, meta.data = meta)
    Idents(sc) <- 'label'
    return(sc)
}