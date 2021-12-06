# 11th November 2021
# Contains all the functions necessary for the processing pipeline

library(scDblFinder)


# Quality control and processing of single cell RNA data
sc_qc_and_filtering <- function(dataset, plot_figures = F){
    
    # Remove doublets
    doublet_cells <- doubletFinder(dataset)
    
    cat(sprintf('Found %s doublets\n', length(doublet_cells)))
    sc <- subset(dataset, cells = WhichCells(dataset, colnames(dataset)[!(colnames(dataset) %in% doublet_cells)]))
    
    if (plot_figures){
        # check potential outliers for counts
        FeatureScatter(sc, feature1 = 'nFeature_RNA', feature2 = 'nCount_RNA',
                       group.by = 'cell_type')
    }
    
    
    # Check quality with mitochondrial genes
    sc <- PercentageFeatureSet(sc, "^MT-", col.name = "percent_mito")
    
    # Check quality with ribosomal protein (could also indicate dead cells)
    sc <- PercentageFeatureSet(sc, "^RP[SL]", col.name = "percent_ribo")
    
    # Check quality with hemoglobin (blood contamination)
    sc <- PercentageFeatureSet(sc, "^HB[^(P)]", col.name = "percent_hb")
    
    
    if (plot_figures){
        VlnPlot(sc, features = c('nCount_RNA', 'nFeature_RNA' ,'percent_mito', 'percent_ribo', 'percent_hb'),
                group.by = 'cell_type')
        
        # Transcript capture efficiency
        smoothScatter(log2(rowSums(GetAssayData(sc))+1), 
                      rowSums(GetAssayData(sc)>0)/ncol(GetAssayData(sc)),
                      pch=19,col="red", ylab = 'Detection probability',
                      xlab = 'Total cell count')
        
        # Library size
        hist(colSums(GetAssayData(sc))/1e3, xlab = 'Library size in thousands', 
             ylab = 'Number of cells', breaks = 10, main = 'Library size')
    }
    
    # Detection based filtering
    gene_level <- 500  # number of minimum genes per cell -> drop cells
    cell_level <- 30  # number of minimum cell with gene -> drop genes
    transcript_level <- 50  # number of minimum transcript level for a gene -> drop gene
    
    dim.prefiltering <- dim(sc)
    filtered_sc <- subset(sc, cells = WhichCells(sc, expression = nFeature_RNA > gene_level))
    filtered_sc <- subset(filtered_sc, 
                          features = rownames(filtered_sc)[rowSums(GetAssayData(filtered_sc) > 0) > cell_level])
    filtered_sc <- subset(filtered_sc,
                          features = rownames(filtered_sc)[rowSums(GetAssayData(filtered_sc)) > transcript_level])
    
    cat(sprintf('Dropped %s genes and %s cells by detection-based filtering\n', 
                nrow(sc)-nrow(filtered_sc), ncol(sc) - ncol(filtered_sc)))
    
    sc <- filtered_sc
    
    if (plot_figures){
        # Top genes
        n = 20
        frac_per_cell <- GetAssayData(sc) / colSums(GetAssayData(sc)) * 100
        most_expressed <- order(apply(frac_per_cell, 1, median), decreasing = T)[n:1]
        boxplot(t(as.matrix(frac_per_cell[most_expressed, ])), horizontal = TRUE, 
                xlab = '% of total count per cell', cex = 0.1, las = 1, 
                col = (scales::hue_pal())(n)[n:1])
    }
    # Mito / ribo filtering
    # Based on observation, remove cells with too high mitochondrial gene fraction 
    # or too low ribosomic gene fraction (proposed 20% and 5%, respectively)
    cat(sprintf('Dropped %s cells based on mitochondrial percentage\n', sum(sc$percent_mito >= 20)))

    
    # Normalization and scaling
    sc <- subset(sc, cells = WhichCells(sc, expression = percent_mito < 20))
    
    return(sc)
}

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



# Normalization for single cell data based on either seurat normalization or manual
# log-norm from DESeq2
NormalizeObject <- function(data, method = 'DESeq2'){
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
    
    print(DimPlot(sc, reduction = 'umap'))
    # CanoGamez2020_memory-Th17: matrix(c(-10, 10, -2, 2), ncol = 2)
    # mouse-lps : matrix(c(-5, -5, 7, 7, -5, 10, 10, 0), ncol = 2)
    # pig-lps: matrix(c(-5, 5, 0, 0), ncol = 2)
    # rat-lps: matrix(c(-5, -5, 7, 7, -5, 5, 5, -5), ncol = 2)
    initial_centers <- matrix(c(-5, -5, 7, 7, -5, 10, 10, 0), ncol = 2)
    clustering <- kmeans(Embeddings(sc, reduction = 'umap'), 
                         centers = initial_centers, 
                         iter.max = 10, nstart = 1)
    Idents(sc) <- clustering$cluster
    sorted_sizes <- sort(clustering$withinss, T)
    new.idents <- c(as.character(which(clustering$withinss == sorted_sizes[1])),
                    as.character(which(clustering$withinss == sorted_sizes[2])),
                    as.character(which(clustering$withinss == sorted_sizes[3])),
                    as.character(which(clustering$withinss == sorted_sizes[4])))
    new.idents.labels <- c('treat_grp1', 'ctrl_grp1', 'ctrl_grp2', 'treat_grp2')
    # Hagai rabbit: c('lsp4_grp1', 'ctrl_grp1', 'ctrl_grp2', 'lps4_grp2')
    # CanoGamaz: c('Th17', 'unstim')
    # Hagai lps: c('ctrl_grp1', 'lps4_grp1', 'ctrl_grp2', 'lps4_grp2')
    names(new.idents.labels) <- new.idents
    sc <- RenameIdents(sc, new.idents.labels)
    print(DimPlot(sc, reduction = 'umap'))
    Idents(singleCell_data) <- Idents(sc)
    return(singleCell_data)
}