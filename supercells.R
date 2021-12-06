#9th November 2021

# Wrapper for the Supercells functionnality to perform QC and DE analysis
library(SuperCell)

#Computes DE genes per cluster and in total for Supercells
superCells_DE <- function(ge_matrix,  # gene expression matrix counts
                          gamma,  # graining level of super cells
                          knn,  # knn used for supercells construction
                          clusters_  # clusters of each cell
)
{
    super <- SCimplify(ge_matrix,
                       k.knn = knn,
                       gamma = gamma,
                       n.var.genes = 1000,
                       directed = FALSE
    )
    
    super$cell_line <- supercell_assign(clusters = clusters_,
                                        supercell_membership = super$membership,
                                        method = "jaccard")
    
    super$GE <- supercell_GE(ge_matrix, super$membership)
    
    
    
    # DE
    markers_super_per_cluster <- supercell_FindAllMarkers(ge = super$GE,
                                                          supercell_size = super$supercell_size, 
                                                          clusters = super$cell_line,  # different from what proposed
                                                          logfc.threshold = 0,
                                                          only.pos = T,
                                                          do.bootstrapping = F)
    
    # create single data.frame() instead of split list for each cluster
    total_super_markers <- data.frame(markers_super_per_cluster[[1]][FALSE, ])
    for(i in seq_along(markers_super_per_cluster)){
        if (class(markers_super_per_cluster[[i]]) == 'data.frame'){
            total_super_markers <- rbind(total_super_markers, markers_super_per_cluster[[i]] %>%
                                             data.frame() %>%
                                             mutate(cluster = names(markers_super_per_cluster)[i],
                                                    gene = rownames(.)))
        }
    }
    
    #filename <- paste('DE_results/superCells_markers_', gamma, '.csv', sep = '')
    #write.csv(total_super_markers, file = filename, row.names = T)
    
    
    total <- total_super_markers <- total_super_markers %>%
        arrange(adj.p.value, 1 / (abs(logFC) + 1))
    
    return(total)
}


# Uses new cell.annotation to perform DE for each cluster according to treatment
# vs control and merge all DEs together for a single gamma
superCells_DE_by_cluster <- function(data,  # gene expression matrix counts
                                     gamma  # graining level of super cells
)
{
    
    super <- SCimplify(GetAssayData(data),
                       cell.annotation = Idents(data),
                       k.knn = 5,
                       gamma = gamma,
                       n.var.genes = 1000,
                       directed = FALSE
    )
        
    super$cell_line <- supercell_assign(clusters = Idents(data),
                                        supercell_membership = super$membership,
                                        method = "jaccard")
        
    super$GE <- supercell_GE(GetAssayData(data), super$membership)
    
    # --------------------------------------------------------------------------
    
    # DE per cluster group
    clusters <- unique(super$cell_line)
    nb_groups <- sapply(clusters, function(x) as.character(str_sub(x, -1, -1)))
    
    DEs <- c()
    for(i in seq_along(max(nb_groups)))
    {
        DE <- supercell_FindMarkers(ge = super$GE,
                                    supercell_size = super$supercell_size,
                                    clusters = super$cell_line,
                                    ident.1 = clusters[grep(paste0('^treat.+', i, '$'), clusters)],
                                    ident.2 = clusters[grep(paste0('^ctrl.+', i, '$'), clusters)],
                                    logfc.threshold = 0,
                                    only.pos = T,
                                    do.bootstrapping = F)
        
        DE <- DE %>%
                mutate(gene = rownames(.))
        DEs <- rbind(DEs, DE)
    }
    DEs <- DEs %>%
        arrange(adj.p.value, 1 / (abs(logFC) + 1), T) %>%
        subset(!duplicated(.)) %>%
        set_rownames(.$gene)
    return(DEs)
}


# Computes DE genes for superCells at different graining levels
superCells_DEs <- function(data,  # normalized logcounts seurat object
                           gammas,  # list of graning levels to use
                           knn,  # number of nearest neighbors to use
                           method = 'standard') 
{
    super_DEs <- list()
    for(gam in gammas){
        if(method == 'standard'){
            super_res <- superCells_DE(GetAssayData(data),
                                       gamma = gam,
                                       knn = knn,
                                       Idents(data))
        }else{
            super_res <- superCells_DE_by_cluster(data,
                                       gamma = gam)
        }
        super_DEs[[as.character(gam)]] <- super_res
    }
    return(super_DEs)
}