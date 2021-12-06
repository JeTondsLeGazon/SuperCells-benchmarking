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
    
    
    # QC
    # supercell_plot(super$graph.supercells, 
    #                group = super$cell_line, 
    #                seed = 1, 
    #                main = "Super-cell colored by cell line assignment")
    # 
    # purity <- supercell_purity(clusters = cell_types, 
    #                            supercell_membership = super$membership)
    # hist(purity, main = "Super-cell purity \nin terms of cell line composition")
    
    
    # Clustering
    # super.PCA <- supercell_prcomp(Matrix::t(super$GE),
    #                               genes.use = super$genes.use,
    #                               supercell_size = super$supercell_size,
    #                               k = 20)
    # D <- dist(super.PCA$x)
    # super.clusters <- supercell_cluster(D = D, k = 5, supercell_size = super$supercell_size) 
    # super$clustering <- super.clusters$clustering
    # 
    # map.cluster.to.cell.line <- supercell_assign(supercell_membership = super$clustering, clusters  = super$cell_line)
    # super$clustering_reordered <- map.cluster.to.cell.line[super$clustering]
    
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

superCells_DE_by_cluster <- function(data,  # gene expression matrix counts
                                     gamma  # graining level of super cells
)
{
    supercells_clusters <- list()
    for (cluster in names(table(Idents(data)))){
        ge_matrix <- GetAssayData(data)[, Idents(data) == cluster]
        
        super <- SCimplify(ge_matrix,
                           k.knn = 5,
                           gamma = gamma,
                           n.var.genes = 1000,
                           directed = FALSE
        )
        
        super$cell_line <- supercell_assign(clusters = rep(cluster, ncol(ge_matrix)),
                                            supercell_membership = super$membership,
                                            method = "jaccard")
        
        super$GE <- supercell_GE(ge_matrix, super$membership)
        supercells_clusters[[cluster]] <- super
    }
    
    # Merging of all supercells together
    supercells <- supercells_clusters[[1]]
    for(cluster in names(supercells_clusters)[2:length(names(supercells_clusters))]){
        supercells$GE <- cbind(supercells$GE, supercells_clusters[[cluster]]$GE)
        supercells$membership <- c(supercells$membership, supercells_clusters[[cluster]]$membership + length(supercells$membership))
        
        supercells$cell_line <- c(supercells$cell_line, supercells_clusters[[cluster]]$cell_line)
        supercells$supercell_size <- c(supercells$supercell_size, supercells_clusters[[cluster]]$supercell_size)
    }
    
    
    # --------------------------------------------------------------------------
    
    # DE
    clusters <- unique(supercells$cell_line)
    
    markers_super_per_cluster <- supercell_FindMarkers(ge = supercells$GE,
                                                          supercell_size = supercells$supercell_size,
                                                          clusters = supercells$cell_line,
                                                          ident.1 = clusters[grep('^treat', clusters)],
                                                          ident.2 = clusters[grep('^ctrl', clusters)],
                                                          logfc.threshold = 0,
                                                          only.pos = F,
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
    

    total <- total_super_markers <- total_super_markers %>%
        arrange(adj.p.value, 1 / (abs(logFC) + 1))
    
    return(total)
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