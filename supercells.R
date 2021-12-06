#9th November 2021

# Wrapper for the Supercells functionnality to perform QC and DE analysis
library(SuperCell)


# Uses new cell.annotation to perform DE for each cluster according to treatment
# vs control and merge all DEs together for a single gamma
superCells_DE <- function(data,  # gene expression matrix counts
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
    supercell_plot(super$graph.supercells, 
                   group = super$cell_line, 
                   seed = seed, 
                   alpha = -pi/2,
                   main  = paste0("Super-cell gamma = ", gamma, " colored by cell line assignment (rotated)"))
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
            super_res <- superCells_DE(data,
                                       gamma = gam)
        super_DEs[[as.character(gam)]] <- super_res
    }
    return(super_DEs)
}