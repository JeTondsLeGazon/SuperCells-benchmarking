#9th November 2021

# Wrapper for the Supercells functionnality to perform QC and DE analysis
library(SuperCell)

# Uses new cell.annotation to perform DE for each cluster according to treatment
# vs control and merge all DEs together for a single gamma
superCells_DE <- function(data,  # gene expression matrix counts
                          gamma,  # graining level of super cells
                          weighted,  # whether to perform weighted t-tests or not
                          test.use,
                          by.group = F
)
{
    
    super <- superCells_GE(data, gamma)
    
    # DE per cluster group
    clusters <- unique(super$cell_line)
    nb_groups <- sapply(clusters, function(x) as.numeric(str_sub(x, -1, -1)))
    
    if (weighted){
        myFunc <- supercell_FindMarkers_weighted
    }else{
        myFunc <- supercell_FindMarkers
    }
    
    if(by.group){
        
        DEs <- c()
        for(i in seq_along(max(nb_groups)))
        {
            DE <- myFunc(ge = super$GE,
                            supercell_size = super$supercell_size,
                            clusters = super$cell_line,
                            ident.1 = clusters[grep(paste0('^treat.+', i, '$'), clusters)],
                            ident.2 = clusters[grep(paste0('^ctrl.+', i, '$'), clusters)],
                            logfc.threshold = 0,
                            only.pos = F,
                            do.bootstrapping = F,
                            test.use = test.use) %>%
                    mutate(gene = rownames(.))
            
            DEs <- rbind(DEs, DE)
        }
    }else{
        DEs <- myFunc(ge = super$GE,
                     supercell_size = super$supercell_size,
                     clusters = super$cell_line,
                     ident.1 = clusters[grep('^treat', clusters)],
                     ident.2 = clusters[grep('^ctrl', clusters)],
                     logfc.threshold = 0,
                     only.pos = F,
                     do.bootstrapping = F,
                     test.use = test.use) %>%
            mutate(gene = rownames(.))
    }
    DEs <- DEs %>%
        arrange(adj.p.value, 1 / (abs(logFC) + 1), T) %>%
        subset(!duplicated(.))
    rownames(DEs) <- DEs$gene
    return(DEs)
}


# Computes DE genes for superCells at different graining levels
superCells_DEs <- function(data,  # normalized logcounts seurat object
                           gammas,  # list of graning levels to use
                           knn,  # number of nearest neighbors to use
                           weighted = F,
                           test.use = 't') 
{
    super_DEs <- list()
    for(gam in gammas){
            super_res <- superCells_DE(data,
                                       gamma = gam,
                                       weighted = weighted,
                                       test.use = test.use)
        super_DEs[[as.character(gam)]] <- super_res
    }
    return(super_DEs)
}


# Computes the gene expression (arithmetic average) matrix for superCell at a 
# specific gamma
superCells_GE <- function(data, gamma){
    super <-  SCimplify(GetAssayData(data),
                        cell.annotation = Idents(data),
                        k.knn = 5,
                        gamma = gamma,
                        n.var.genes = 1000,
                        directed = FALSE
    )
    
    super$cell_line <- supercell_assign(clusters = Idents(data),
                                        supercell_membership = super$membership,
                                        method = "jaccard")
    
    super$GE <- log(supercell_GE(expm1(GetAssayData(data)), super$membership) + 1)
    return(super)
}
