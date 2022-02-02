#9th November 2021

# Wrapper for the Supercells functionnality to perform QC and DE analysis
library(SuperCell)

# Uses new cell.annotation to perform DE for each cluster according to treatment
# vs control and merge all DEs together for a single gamma
superCells_DE <- function(data,  # gene expression matrix counts
                          gamma,  # graining level of super cells
                          weighted,  # whether to perform weighted t-tests or not
                          test.use,
                          by.group = F,
                          split.by = 'sample',
                          bm = T
)
{
    super <- superCellWrapper(data = data, 
                              gamma = gamma, 
                              split.by = split.by,
                              arithmetic = T,
                              SC.type = 'Exact',
                              bm = bm)
    
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
                            test.use = test.use)
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
                     test.use = test.use)
    }
    return(arrangeDE(DEs))
}


# Computes DE genes for superCells at different graining levels
superCells_DEs <- function(data,  # normalized logcounts seurat object
                           gammas,  # list of graning levels to use
                           knn,  # number of nearest neighbors to use
                           weighted = F,
                           test.use = 't', 
                           split.by = 'sample',
                           bm = T) 
{
    super_DEs <- list()
    for(gam in gammas){
            super_res <- superCells_DE(data,
                                       gamma = gam,
                                       weighted = weighted,
                                       test.use = test.use, 
                                       split.by = split.by,
                                       bm = bm)
        super_DEs[[as.character(gam)]] <- super_res
    }
    return(super_DEs)
}


# Create super cells using arithmetic or geometric average for
# the gene expression matrix
createSuperCells <- function(data, 
                             gamma, 
                             arithmetic = TRUE, 
                             split.by = 'sample'){
    
    split.condition <- data[[split.by]][[1]]
    super <-  SCimplify(GetAssayData(data),
                        cell.annotation = split.condition,
                        k.knn = 5,
                        gamma = gamma,
                        n.var.genes = 1000,
                        directed = FALSE
    )
    
    super$cell_line <- supercell_assign(clusters = data$label,
                                        supercell_membership = super$membership,
                                        method = "jaccard")
    
    # arithmetic average
    if(arithmetic){
        super$GE <- log(supercell_GE(expm1(GetAssayData(data)), super$membership) + 1)
    }else{  # geometric average
        super$GE <- supercell_GE(GetAssayData(data), super$membership)
    }
    return(super)
}


# Create supercells following benchmarking functions from Mariia
# Uses package supercellsBM
createSuperCellsBM <- function(data, 
                               gamma, 
                               arithmetic = TRUE, 
                               split.by = 'sample',
                               SC.type = 'Exact'){
    SC.list <- compute_supercells(
        sc = data,
        ToComputeSC = T,
        data.folder = 'data/',
        filename = paste0('superCells', gamma),
        gamma.seq = c(gamma),
        n.var.genes = 1000,
        k.knn = 5,
        n.pc = 10,
        approx.N = 1000,
        fast.pca = TRUE,
        genes.use = NULL, 
        genes.exclude = NULL,
        seed.seq = c(0),
        split.by = split.by
    )
    super <- SC.list[[SC.type]][[as.character(gamma)]][[1]]
    super$cell_line <- supercell_assign(clusters = data$label,
                                        supercell_membership = super$membership,
                                        method = "jaccard")
    if(arithmetic){
        super$GE <- log(supercell_GE(expm1(GetAssayData(data)), super$membership) + 1)
    }else{  # geometric average
        super$GE <- supercell_GE(GetAssayData(data), super$membership)
    }
    return(super)
}


# Wrapper function for createSuperCells and createSuperCellsBM
superCellWrapper <- function(data, 
                             gamma, 
                             arithmetic = TRUE, 
                             split.by = 'sample',
                             SC.type = 'Exact',
                             bm = T,
                             norm = T){
    res <- load_superCell(gamma = gamma,
                          arithmetic = arithmetic,
                          split.by = split.by,
                          SC.type = SC.type,
                          bm = bm,
                          norm = norm)
    
    if(!is.null(res)){
    #    return(res)
    }
    if(bm){  # benchmarking from Mariia
        super <- createSuperCellsBM(data = data,
                                    gamma = gamma,
                                    arithmetic = arithmetic,
                                    split.by = split.by,
                                    SC.type = SC.type)
    }else{
        super <- createSuperCells(data = data,
                                  gamma = gamma,
                                  arithmetic = arithmetic,
                                  split.by = split.by)
    }
    save_superCell(super = super,
                   gamma = gamma,
                   arithmetic = arithmetic,
                   split.by = split.by,
                   SC.type = SC.type,
                   bm = bm,
                   norm = norm)
    return(super)
}


# Save supercell object to corresponding filename
save_superCell <- function(super,
                           gamma,
                           arithmetic,
                           split.by,
                           SC.type,
                           bm, norm){
    filename <- createFilename(gamma = gamma,
                               arithmetic = arithmetic,
                               split.by = split.by,
                               SC.type = SC.type,
                               bm = bm,
                               norm = norm)
    path <- file.path('data', 'SC')
    if(!dir.exists(path)){
        dir.create(path, recursive = T)
    }
    saveRDS(super, file.path(path, filename))
}


# Check if corresponding superCell already exist in data/SC/ and load it
load_superCell <- function(gamma,
                           arithmetic,
                           split.by,
                           SC.type,
                           bm,
                           norm){
   filename <- createFilename(gamma = gamma,
                              arithmetic = arithmetic,
                              split.by = split.by,
                              SC.type = SC.type,
                              bm = bm,
                              norm = norm)
    path <- file.path('data', 'SC')
    if(!dir.exists(path)){
        dir.create(path, recursive = T)
    }
    if(file.exists(file.path(path, filename))){
        return(readRDS(file.path(path, filename)))
    }else{
        return(NULL)
    }
}


# Create corresponding filename for different kinds of supercells
createFilename <- function(gamma,
                           arithmetic,
                           split.by,
                           SC.type,
                           bm,
                           norm){
    gamma.ch <- as.character(gamma)
    a_g <- ifelse(arithmetic, 'a', 'g')
    bm_std <- ifelse(bm, 'bm', 'std')
    norm_raw <- ifelse(norm, 'norm', 'raw')
    SC.type <- tolower(SC.type)
    filename <- paste(gamma.ch, a_g, bm_std, SC.type, split.by, norm_raw, sep = '_')
    filename <- paste0(filename, '.rds')
    return(filename)
}