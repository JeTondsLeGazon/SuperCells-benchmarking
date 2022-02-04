#9th November 2021

# Wrapper for the Supercells functionnality to perform QC and DE analysis
library(SuperCell)


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
    
    # Arithmetic average based on normalized (non-log) counts, then take logcount
    super$GE <- log1p(supercell_GE(expm1(data@assays$RNA@data), super$membership))
    
    # Same as above but non-log to input into DESeq2 and EdgeR after multiplication
    # by supercell_size
    super$counts <- supercell_GE(expm1(data@assays$RNA@data), super$membership)
    return(super)
}


# Wrapper function for createSuperCells and createSuperCellsBM
superCellWrapper <- function(data, 
                             gamma, 
                             arithmetic = TRUE, 
                             split.by = 'sample',
                             SC.type = 'Exact',
                             norm = T){
    res <- load_superCell(gamma = gamma,
                          arithmetic = arithmetic,
                          split.by = split.by,
                          SC.type = SC.type,
                          norm = norm)
    
    if(!is.null(res)){
    #    return(res)
    }
    super <- createSuperCellsBM(data = data,
                                gamma = gamma,
                                arithmetic = arithmetic,
                                split.by = split.by,
                                SC.type = SC.type)
   
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
                           norm){
   filename <- createFilename(gamma = gamma,
                              arithmetic = arithmetic,
                              split.by = split.by,
                              SC.type = SC.type,
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
                           norm){
    gamma.ch <- as.character(gamma)
    a_g <- ifelse(arithmetic, 'a', 'g')
    norm_raw <- ifelse(norm, 'norm', 'raw')
    SC.type <- tolower(SC.type)
    filename <- paste(gamma.ch, a_g, SC.type, split.by, norm_raw, sep = '_')
    filename <- paste0(filename, '.rds')
    return(filename)
}