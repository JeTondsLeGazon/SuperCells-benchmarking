#9th November 2021

# Wrapper for the Supercells functionnality to perform QC and DE analysis


# Create supercells following benchmarking functions from Mariia
# Uses package supercellsBM
createSuperCellsBM <- function(data, 
                               gamma, 
                               arithmetic = TRUE, 
                               split.by = 'sample',
                               SC.type = 'Exact',
                               force_compute = TRUE){
    
    filename_no_extension <- paste('superCells', gamma, split.by, sep = '_')
    filename <- paste0(filename_no_extension, '.Rds')
    ToComputeSC <- force_compute | !file.exists(filename)
    SC.list <- compute_supercells(
        sc = data,
        ToComputeSC = ToComputeSC,
        data.folder = 'data/SC',
        filename = filename_no_extension,
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
