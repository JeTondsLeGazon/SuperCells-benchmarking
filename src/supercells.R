#9th November 2021

# Wrapper for the Supercells functionnality to perform QC and DE analysis


# Create supercells following benchmarking functions from Mariia
# Uses package supercellsBM
createSuperCellsBM <- function(data, 
                               gamma,
                               data_folder,
                               arithmetic = TRUE, 
                               split.by = 'sample',
                               SC.type = 'Exact',
                               force_compute = TRUE){
    
    filename_no_extension <- paste('superCells', gamma, split.by, sep = '_')
    filename <- paste0(filename_no_extension, '.Rds')
    ToComputeSC <- force_compute | !file.exists(file.path(data_folder, 'SC', filename))
    SC.list <- compute_supercells(
        sc = data,
        ToComputeSC = ToComputeSC,
        data.folder = data_folder,
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
    super$label <- supercell_assign(clusters = data$label,
                                    supercell_membership = super$membership,
                                    method = "jaccard")
    
    # Arithmetic average based on normalized (non-log) counts, then take logcount
    super$GE <- log1p(supercell_GE(expm1(data@assays$RNA@data), super$membership))
    
    # Counts for DESeq2 and EdgeR (raw counts as normalization will happen in
    # both algorithms later)
    super$counts <- supercell_GE(data@assays$RNA@counts, super$membership)
    return(super)
}
