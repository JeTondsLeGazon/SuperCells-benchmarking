# 3 February 2022

# Contains all differential expression computation functions for pseudo-bulk

computeDE_pseudobulk <- function(data, test, save.path){
    message('Computing Pseudo-bulk DE genes')
    pseudo_markers <- compute_DE_bulk(data)
    saveRDS(pseudo_markers, file.path(results_folder, "pseudoMarkers.rds"))
    message('Done computing Pseudo-bulk DE genes')
    saveMarkers(markers = super_markers, 
                algo = test,
                split.by = NULL,
                base.path = save.path,
                kind = 'pseudo')
}


message('Computing Pseudo-bulk DE genes manually')
pseudo_markers_manual <- find_markers(pseudobulk_norm, stat.test) %>%
    arrange(adj.p.value, 1 / (abs(logFC) + 1)) %>%
    mutate(gene = row.names(.)) %>%
    subset(logFC > 0)
saveRDS(pseudo_markers_manual, file.path(results_folder, "pseudoMarkersManual.rds"))
message('Done computing Pseudo-bulk DE genes manually')