library(metacell)
library(Seurat)
mc_folder <- 'data/mc/'
files <- list.files(mc_folder)
samples <- sapply(list.files(mc_folder), function(x) strsplit(x, '_')[[1]][1])
sc_data <- readRDS('data/hagai_mouse_lps_data/singleCellFiltered.rds')
N.sc <- ncol(sc_data)

# Default
MC.default <- sapply(seq_len(20), function(x) list())
genes.per.gamma <- sapply(seq_along(files), function(x) list())
for(samp in samples){
  id <- grep(samp, files)
  data <- readRDS(file.path(mc_folder, files[id]))
  gammas <- as.numeric(names(data$Metacell_default))
  for(i in seq_along(gammas)){
    ge <- data$Metacell_default[[i]][[1]]$mc_info$mc@mc_fp
    
    #find intersection of all genes
    if(length(genes.per.gamma[[i]]) == 0){
      genes.per.gamma[[i]] <- rownames(ge)
    }else{
      genes.per.gamma[[i]] <- intersect(genes.per.gamma[[i]],
                                                          rownames(ge))
    }
  }
}
for(samp in samples){
  id <- grep(samp, files)
  data <- readRDS(file.path(mc_folder, files[id]))
  gammas <- as.numeric(names(data$Metacell_default))
  for(i in seq_along(gammas)){
    keep.genes <- genes.per.gamma[[i]]
    ge <- data$Metacell_default[[i]][[1]]$mc_info$mc@mc_fp[keep.genes, ]
    if(length(MC.default[[i]]) ==  0){
      MC.default[[i]]$ge <- ge
      MC.default[[i]]$sample <- rep(samp, ncol(ge))
    }else{
      MC.default[[i]]$ge <- cbind(MC.default[[i]]$ge,
                                        ge)
      MC.default[[i]]$sample <- c(MC.default[[i]]$sample,
                                        rep(samp, ncol(ge)))
    }
  }
}
MC.default <- MC.default[sapply(MC.default, function(x) length(x) > 0)]
for(i in seq_along(MC.default)){
  ge <- MC.default[[i]]$ge
  new.colnames <- paste(MC.default[[i]]$sample, seq_len(ncol(ge)), sep = '_')
  colnames(MC.default[[i]]$ge) <- new.colnames
  MC.default[[i]]$gamma <- round(N.sc / ncol(ge))
}

names(MC.default) <- sapply(MC.default, function(x) x$gamma)
saveRDS(MC.default, 'data/mc/mc_default.rds')
