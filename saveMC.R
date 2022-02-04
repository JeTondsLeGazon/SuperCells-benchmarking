library(metacell)
library(config)

find_common_genes <- function(files, mc_folder, samples, mc.type){
  genes.per.gamma <- list()
  for(file in files){
    data <- readRDS(file.path(mc_folder, file))[[mc.type]]
    gammas <- as.numeric(names(data))
    if(length(genes.per.gamma) == 0){
      genes.per.gamma <- sapply(seq_along(gammas), function(x) list())
    }
    for(i in seq_along(gammas)){
        ge <- data[[i]][[1]]$mc_info$mc@mc_fp

      
      #find intersection of all genes
      if(length(genes.per.gamma[[i]]) == 0){
        genes.per.gamma[[i]] <- rownames(ge)
      }else{
        genes.per.gamma[[i]] <- intersect(genes.per.gamma[[i]],
                                          rownames(ge))
      }
    }
  }
  return(genes.per.gamma)
}


create_metacell <- function(genes, samples, files, mc_folder, mc.type){
  MC <- sapply(seq_len(20), function(x) list())
  MC.mcClass <- sapply(seq_len(20), function(x) list())
  
  # Rearrange data
  for(samp in samples){
    id <- grep(samp, files)
    data <- readRDS(file.path(mc_folder, files[id]))[[mc.type]]
    gammas <- as.numeric(names(data))
    for(i in seq_along(gammas)){
      keep.genes <- genes[[i]]
      ge <- data[[i]][[1]]$mc_info$mc@mc_fp[keep.genes, ]
      e_gc <- data[[i]][[1]]$mc_info$mc@e_gc[keep.genes, ]
      if(length(MC[[i]]) ==  0){
        MC[[i]]$ge <- ge
        MC[[i]]$sample <- rep(samp, ncol(ge))
        MC[[i]]$size <- table(data[[i]][[1]]$mc_info$mc@mc)
        MC[[i]]$membership <- data[[i]][[1]]$mc_info$mc@mc
        MC[[i]]$e_gc <- e_gc
        MC.mcClass[[i]] <- data[[i]][[1]]$mc_info$mc
      }else{
        MC[[i]]$ge <- cbind(MC[[i]]$ge,
                                    ge)
        MC[[i]]$sample <- c(MC[[i]]$sample,
                                    rep(samp, ncol(ge)))
        MC[[i]]$size <- c(MC[[i]]$size, 
                                  table(data[[i]][[1]]$mc_info$mc@mc))
        names(MC[[i]]$size) <- seq_along(MC[[i]]$size)
        MC[[i]]$membership <- c(MC[[i]]$membership,
                                        data[[i]][[1]]$mc_info$mc@mc + max(MC[[i]]$membership))
        MC[[i]]$e_gc <- cbind(MC[[i]]$e_gc,
                              e_gc)
        MC.mcClass[[i]]@mc <- c(MC.mcClass[[i]]@mc,
                                        
                                        data[[i]][[1]]$mc_info$mc@mc + max(MC.mcClass[[i]]@mc))
      }
    }
  }
  MC <- MC[sapply(MC, function(x) length(x) > 0)]
  
  # Rename with gammas
  for(i in seq_along(MC)){
    ge <- MC[[i]]$ge
    new.colnames <- paste(MC[[i]]$sample, seq_len(ncol(ge)), sep = '_')
    colnames(MC[[i]]$ge) <- new.colnames
    MC[[i]]$gamma <- round(N.sc / ncol(ge))
    
    MC.mcClass[[i]] <- mc_compute_fp(mc = MC.mcClass[[i]],
                                             us = data_raw@assays$RNA@counts[,names(MC.mcClass[[i]]@mc)])
    
    MC[[i]]$mc_fp_updated <- MC.mcClass[[i]] 
  }
  
  names(MC) <- sapply(MC, function(x) x$gamma)
  return(MC)
}

# ------------------------------------------------------------------------------
# ----------------------     MAIN    -------------------------------------------
# ------------------------------------------------------------------------------
mc.type <- 'Metacell_default'
config <- config::get(file = 'configs/hagai_mouse_lps_config.yml')
sample_vs_condition <- paste('mc', config$splitBy, sep = '_')
mc_folder <- file.path('data', config$intermediaryDataFile, sample_vs_condition)
files <- list.files(mc_folder)
files.used <- files[grep('MC', files)]
samples_condition <- sapply(files.used, function(x) strsplit(x, '_')[[1]][1])
#sc_data <- readRDS('data/hagai_mouse_lps_data/singleCellFiltered.rds')
N.sc <- 16882  # TODO: change this

# Default
genes <- find_common_genes(files = files.used,
                  mc_folder = mc_folder,
                  samples = samples_condition,
                  mc.type = mc.type)

MC <- create_metacell(genes = genes,
                              samples = samples_condition,
                              files = files.used,
                              mc_folder = mc_folder,
                              mc.type = mc.type)
if(mc.type == 'Metacell_default'){
  filename <- 'mc_default.rds'
}else{
  filename <- 'mc_SC_like.rds'
}
saveRDS(MC, file.path(mc_folder, filename))
