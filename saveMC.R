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
  MC.default <- sapply(seq_len(20), function(x) list())
  
  # Rearrange data
  for(samp in samples){
    id <- grep(samp, files)
    data <- readRDS(file.path(mc_folder, files[id]))[[mc.type]]
    gammas <- as.numeric(names(data))
    for(i in seq_along(gammas)){
      keep.genes <- genes[[i]]
      ge <- data[[i]][[1]]$mc_info$mc@mc_fp[keep.genes, ]
      if(length(MC.default[[i]]) ==  0){
        MC.default[[i]]$ge <- ge
        MC.default[[i]]$sample <- rep(samp, ncol(ge))
        MC.default[[i]]$size <- table(data[[i]][[1]]$mc_info$mc@mc)
        MC.default[[i]]$membership <- data[[i]][[1]]$mc_info$mc@mc
      }else{
        MC.default[[i]]$ge <- cbind(MC.default[[i]]$ge,
                                    ge)
        MC.default[[i]]$sample <- c(MC.default[[i]]$sample,
                                    rep(samp, ncol(ge)))
        MC.default[[i]]$size <- c(MC.default[[i]]$size, 
                                  table(data[[i]][[1]]$mc_info$mc@mc))
        names(MC.default[[i]]$size) <- seq_along(MC.default[[i]]$size)
        MC.default[[i]]$membership <- c(MC.default[[i]]$membership,
                                        data[[i]][[1]]$mc_info$mc@mc + max(MC.default[[i]]$membership))
      }
    }
  }
  MC.default <- MC.default[sapply(MC.default, function(x) length(x) > 0)]
  
  # Rename with gammas
  for(i in seq_along(MC.default)){
    ge <- MC.default[[i]]$ge
    new.colnames <- paste(MC.default[[i]]$sample, seq_len(ncol(ge)), sep = '_')
    colnames(MC.default[[i]]$ge) <- new.colnames
    MC.default[[i]]$gamma <- round(N.sc / ncol(ge))
    print(MC.default[[i]]$gamma)
    print(round(N.sc / max(MC.default[[i]]$membership)))
  }
  
  names(MC.default) <- sapply(MC.default, function(x) x$gamma)
  return(MC.default)
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

MC.default <- create_metacell(genes = genes,
                              samples = samples_condition,
                              files = files.used,
                              mc_folder = mc_folder,
                              mc.type = mc.type)
saveRDS(MC.default, file.path(mc_folder, mc.type))
