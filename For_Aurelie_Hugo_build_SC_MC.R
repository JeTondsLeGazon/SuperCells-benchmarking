#Building MC and SC for DEA
ToComputeSC <- TRUE

SC.list <- compute_supercells(
  sc.GE, # log normalized gene expression here
  ToComputeSC = ToComputeSC,
  data.folder = data.folder,
  filename = filename,
  gamma.seq = .gamma.seq,
  n.var.genes = .N.var.genes,
  k.knn = .k.knn,
  n.pc = .N.comp,
  approx.N = .approx.N,
  fast.pca = TRUE,
  genes.use = .genes.use, 
  genes.exclude = .genes.omit,
  seed.seq = .seed.seq
) 

min_mc_size_seq <- c(1, 10, 20, 30, 50)

SC.mc <- compute_supercells_metacells_with_min_mc_size(
  sc.counts = sc.counts, # single-cell count matix
  gamma.seq = .gamma.seq,
  SC.list = SC.list,
  min_mc_size_seq = min_mc_size_seq,
  proj.name = proj.name,
  ToComputeSC = ToComputeSC, 
  mc.k.knn = 100,
  T_vm_def = 0.08,
  MC.folder = "MC", 
  MC_gene_settings = c('Metacell_default', 'Metacell_SC_like') # do not change
) 

# if  you need to compute SC with additional gammas corresponding to gammas obtained with metacell

additional_gamma_seq <- get_actual_gammas_metacell(SC.mc[1])
SC.list <- compute_supercells_additional_gammas(
  SC.list,
  additional_gamma_seq = additional_gamma_seq,
  ToComputeSC = ToComputeSC,
  data.folder = data.folder,
  filename = filename,
  approx.N = .approx.N,
  fast.pca = TRUE
)

# a small trick with MC, just run it :)
SC.mc.fp <- SC.mc
names(SC.mc.fp) <- sapply(names(SC.mc), FUN = function(x){paste0(x, '_fp')})

SC.mc.av <- SC.mc
names(SC.mc.av) <- sapply(names(SC.mc), FUN = function(x){paste0(x, '_av')})

SC.mc.expanded <- c(SC.mc.fp, SC.mc.av)

names(SC.mc.expanded)
SC.mc.expanded <- SC.mc.expanded[c(1,4)]
names(SC.mc.expanded)

rm(SC.mc.fp, SC.mc.av, SC.mc)
# done with a small trick 

# merge SC and MC 
SC.list <- c(SC.list, SC.mc.expanded)
rm(SC.mc.expanded)

## Get GE for simplified data
# this function will take care of averaging gene expression for all methods, except MC_def_fp. For MC_def_fp, it will keep MC footprint computed woth metacell

SC.GE.list <- compute_supercells_GE(
  sc.GE = sc.GE, 
  SC.list = SC.list,
  ToComputeSC_GE = ToComputeSC_GE, 
  data.folder = data.folder,
  filename = filename
)