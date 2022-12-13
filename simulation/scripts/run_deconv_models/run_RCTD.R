#setwd("/houston_20t//alexw/ST/iSort/simulation_finer_resolution/")
library(spacexr)
library(Matrix)
library(testit)
#library(devtools)
# devtools::install_github("dmcable/RCTD", build_vignettes = TRUE)
create_RCTD_input <- function(mixtures, sig_matrix, meta){
  #' Read and return mixture object and cell type signature object from file names for RCTD.
  #' 
  #' @param mixtures data frame of the mixtures, row names are genes, columns: sample IDs.
  #' @param sig_fn data frame of signature matrix, row names are genes, columns: cell types (can be overlapped).
  #' @param meta data frame of the meta data of the mixtures, first column ID, columns denotes coordinates.
  #' @return RCTD input object containing the puck and reference.
  
  # Make sure all spots in the mixtures are in the meta data
  assert("Not all spot IDs are in the meta index", {
    mix_spots = colnames(mixtures)
    meta_spots = row.names(meta)
    all(mix_spots %in% meta_spots)})
  
  # Reorder the meta data
  meta <- meta[colnames(mixtures),]
  colnames(meta) <- c("x", "y")
  
  # Read in the cell type signatures
  overlap_genes <- intersect(row.names(mixtures), row.names(sig_matrix))
  sig_matrix <- sig_matrix[overlap_genes, ]
  
  # Create reference 
  clusters <- colnames(sig_matrix)
  colnames(sig_matrix) <- paste("cell", 1:dim(sig_matrix)[2], sep="")
  cell_types <- clusters
  names(cell_types) <- colnames(sig_matrix)
  cell_types <- as.factor(cell_types)
  sig.nUMIs <- apply(sig_matrix, 2, sum)
  names(sig.nUMIs) <- colnames(sig_matrix)
  reference <- Reference(sig_matrix, cell_types, sig.nUMIs)
  
  mix.nUMIs <- colSums(mixtures)
  puck <- SpatialRNA(meta, mixtures, mix.nUMIs)
  barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 
  puck <- restrict_puck(puck, barcodes)
  
  return(list(puck=puck, reference=reference))
}
run_RCTD <- function(puck, reference, max_cores = 10, CELL_MIN_INSTANCE = 5){
  myRCTD <- create.RCTD(puck, reference, max_cores = max_cores, CELL_MIN_INSTANCE = 5)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  results <- myRCTD@results
  rctd.props = as.data.frame(as.matrix(results$weights))
  return(rctd.props)
}
createRef <- function(count_fn, meta_fn, col='bio_celltype') {
  refA <-read.csv(count_fn, row.names=1)
  refA <- t(refA)
  metaA <- read.csv(meta_fn, row.names=1)
  ctA <- factor(metaA[, col]); names(ctA) <- row.names(metaA)
  ctA <- as.factor(ctA)
  nUMIA <- colSums(refA)
  all(names(ctA) == names(refA))
  reference <- Reference(refA, ctA, nUMIA)
  return(reference)
}
make_mixture_RCTD <- function(mixture_count_fn){
  mixture = read.csv(mixture_count_fn,row.names=1)
  xs <- c()
  ys <- c()
  inds <- row.names(mixture)
  for (i in inds) {
    x <- as.integer(strsplit(i, "x")[[1]][1])
    y <- as.integer(strsplit(i, "x")[[1]][2])
    xs <- c(xs, x)
    ys <- c(ys, y)
  }
  
  mixture_meta <- data.frame(x=xs, y=ys, row.names=inds)
  mixture = t(mixture)
  mixture.nUMI = colSums(mixture)
  mixture.spatial <- SpatialRNA(mixture_meta, mixture, mixture.nUMI)
  mixture.barcodes <- colnames(mixture.spatial@counts)
  mixture.spatial <- restrict_puck(mixture.spatial, mixture.barcodes)
  return(mixture.spatial)
}

# run_RCTD_project <- function(projDir, outDir, col){
#   # Create folder to save results
#   # if (col == 'bio_celltype'){
#   #   outDir <- file.path(projDir, 'RCTD_minor_results')
#   # } else {
#   #   outDir <- file.path(projDir, 'RCTD_major_results')
#   # }
#   
#   dir.create(outDir, showWarnings = FALSE)
#   
#   curwd <- getwd()
#   
#   mixture_count_fn <- file.path(projDir, 'simulated_mixture_raw_counts.csv')
#   
#   # refM_count_fn <- file.path(projDir, 'MIST_reference_count.csv')
#   # refM_meta_fn <- file.path(projDir, 'MIST_reference_meta.csv')
#   # 
#   # refA_count_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_count_internal.csv')
#   # refA_meta_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_meta_internal.csv')
#   # 
#   # refB_count_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_count_external.csv')
#   # refB_meta_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_meta_external.csv')
#   
#   # refI_count_fn <- file.path(curwd, 'single_cell_ref/immune_only_ref_count_external.csv')
#   # refI_meta_fn <- file.path(curwd, 'single_cell_ref/immune_only_ref_meta_external.csv')
#   # 
#   # refH_count_fn <- file.path(curwd, 'single_cell_ref/MIST_hybrid_immuneB_cpm.csv')
#   # refH_meta_fn <- file.path(curwd, 'single_cell_ref/MIST_hybrid_immuneB_meta.csv')
#   
#   # refI_count_fn <- file.path(curwd, 'single_cell_ref/immune_corrected_external_ref_count.csv')
#   # refI_meta_fn <- file.path(curwd, 'single_cell_ref/immune_corrected_external_ref_meta.csv')
#   refH_count_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_count_corrected_external.csv')
#   refH_meta_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_meta_corrected_external.csv')
#   ### END HERE ###
#   
#   
#   ## Read and process mixtures
#   mixture = read.csv(mixture_count_fn,row.names=1)
#   xs <- c()
#   ys <- c()
#   inds <- row.names(mixture)
#   for (i in inds) {
#       x <- as.integer(strsplit(i, "x")[[1]][1])
#       y <- as.integer(strsplit(i, "x")[[1]][2])
#       xs <- c(xs, x)
#       ys <- c(ys, y)
#   }
#   
#   mixture_meta <- data.frame(x=xs, y=ys, row.names=inds)
#   mixture = t(mixture)
#   mixture.nUMI = colSums(mixture)
#   mixture.spatial <- SpatialRNA(mixture_meta, mixture, mixture.nUMI)
#   mixture.barcodes <- colnames(mixture.spatial@counts)
#   mixture.spatial <- restrict_puck(mixture.spatial, mixture.barcodes)
#   
#   # reference.A <- createRef(refA_count_fn, refA_meta_fn, col)
#   # st = Sys.time()
#   # props.A <- run_RCTD(mixture.spatial,reference.A, max_cores = 20, CELL_MIN_INSTANCE = 5)
#   # et = Sys.time()
#   # print(et-st)
#   # write.csv(props.A, file.path(outDir, 'estimated_proportions_ref_internal.csv'))
#   # print('Internal reference RCTD done.')
#   # 
#   # reference.B <- createRef(refB_count_fn, refB_meta_fn, col)
#   # st = Sys.time()
#   # props.B <- run_RCTD(mixture.spatial,reference.B, max_cores = 20, CELL_MIN_INSTANCE = 5)
#   # et = Sys.time()
#   # print(et-st)
#   # write.csv(props.B, file.path(outDir, 'estimated_proportions_ref_external.csv'))
#   # print('External reference RCTD done.')
#   
#   ### Use MIST as reference
#   # if (col == 'bio_celltype_major'){
#   #   reference.MIST <- createRef(refM_count_fn, refM_meta_fn, col)
#   #   st = Sys.time()
#   #   props.MIST <- run_RCTD(mixture.spatial,reference.MIST, max_cores = 20, CELL_MIN_INSTANCE = 5)
#   #   et = Sys.time()
#   #   print(et-st)
#   #   write.csv(props.MIST, file.path(outDir, 'estimated_proportions_ref_MIST.csv'))
#   #   print('MIST reference RCTD done.')
#   # }
#   
#   reference.H <- createRef(refH_count_fn, refH_meta_fn, 'bio_celltype')
#   st = Sys.time()
#   props.H <- run_RCTD(mixture.spatial,reference.H, max_cores = 20, CELL_MIN_INSTANCE = 5)
#   et = Sys.time()
#   print(et-st)
#   write.csv(props.H, file.path(outDir, "estimated_proportions_external_corrected.csv"))
#   print('Immune only reference RCTD done.')
#   
#   # reference.H <- createRef('/houston_20t/alexw/ST/iSort/simulation_finer_resolution/scripts/integration/hybrid_reference_with_external_corrected_all_cell_types.csv',
#   #                          '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/scripts/integration/hybrid_reference_with_external_corrected_all_cell_types_meta.csv', 
#   #                          "bio_celltype")
#   # st = Sys.time()
#   # props.H <- run_RCTD(mixture.spatial,reference.H, max_cores = 20, CELL_MIN_INSTANCE = 5)
#   # et = Sys.time()
#   # print(et-st)
#   # write.csv(props.H, '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/scripts/integration/hybrid_reference_all_celltypes_results.csv')
#   #write.csv(props.H, file.path(projDir, 'RCTD_estimated_proportions_external_mist_hybrid_external_immune.csv'))
#   #print('Hybrid reference RCTD done.')
#   
# }
# ### CHANGE FOR DIFFERENT SCINERIO ####
# projDir <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/data/cancer_inf_0.1'
# outDir <- file.path(projDir, 'RCTD_minor_results/')
# 
# ref_count_fn <- file.path(workDir, 'single_cell_ref/single_cell_ref_count_corrected_augmented_external.csv')
# ref_meta_fn <- file.path(workDir, 'single_cell_ref/single_cell_ref_meta_corrected_augmented_external.csv')
# mixture_count_fn <- file.path(projDir, 'simulated_mixture_raw_counts.csv')
# mixture <- make_mixture_RCTD(mixture_count_fn)
# reference <- createRef(ref_count_fn, ref_meta_fn)
# mixture.props <- run_RCTD(mixture, reference)
# write.csv(mixture.props, file.path(projDir, 'RCTD_minor_results/estimated_proportions_external_corrected_augmented.csv'))
# #system.time(run_RCTD_project(projDir2, 'bio_celltype'))
# #system.time(run_RCTD_project(projDir, outDir, 'bio_celltype'))
# mixture_count_fn
