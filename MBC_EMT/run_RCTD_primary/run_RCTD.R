setwd("/home/humble_local_25t/alexw/projects/ShawnLab/EMS_deconv_0307/scripts/run_models/")

#library(devtools)
# devtools::install_github("dmcable/RCTD", build_vignettes = TRUE)
.libPaths(c("/home/humble_local_25t/alexw/R/x86_64-pc-linux-gnu-library/4.0",
            "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/home/humble_local2_25t/shared_R_libraries"))
library(spacexr)
library(Matrix)
library(testit)

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


fn <- '/home/humble_local_25t/alexw/projects/ShawnLab/EMS_deconv_0307/output_files/24679/ReSort_mixture_raw.csv'
mixture <- make_mixture_RCTD(fn)
ref_count_fn <-  '/home/humble_local_25t/alexw/projects/ShawnLab/EMS_deconv_0307/output_files/24679/ReSort_reference_raw_count.csv'
ref_meta_fn <-  '/home/humble_local_25t/alexw/projects/ShawnLab/EMS_deconv_0307/output_files/24679/ReSort_reference_meta.csv'
reference <- createRef(ref_count_fn, ref_meta_fn)
ems.props <- run_RCTD(mixture, reference)
#write.csv(ems.props, file.path(projDir, 'RCTD_props.csv'))
write.csv(ems.props, file.path('../data/ReSort_primary/primary_props.csv'))