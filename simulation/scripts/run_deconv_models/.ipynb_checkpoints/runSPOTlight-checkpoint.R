#devtools::install_github("https://github.com/MarcElosua/SPOTlight")
#install.packages("dplyr")
#install.packages('Seurat')
workDir <- "/houston_20t/alexw/ST/iSort/simulation_finer_resolution/"
setwd(workDir)

library(Seurat)
library(SPOTlight)
library(dplyr)

run_spotlight <- function(mixtures, ref_count_fn, ref_meta_fn, col) {
  ## PROCESS REFERENCE DATA
  ref <-  t(as.matrix(read.csv(ref_count_fn, row.names=1, 
                               check.names = F, stringsAsFactors = F)))
  ref_meta <- read.csv(ref_meta_fn, row.names=1, 
                       check.names = F, stringsAsFactors = F)
  cts <- ref_meta[, col]
  colnames(ref) <- paste("spot", 1:dim(ref)[2], sep="")
  pheno_sigs <- data.frame(row.names=colnames(ref), spotID = colnames(ref), class=cts)
  sc_obj <- CreateSeuratObject(count=ref, project='sc_ref', meta.data = pheno_sigs)
  sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 1e6, margin = 1)
  Seurat::Idents(object = sc_obj) <- sc_obj@meta.data$class
  cluster_markers_all <- Seurat::FindAllMarkers(object = sc_obj, assay = "RNA",
                                                slot = "data", verbose = TRUE, only.pos = TRUE)
  cluster_markers_all <- cluster_markers_all[cluster_markers_all$avg_log2FC >= 0.59 &
                                               cluster_markers_all$p_val_adj <= 0.01, ]
  ## PROCESS SPATIAL DATA
  spatial_counts <- CreateSeuratObject(count=mixtures, project='st_mixture')
  spatial_counts <- NormalizeData(spatial_counts, normalization.method = 'LogNormalize', 
                                  scale.factor = 1e6, margin=1)
  
  ## RUN SPOTLIGHT
  spotlight_ls <- spotlight_deconvolution(
    se_sc = sc_obj,
    counts_spatial = spatial_counts@assays$RNA@counts,
    clust_vr = "class", # Variable in sc_seu containing the cell-type annotation
    cluster_markers = cluster_markers_all, # Dataframe with the marker genes
    ntop = NULL, # How many of the marker genes to use (by default all)
    transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
    method = "nsNMF", # Factorization method
    min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
  )
  
  ## Post-process results
  spotlight.props <- data.frame(spotlight_ls[[2]], row.names=colnames(mixtures))
  xs <- c()
  ys <- c()
  for (spot in row.names(spotlight.props)){
    x <- as.integer(strsplit(spot, split = "x")[[1]][1])
    y <- as.integer(strsplit(spot, split = "x")[[1]][2])
    xs <- c(xs, x)
    ys <- c(ys, y)
  }
  
  spotlight.props$coordX <- as.numeric(xs)
  spotlight.props$coordY <- as.numeric(ys)
  spotlight.props <- spotlight.props[, c(unique(cts), c('coordX', 'coordY'))]
  return(spotlight.props)
}
run_spotlight_proj <- function(projDir, col){
  ## Create folder to save results if not existed
  if (col == 'bio_celltype'){
    outDir <- file.path(projDir, 'spotlight_minor_results')
  } else {
    outDir <- file.path(projDir, 'spotlight_major_results')
  }
  dir.create(outDir, showWarnings = FALSE)
  
  mixture_count_fn <- paste(projDir, 'simulated_mixture_raw_counts.csv', sep='')
  refM_count_fn <- paste(projDir, 'MIST_reference_count.csv', sep='')
  refM_meta_fn <- paste(projDir, 'MIST_reference_meta.csv', sep='')
  
  curwd <- getwd()
  refA_count_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_count_internal.csv')
  refA_meta_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_meta_internal.csv')  
  refB_count_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_count_external.csv')
  refB_meta_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_meta_external.csv')
  
  ## Read and process mixtures
  mixture = as.matrix(read.csv(mixture_count_fn,row.names=1))
  mixture = t(mixture)
  
  # PDACA as reference
  st = Sys.time()
  props.refA <- run_spotlight(mixtures = mixture, ref_count_fn = refA_count_fn, ref_meta_fn = refA_meta_fn, col = col)
  et = Sys.time()
  print(et-st)
  write.csv(props.refA, file.path(outDir, 'estimated_proportions_ref_internal.csv'))
  print("External reference done.")
  
  # PDACB as reference
  st = Sys.time()
  props.refB <- run_spotlight(mixtures = mixture, ref_count_fn = refB_count_fn, ref_meta_fn = refB_meta_fn, col=col)
  et = Sys.time()
  print(et-st)
  write.csv(props.refB, file.path(outDir, 'estimated_proportions_ref_external.csv'))
  print("External reference done.")
  
  # MIST as reference
  
  if (col == 'bio_celltype_major'){
    st = Sys.time()
    props.refM <- run_spotlight(mixtures = mixture, ref_count_fn = refM_count_fn, ref_meta_fn = refM_meta_fn, col=col)
    et = Sys.time()
    print(et-st)
    write.csv(props.refM, file.path(outDir, 'estimated_proportions_ref_MIST.csv'))
    print("MIST reference done.")
  }

}
### INPUT FILES, TO BE CHANGED ####

projDir2 <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/data/cancer_inf_0.1/'
system.time(run_spotlight_proj(projDir2, 'bio_celltype'))
system.time(run_spotlight_proj(projDir2, 'bio_celltype_major'))
