#workDir <- "/houston_20t/alexw/ST/iSort/simulation_finer_resolution/"
setwd(workDir)

library(Seurat)
library(Giotto)

### Helper functions
make_signature_matrix <- function(ref_count_fn, ref_meta_fn, col){
  ref <- t(as.matrix(read.csv(ref_count_fn, row.names=1)))
  ref_meta <- read.csv(ref_meta_fn, row.names=1)
  cell_types <- ref_meta[, col]
  colnames(ref) <- paste("spot", 1:dim(ref)[2], sep="")
  pheno_sigs <- data.frame(row.names=colnames(ref), cellID = colnames(ref), class=cell_types)
  sc_obj <- CreateSeuratObject(count=ref, project='sc_ref', meta.data = pheno_sigs)
  sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 1e6, margin = 1)
  Seurat::Idents(object = sc_obj) <- sc_obj@meta.data$class
  cluster_markers_all <- Seurat::FindAllMarkers(object = sc_obj, assay = "RNA",
                                                slot = "data", verbose = TRUE, only.pos = TRUE)
  
  cluster_markers_all <- cluster_markers_all[(cluster_markers_all$avg_log2FC >= 0.26) &
                                               (cluster_markers_all$p_val_adj <= 0.05), ]
  sign_genes <- c()
  for (ct in unique(cluster_markers_all$cluster)) {
    ct_marker_df <- cluster_markers_all[cluster_markers_all$cluster==ct, ]
    k <- min(1000, dim(ct_marker_df)[1])
    ct_marker_df <- ct_marker_df[1:k,]
    sign_genes <- c(sign_genes, ct_marker_df$gene)
  }
  sign_genes <- unique(sign_genes)
  print(paste('Length of signature genes = ', length(sign_genes), sep=''))
  expr_matrix <- as.matrix(sc_obj@assays$RNA@data)
  sig_matrix <- makeSignMatrixDWLSfromMatrix(expr_matrix, sign_genes, cell_types)
  return(sig_matrix)
}

make_mixture_giotto <- function(mixture_count_fn){
  mixture_expr <- t(as.matrix(read.csv(mixture_count_fn, row.names = 1)))
  locations <- colnames(mixture_expr)
  xs <- c()
  ys <- c()
  for (spot in locations){
    x <- as.integer(strsplit(spot, split = "x")[[1]][1])
    y <- as.integer(strsplit(spot, split = "x")[[1]][2])
    xs <- c(xs, x)
    ys <- c(ys, y)
  }
  locations <- data.frame(row.names=locations, x=xs, y=ys)
  
  mixture_giotto <- createGiottoObject(mixture_expr, spatial_locs = locations)
  ## QC of giotto 
  mixture_giotto <- filterGiotto(gobject = mixture_giotto, 
                                 expression_threshold = 1, 
                                 gene_det_in_min_cells = 20, 
                                 min_det_genes_per_cell = 5)
  ## Normalize mixture
  mixture_giotto <- normalizeGiotto(gobject = mixture_giotto, scalefactor = 1e6, verbose = T)
  
  ## Cluster Giotto
  mixture_giotto <- calculateHVG(gobject = mixture_giotto)
  mixture_giotto <- runPCA(gobject = mixture_giotto)
  mixture_giotto <- runUMAP(mixture_giotto, dimensions_to_use = 1:5)
  mixture_giotto <- createNearestNetwork(gobject = mixture_giotto, dimensions_to_use = 1:30, k = 5)
  print("Network created. Now doing leiden clustering.")
  mixture_giotto =  doLeidenCluster(mixture_giotto, name = 'leiden_clus')
  return(mixture_giotto)
}

run_spatialDWLS <- function(mixture_count_fn, ref_count_fn, ref_meta_fn, col){
  mixture_giotto <- make_mixture_giotto(mixture_count_fn)
  sig_matrix <- make_signature_matrix(ref_count_fn, ref_meta_fn, col)
  start_time <- Sys.time()
  mixture_giotto<- runSpatialDeconv(mixture_giotto, deconv_method = c("DWLS"),
                                       expression_values = c("normalized"),logbase = 2,
                                       cluster_column = "leiden_clus", sig_matrix,
                                       n_cell = 10, cutoff = 2, return_gobject = F)
  end_time <- Sys.time()
  return(mixture_giotto)
}

# ### Helper functions end here
# 
# ## INPUT PARAMETERS
# mixture_fn <- file.path(workDir, 'data/cancer_inf_0.1/simulated_mixture_raw_counts.csv')
# 
# ### EXTERNAL REFERENCE CORRECTED
# ref_count_fn1 <- file.path(workDir, 'single_cell_ref/single_cell_ref_count_corrected_external.csv')
# ref_meta_fn1 <- file.path(workDir, 'single_cell_ref/single_cell_ref_meta_corrected_external.csv')
# 
# res1 <- run_spatialDWLS(mixture_fn, ref_count_fn1, ref_meta_fn1, 'bio_celltype')
# write.csv(res1, file.path(workDir, 'data/cancer_inf_0.1/spatialDWLS_minor_results/estimated_proportions_ref_external_corrected.csv'))
