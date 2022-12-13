#setwd("/houston_20t//alexw/ST/iSort/simulation_finer_resolution/")
# install the MuSiC package
# install.packages('devtools')
# library(devtools)
# devtools::install_github('xuranw/MuSiC')
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biobase")
library(MuSiC)
library(ggplot2)
library(ggpubr)
library(Biobase)

runMusic <- function(mixture_count_fn, ref_count_fn, ref_meta_fn, col){
  mixture = read.csv(mixture_count_fn, row.names=1)
  mixture = t(mixture)
  mixture.eset <- ExpressionSet(assayData = mixture)
  
  ref_count <-  t(as.matrix(read.csv(ref_count_fn, row.names=1, 
                                      check.names = F, stringsAsFactors = F)))
  ref_meta <- read.csv(ref_meta_fn, row.names=1, 
                        check.names = F, stringsAsFactors = F)
  
  
  ref_meta$spotID <- row.names(ref_meta)
  ref_meta <- ref_meta[, c(col, 'spotID')]
  ref_meta.meta <- data.frame(labelDescription=c(col, "spotID"),
                               row.names=c(col, "spotID"))
  
  phenoData <- new("AnnotatedDataFrame",data=ref_meta, varMetadata=ref_meta.meta)
  
  ref.eset <- ExpressionSet(assayData = ref_count, phenoData = phenoData)
  Est.prop.ref = music_prop(bulk.eset = mixture.eset, sc.eset = ref.eset,
                             clusters = col, samples = 'spotID',
                             select.ct = NULL, 
                             verbose = T)
  props <- as.data.frame(Est.prop.ref$Est.prop.weighted)
  return (props)
}

runMusic_project <- function(projDir, outDir, col='bio_celltype'){
  ## Create folder to save results if not existed
  # if (col == 'bio_celltype'){
  #   outDir <- file.path(projDir, 'music_minor_results')
  # } else {
  #   outDir <- file.path(projDir, 'music_major_results')
  # }
  
  dir.create(outDir, showWarnings = FALSE)
  
  mixture_count_fn <- file.path(projDir, 'simulated_mixture_raw_counts.csv')
  # refM_count_fn <- file.path(projDir, 'MIST_reference_count.csv')
  # refM_meta_fn <- file.path(projDir, 'MIST_reference_meta.csv')
  
  curwd <- getwd()
  
  # refA_count_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_count_internal.csv')
  # refA_meta_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_meta_internal.csv')
  # 
  # refB_count_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_count_external.csv')
  # refB_meta_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_meta_external.csv')
  # 
  # refI_count_fn <- file.path(curwd, 'single_cell_ref/immune_only_ref_count_external.csv')
  # refI_meta_fn <- file.path(curwd, 'single_cell_ref/immune_only_ref_meta_external.csv')
  # 
  refH_count_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_count_corrected_external.csv')
  refH_meta_fn <- file.path(curwd, 'single_cell_ref/single_cell_ref_meta_corrected_external.csv')
  
  # refI_count_fn <- file.path(curwd, 'single_cell_ref/immune_corrected_external_ref_count.csv')
  # refI_meta_fn <- file.path(curwd, 'single_cell_ref/immune_corrected_external_ref_meta.csv')
  
  
  
  ## Read and process mixtures
  mixture = read.csv(mixture_count_fn,row.names=1)
  mixture = t(mixture)
  mixture.eset <- ExpressionSet(assayData = mixture)

  # # PDACA as reference
  # st = Sys.time()
  # props.refA <- runMusic(mixture.eset = mixture.eset, ref_count_fn = refA_count_fn, ref_meta_fn = refA_meta_fn, col = col)
  # et = Sys.time()
  # print(et-st)
  # write.csv(props.refA, file.path(outDir, 'estimated_proportions_ref_internal.csv'))
  # print('Internal reference done.')
  #   
  # # PDACB as reference
  # st = Sys.time()
  # props.refB <- runMusic(mixture.eset = mixture.eset, ref_count_fn = refB_count_fn, ref_meta_fn = refB_meta_fn, col=col)
  # et = Sys.time()
  # print(et-st)
  # write.csv(props.refB, file.path(outDir, 'estimated_proportions_ref_external.csv'))
  # print('External reference done.')
  # 
  # # MIST as reference
  # if (col == 'bio_celltype_major'){
  #   st = Sys.time()
  #   props.refM <- runMusic(mixture.eset = mixture.eset, ref_count_fn = refM_count_fn, ref_meta_fn = refM_meta_fn, col=col)
  #   et = Sys.time()
  #   print(et-st)
  #   write.csv(props.refM, paste(outDir, 'estimated_proportions_ref_MIST.csv'))
  #   print('MIST reference done.')
  # }
  
  # Immune only as reference
  st = Sys.time()
  props.refH <- runMusic(mixture.eset = mixture.eset, ref_count_fn = refH_count_fn, ref_meta_fn = refH_meta_fn, col='bio_celltype')
  et = Sys.time()
  print(et-st)
  write.csv(props.refH, file.path(outDir, 'estimated_proportions_external_corrected.csv'))
  print('Immune only reference done.')
  
  # # Hybrid reference
  # st = Sys.time()
  # props.refH <- runMusic(mixture.eset = mixture.eset, ref_count_fn = refH_count_fn, ref_meta_fn = refH_meta_fn, col='bio_celltype')
  # et = Sys.time()
  # print(et-st)
  # write.csv(props.refH, file.path(projDir, 'music_estimated_proportions_external_mist_hybrid_external_immune.csv'))
  # print('Hybrid reference done.')
}

# ### INPUT FILES, TO BE CHANGED ####
# projDir <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/data/cancer_inf_0.1/'
# outDir <- file.path(projDir, 'music_minor_results/')
# workDir <- getwd()
# ref_count_fn <- file.path(workDir, 'single_cell_ref/single_cell_ref_count_corrected_augmented_external.csv')
# ref_meta_fn <- file.path(workDir, 'single_cell_ref/single_cell_ref_meta_corrected_augmented_external.csv')
# mixture_count_fn <- file.path(projDir, 'simulated_mixture_raw_counts.csv')
# mixture.props <- runMusic(mixture_count_fn, ref_count_fn, ref_meta_fn, 'bio_celltype')
# write.csv(mixture.props, file.path(projDir, 'music_minor_results/estimated_proportions_external_corrected_augmented.csv'))
# # system.time(runMusic_project(projDir, outDir, col = 'bio_celltype'))
# # system.time(runMusic_project(projDir, col = 'bio_celltype_major'))
