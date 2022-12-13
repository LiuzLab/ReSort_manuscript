workDir <- "/houston_20t/alexw/ST/iSort/simulation_finer_resolution/scripts/run_models/"
setwd(workDir)

source('run_MuSiC.R')
source('run_RCTD.R')
source('runSPOTlight.R')
source('run_Giotto_spatialDWLS.R')

mixture_fn <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/data/cancer_inf_0.1/simulated_mixture_raw_counts.csv'
res_folder <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/results/'
dir.create(res_folder, 'None')
  
#### 1. Reference internal ####
inter.count.fn <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/ref_cp10k/count_internal.csv'
inter.meta.fn <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/ref_cp10k/meta_internal.csv'

#### 2. Reference external ####
exter.count.fn <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/ref_cp10k/count_external.csv'
exter.meta.fn <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/ref_cp10k/meta_external.csv'

#### 3. Reference MIST ####
mist.count.fn <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/ref_cp10k/count_MIST.csv'
mist.meta.fn <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/ref_cp10k/meta_MIST.csv'

#### 4. Reference Immune ####
immune.count.fn <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/ref_cp10k/count_external_immune.csv'
immune.meta.fn <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/ref_cp10k/meta_external_immune.csv'

#### 5. Reference Immune Corrected ####
ic.count.fn <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/ref_cp10k/count_external_immune_corrected.csv'
ic.meta.fn <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/ref_cp10k/meta_external_immune.csv'

#### 6. Reference All corrected ####
ac.count.fn <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/ref_cp10k/count_external_corrected.csv'
ac.meta.fn <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/ref_cp10k/meta_external.csv'


########################################### 1. RCTD ###########################################
# mixture.RCTD <- make_mixture_RCTD(mixutre_fn)
# RCTD.result.folder <- file.path(res_folder, 'RCTD')
# dir.create(RCTD.result.folder, 'None')
# 
# ## 1
# rctd.ref1 <- createRef(inter.count.fn, inter.meta.fn)
# rctd.props1 <- run_RCTD(mixture.RCTD, rctd.ref1)
# write.csv(rctd.props1, file.path(RCTD.result.folder, 'props_internal.csv'))
# 
# ## 2
# rctd.ref2 <- createRef(exter.count.fn, exter.meta.fn)
# rctd.props2 <- run_RCTD(mixture.RCTD, rctd.ref2)
# write.csv(rctd.props2, file.path(RCTD.result.folder, 'props_external.csv'))
# 
# ## 3
# rctd.ref3 <- createRef(mist.count.fn, mist.meta.fn)
# rctd.props3 <- run_RCTD(mixture.RCTD, rctd.ref3)
# write.csv(rctd.props3, file.path(RCTD.result.folder, 'props_MIST.csv'))
# 
# ## 4
# rctd.ref4 <- createRef(immune.count.fn, immune.meta.fn)
# rctd.props4 <- run_RCTD(mixture.RCTD, rctd.ref4)
# write.csv(rctd.props4, file.path(RCTD.result.folder, 'props_immune.csv'))
# 
# ## 5
# rctd.ref5 <- createRef(ic.count.fn, ic.meta.fn)
# rctd.props5 <- run_RCTD(mixture.RCTD, rctd.ref5)
# write.csv(rctd.props5, file.path(RCTD.result.folder, 'props_immune_corrected.csv'))
# 
# ## 6
# rctd.ref6 <- createRef(ac.count.fn, ac.meta.fn)
# rctd.props6 <- run_RCTD(mixture.RCTD, rctd.ref6)
# write.csv(rctd.props6, file.path(RCTD.result.folder, 'props_all_corrected.csv'))

########################################### 2. MUSIC ###########################################
# music.result.folder <- file.path(res_folder, 'music')
# dir.create(music.result.folder, 'None')
# 
# ## 1
# music.props1 <- runMusic(mixutre_fn, inter.count.fn, inter.meta.fn, 'bio_celltype')
# write.csv(music.props1, file.path(music.result.folder, 'props_internal.csv'))
# 
# ## 2
# music.props2 <- runMusic(mixutre_fn, exter.count.fn, exter.meta.fn, 'bio_celltype')
# write.csv(music.props2, file.path(music.result.folder, 'props_external.csv'))
# 
# ## 3
# music.props3 <- runMusic(mixutre_fn, mist.count.fn, mist.meta.fn, 'bio_celltype')
# write.csv(music.props3, file.path(music.result.folder, 'props_MIST.csv'))
# 
# ## 4
# music.props4 <- runMusic(mixutre_fn, immune.count.fn, immune.meta.fn, 'bio_celltype')
# write.csv(music.props4, file.path(music.result.folder, 'props_immune.csv'))
# 
# ## 5
# music.props5 <- runMusic(mixutre_fn, ic.count.fn, ic.meta.fn, 'bio_celltype')
# write.csv(music.props5, file.path(music.result.folder, 'props_immune_corrected.csv'))
# 
# ## 6
# music.props6 <- runMusic(mixutre_fn, ac.count.fn, ac.meta.fn, 'bio_celltype')
# write.csv(music.props6, file.path(music.result.folder, 'props_all_corrected.csv'))
# 

########################################### 3. spatialDWLS ###########################################
spatialDWLS.result.folder <- file.path(res_folder, 'spatialDWLS')
dir.create(spatialDWLS.result.folder, 'None')

## 1
spatialDWLS.props1 <- run_spatialDWLS(mixture_fn, inter.count.fn, inter.meta.fn, 'bio_celltype')
write.csv(spatialDWLS.props1, file.path(spatialDWLS.result.folder, 'props_internal.csv'))

## 2
spatialDWLS.props2 <- run_spatialDWLS(mixture_fn, exter.count.fn, exter.meta.fn, 'bio_celltype')
write.csv(spatialDWLS.props2, file.path(spatialDWLS.result.folder, 'props_external.csv'))

## 3
spatialDWLS.props3 <- run_spatialDWLS(mixture_fn, mist.count.fn, mist.meta.fn, 'bio_celltype')
write.csv(spatialDWLS.props3, file.path(spatialDWLS.result.folder, 'props_mist.csv'))

## 4
spatialDWLS.props4 <- run_spatialDWLS(mixture_fn, immune.count.fn, immune.meta.fn, 'bio_celltype')
write.csv(spatialDWLS.props4, file.path(spatialDWLS.result.folder, 'props_immune.csv'))

## 5
spatialDWLS.props5 <- run_spatialDWLS(mixture_fn, ic.count.fn, icmeta.fn, 'bio_celltype')
write.csv(spatialDWLS.props5, file.path(spatialDWLS.result.folder, 'props_immune_corrected.csv'))

## 6
spatialDWLS.props6 <- run_spatialDWLS(mixture_fn, ac.count.fn, ac.meta.fn, 'bio_celltype')
write.csv(spatialDWLS.props6, file.path(spatialDWLS.result.folder, 'props_all_corrected.csv'))

########################################### 4. spotLight ###########################################
# spotlight.result.folder <- file.path(res_folder, 'spotlight')
# dir.create(spotlight.result.folder, 'None')
# ## 1
# spotlight.props1 <- run_spotlight(mixture_fn, inter.count.fn, inter.meta.fn, 'bio_celltype')
# write.csv(spotlight.props1, file.path(spotlight.result.folder, 'props_internal.csv'))
# ## 2
# spotlight.props2 <- run_spotlight(mixture_fn, exter.count.fn, exter.meta.fn, 'bio_celltype')
# write.csv(spotlight.props2, file.path(spotlight.result.folder, 'props_external.csv'))
# ## 3
# spotlight.props3 <- run_spotlight(mixture_fn, mist.count.fn, mist.meta.fn, 'bio_celltype')
# write.csv(spotlight.props3, file.path(spotlight.result.folder, 'props_MIST.csv'))
# ## 4
# spotlight.props4 <- run_spotlight(mixture_fn, immune.count.fn, immune.meta.fn, 'bio_celltype')
# write.csv(spotlight.props4, file.path(spotlight.result.folder, 'props_immune.csv'))
# ## 5
# spotlight.props5 <- run_spotlight(mixture_fn, ic.count.fn, ic.meta.fn, 'bio_celltype')
# write.csv(spotlight.props5, file.path(spotlight.result.folder, 'props_immune_corrected.csv'))
# ## 6
# spotlight.props6 <- run_spotlight(mixture_fn, ac.count.fn, ac.meta.fn, 'bio_celltype')
# write.csv(spotlight.props6, file.path(spotlight.result.folder, 'props_all_corrected.csv'))
