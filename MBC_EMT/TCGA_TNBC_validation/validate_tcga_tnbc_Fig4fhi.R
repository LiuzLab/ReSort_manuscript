library(pheatmap)
library(Seurat)
library(ggplot2)
### Reproduce Figure 4f, h, i
setwd("/mnt/humble_50t/alexw/ReSort_manuscript/LingEMT_Dec22/TCGA_TNBC_validation/")
tnbrca.count <- read.table("nature11412-s2/TCGA_BRCA_TN_RPKM_libnorm.txt")

em_markers <- read.csv("top100_em_markers.csv")
em_markers <-em_markers[order(em_markers$marker),]

em_markers <- em_markers[,-1]
row.names(em_markers) <- em_markers$gene
View(em_markers)

genes <- row.names(em_markers)
length(genes)
genes <- genes[(genes %in% row.names(tnbrca.count))]
length(genes)


tnbrca.count <- tnbrca.count[genes,]
dim(tnbrca.count)
em_markers <- em_markers[row.names(tnbrca.count),]

pt = pheatmap(tnbrca.count, scale='row', 
              annotation_row  = data.frame(row.names = row.names(em_markers),
                                           marker = em_markers$marker), cluster_rows = F)
print(pt)
sample_clust <- factor(cutree(pt$tree_col, k=13))
pheatmap(tnbrca.count, scale='row',
         annotation_row  = data.frame(row.names = row.names(em_markers), marker = em_markers$marker),
         annotation_col = data.frame(row.names=names(sample_clust), clust=sample_clust),cluster_rows = F)

sample_clust_df <- ifelse(sample_clust==3,"Mlike", ifelse(sample_clust==6, "Elike", "Others"))
table(sample_clust_df)
sample_clust_df <- data.frame(row.names=names(sample_clust), 'EM'=sample_clust_df)

sample_cluster_df.1 <- data.frame(row.names = row.names(sample_clust_df)[sample_clust_df$EM %in% c("Elike", 'Mlike')],
                                  EM = sample_clust_df[sample_clust_df$EM %in% c("Elike", 'Mlike'),])
dim(sample_cluster_df.1)
tnbrca.count.1 <- tnbrca.count[,row.names(sample_cluster_df.1)]
tnbrca.count.2 <- t(scale(t(tnbrca.count.1)))
tnbrca.count.2[tnbrca.count.2 > 3] <- 3
tnbrca.count.2[tnbrca.count.2 < -3] <- 3

pheat <- pheatmap(tnbrca.count.2, scale='none', 
         annotation_row  = data.frame(row.names = row.names(em_markers),
                                      marker = em_markers$marker), 
         annotation_col = sample_cluster_df.1,
         cluster_rows = F,
         show_rownames = F)

pdf("heatmap_w_EM-like_annotations.pdf", width=7, height=6)
print(pheat)
dev.off()

cibersort.res <- read.csv("CIBERSORTx_TCGA_BRCA_TN_RPKM_libnormcsv.csv", row.names = 1)
rnames <- c()
for (rn in row.names(cibersort.res)){
  rnames <- c(rnames, gsub( '-', '.', rn))
}

row.names(cibersort.res) <- rnames
cibersort.res <- cibersort.res[, 1:22]
cibersort.res <- cibersort.res[, c('Macrophages.M2', 'Macrophages.M1', 'Macrophages.M0')]
cibersort.res <- cibersort.res[row.names(sample_cluster_df.1),]
#cibersort.res[cibersort.res<0] <- 0
#cibersort.res <- as.data.frame(t(apply(cibersort.res, 1, function(x)(x/sum(x)))))
cibersort.res['EM'] <- sample_cluster_df.1$EM
cibersort.res <- cibersort.res[complete.cases(cibersort.res),]


library(reshape2)
library(ggpubr)
cibersort.res.long <- cibersort.res
cibersort.res.long$identifier <- row.names(cibersort.res.long)
cibersort.res.long <- melt(cibersort.res.long, id.vars = c("identifier", "EM"))
cibersort.res.long$variable <- factor(cibersort.res.long$variable, levels=c('Macrophages.M1', 'Macrophages.M0', 'Macrophages.M2'))
View(cibersort.res.long)

ggplot(cibersort.res.long, aes(x=value, y=EM, fill=EM)) + 
  geom_violin(position=position_dodge(1)) + 
  stat_summary(fun=median, geom="point", size=3, color="black",position=position_dodge(1)) + 
  stat_compare_means(mapping=aes(x=value, y=EM, fill=EM), method = "t.test", size=5, 
                     method.args = list(alternative = "greater"),
                     label='p.format') + 
  labs(x= "Abundance", y = "Tumor") +
  # scale_x_discrete(labels= c('M1', 'M0', 'M2')) +
  theme_tufte(base_size = 20) + geom_rangeframe() + 
  facet_wrap(~variable,scales = "free", ncol = 1) 

t.test(as.numeric(cibersort.res[cibersort.res$EM == 'Elike', 'Macrophages.M0']),
       as.numeric(cibersort.res[cibersort.res$EM == 'Mlike', 'Macrophages.M0']))
# p=0.027
t.test(as.numeric(cibersort.res[cibersort.res$EM == 'Elike', 'Macrophages.M2']),
       as.numeric(cibersort.res[cibersort.res$EM == 'Mlike', 'Macrophages.M2']))
# p=0.0005
t.test(as.numeric(cibersort.res[cibersort.res$EM == 'Elike', 'Macrophages.M1']),
       as.numeric(cibersort.res[cibersort.res$EM == 'Mlike', 'Macrophages.M1']))
# p=0.069

mouse_res_fn <- "/mnt/humble_50t/alexw/ReSort_manuscript/LingEMT_Dec22/data/ReSort_secondary/ReSort_props.csv"
#ori_fn <- "/home/humble_local_25t/alexw/projects/ShawnLab/EMS_deconv_0307/scripts/s24679_LM22_props.csv"
mouse_res <- read.csv(mouse_res_fn, row.names=1)
mouse_res <- mouse_res[mouse_res$MixtureType %in% c("ES Mixture", "MS Mixture"), c('Macrophages.M1', 'Macrophages.M0', 'Macrophages.M2', 'MixtureType')]

mouse_res.long <- melt(mouse_res, id.vars = c("MixtureType"))

ggplot(mouse_res.long, aes(x=value, y=MixtureType, fill=MixtureType)) + 
  geom_violin(position=position_dodge(1)) +
  stat_summary(fun=median, geom="point", size=2, 
               color="black",position=position_dodge(1)) +
  facet_wrap(~variable,scales = "free_x") + 
  theme_tufte(base_size = 20) + geom_rangeframe() + 
  labs(x= "Cell type proportion (%)", y = "Tumor")  
  # theme_tufte(base_size = 20) + geom_rangeframe()

t.test(as.numeric(mouse_res[mouse_res$MixtureType == 'ES Mixture', 'Macrophages.M0']),
       as.numeric(mouse_res[mouse_res$MixtureType == 'MS Mixture', 'Macrophages.M0']),
       alternative='greater')
# p=3e-5

t.test(as.numeric(mouse_res[mouse_res$MixtureType == 'ES Mixture', 'Macrophages.M2']),
       as.numeric(mouse_res[mouse_res$MixtureType == 'MS Mixture', 'Macrophages.M2']),
       alternative='greater')
# p=0.0004

t.test(as.numeric(mouse_res[mouse_res$MixtureType == 'MS Mixture', 'Macrophages.M1']),
       as.numeric(mouse_res[mouse_res$MixtureType == 'ES Mixture', 'Macrophages.M1']),
       alternative='greater')
# p=0.0003