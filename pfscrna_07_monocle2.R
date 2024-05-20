rm(list=ls())
#############
library(remotes)
library(pheatmap)
library(dplyr)
library(monocle)
library(clusterProfiler)
library(Seurat)
library(reticulate)
library(sceasy)
library(SingleCellExperiment)
library(tidyverse)
library(patchwork)
library(reshape2)
library(ggpubr)
library(colorRamps)
###################
###################
seurat_mes <- readRDS('seurat_mes.rds')
seurat_mes[['percent.mt']] <- PercentageFeatureSet(seurat_mes, pattern = '^MT-')
seurat_mes <- NormalizeData(seurat_mes, normalization.method = 'LogNormalize', scale.factor = 10000)
seurat_mes <- FindVariableFeatures(seurat_mes, selection.method = 'vst', nfeatures = 2000)
all_genes <- rownames(seurat_mes)
seurat_mes <- ScaleData(seurat_mes, features = all_genes)
#########################################################
cds_mes <- readRDS('cds_01_mes.rds')
cds_mes <- estimateSizeFactors(cds_mes)
cds_mes
cds_mes <- estimateDispersions(cds_mes, fitType='locfit')#
cds_mes
colnames(pData(cds_mes))
cds_mes@phenoData@data[['Cluster']] <- cds_mes@phenoData@data[['celltype']]
table(seurat_mes$celltype)
Idents(seurat_mes)='celltype'
deg_cluster <- FindAllMarkers(seurat_mes)
table(deg_cluster$cluster)
express_genes <- subset(deg_cluster,p_val_adj<0.05)$gene

cds_mes <- setOrderingFilter(cds_mes,express_genes)
cds_mes
plot_ordering_genes(cds_mes)
colnames(pData(cds_mes))
#############
cds_mes <- reduceDimension(cds_mes,max_components = 2,method = 'DDRTree')
colnames(pData(cds_mes))
cds_mes <- orderCells(cds_mes)
colnames(pData(cds_mes))
pData(cds_mes)$State <- factor(pData(cds_mes)$Cluster)
colnames(cds_mes)
pseudo_data <- as.data.frame(cds_mes@phenoData@data[['Pseudotime']],row.names = colnames(cds_mes))
pseudo_pheno_data <- cds_mes@phenoData@data
write.csv(pseudo_pheno_data,'pseudotime_mes.csv')
plot_cell_trajectory(cds_mes,color_by='Pseudotime',size=1,show_backbone=TRUE)
plot_cell_trajectory(cds_mes,color_by='Cluster',size=1,show_backbone=TRUE)
plot_cell_trajectory(cds_mes, color_by = 'State', size=1,show_backbone=TRUE)
plot_cell_trajectory(cds_mes, color_by = 'State') + facet_wrap('~State', nrow = 1)
plot_complex_cell_trajectory(cds_mes, x = 1, y = 2,color_by = 'Cluster')
############################
pheno_data <- pData(cds_mes)

ggplot(pheno_data, aes(Pseudotime, colour = Cluster, fill=Cluster)) +
  
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
#############################
ggplot(pheno_data, aes(Pseudotime, colour =Cluster, fill=Cluster)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+
  theme_classic2()+
  scale_color_manual(name= '', values = ClusterName_color_panel)