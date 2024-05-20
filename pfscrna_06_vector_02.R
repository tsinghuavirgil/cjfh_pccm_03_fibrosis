rm(list=ls())
#############
library(devtools)
library(BiocManager)
library(data.table)
library(circlize)#
library(gatepoints)#
library(stringr)
library(igraph)
library(gmodels)
library(Seurat)
library(reticulate)
library(SingleCellExperiment)
library(monocle)
library(slingshot)
##################
pbmc <- readRDS('cds_mes.rds')
pbmc
pbmc <- NormalizeData(pbmc,normalization.method="LogNormalize",
                      scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc,selection.method="vst",
                             nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),
               npcs=50)
pbmc <- RunUMAP(pbmc, dims = 1:50)
DimPlot(pbmc, reduction = "umap")
#####################################
#source('https://raw.githubusercontent.com/jumphone/Vector/master/Vector.R')
pbmc_umap <- as.matrix(read.table('mes_umapcord.tsv',sep='\t',
                                  header=T,row.names=1))
pbmc@reductions$umap@cell.embeddings<-pbmc_umap
pbmc_pca <- as.matrix(read.table('mes_pcacord.tsv',sep='\t',
                                 header=T,row.names=1))
pbmc@reductions$pca@cell.embeddings<-pbmc_pca
#######
vec_pca=vector.rankPCA(pbmc_pca)
pbmc_out=vector.buildGrid(pbmc_umap, N=30,SHOW=TRUE)
# Build network
pbmc_out=vector.buildNet(pbmc_out, CUT=1, SHOW=TRUE)

# Calculate Quantile Polarization (QP) score
pbmc_out=vector.getValue(pbmc_out, vec_pca, SHOW=TRUE)

# Get pixel's QP score
pbmc_out=vector.gridValue(pbmc_out,SHOW=TRUE)

# Find starting point
pbmc_out=vector.autoCenter(pbmc_out,UP=0.9,SHOW=TRUE)

# Infer vector
pbmc_out<-vector.drawArrow(pbmc_out,P=0.9,SHOW=TRUE,
                           COL=pbmc_out$COL,
                           SHOW.SUMMIT=TRUE)
#pbmc_out$P.PS : Peseudotime Score (PS) of each cell
pbmc_ps<-cbind(pbmc_out$P.PS,pbmc_out$USED_INDEX)
colnames(pbmc_ps) = c('pseudotime_score','old_index')
######
#gene#
######
target_gene = pbmc@assays$RNA@data[which(rownames(pbmc) =='COL1A1'),]
pbmc_out=vector.buildGrid(pbmc_umap, N=30,SHOW=TRUE)
pbmc_out=vector.buildNet(pbmc_out, CUT=1, SHOW=TRUE)
pbmc_out=vector.getValue(pbmc_out, vec_pca, SHOW=TRUE)

pbmc_out$VALUE=target_gene

pbmc_out=vector.showValue(pbmc_out)
pbmc_out=vector.gridValue(pbmc_out, SHOW=TRUE)
pbmc_out=vector.autoCenter(pbmc_out,UP=0.9,SHOW=TRUE)
pbmc_out=vector.drawArrow(pbmc_out,P=0.9,SHOW=TRUE, COL=pbmc_out$COL)
###############
pbmc_out=vector.buildGrid(pbmc_umap, N=30,SHOW=TRUE)
pbmc_out=vector.buildNet(pbmc_out, CUT=1, SHOW=TRUE)
pbmc_out=vector.getValue(pbmc_out, vec_pca, SHOW=TRUE)
pbmc_out=vector.gridValue(pbmc_out,SHOW=TRUE)

pbmc_out=vector.selectCenter(pbmc_out)

pbmc_out=vector.drawArrow(pbmc_out,P=0.9,SHOW=TRUE, COL=pbmc_out$COL)