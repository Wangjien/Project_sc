library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(SeuratWrappers)

# ----------------------------------------------------------------
# 使用fastmnn 的大群结果进行分分析
# ----------------------------------------------------------------
# load data 
big.fastmnn = '/root/wangje/wangjien_new/liu-big/scRNA_FastMNN_new.RData'
load(big.fastmnn)
scRNA_sub = scRNA[,scRNA$celltype == 'Epithelials']

scRNAlist <- SplitObject(scRNA_sub,split.by="sample")
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x))

scRNA <- RunFastMNN(object.list = scRNAlist)
scRNA <- RunUMAP(scRNA, reduction = "mnn", dims = 1:50)
scRNA <- RunTSNE(scRNA, reduction = "mnn", dims = 1:50)
scRNA <- FindNeighbors(scRNA, reduction = "mnn", dims = 1:50)
scRNA <- FindClusters(scRNA,resolution=seq(0.05,1,0.05))



plist = list()
for(i in names(flist)){
    plist = DimPlot(flist[i], raster = F) + labs(title = i)
}

