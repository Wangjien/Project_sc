library(seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(SeuratWrappers)

# ----------------------------------------------------------------
# 使用fastmnn 的大群结果进行分分析
# ----------------------------------------------------------------

scRNAlist <- SplitObject(NKT_FastMNN,split.by="Patient")
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x))

scRNA <- RunFastMNN(object.list = scRNAlist)
scRNA <- RunUMAP(scRNA, reduction = "mnn", dims = 1:30)
scRNA <- RunTSNE(scRNA, reduction = "mnn", dims = 1:30)
scRNA <- FindNeighbors(scRNA, reduction = "mnn", dims = 1:30)
scRNA <- FindClusters(scRNA,resolution=seq(0.05,1,0.05))




plist = list()
for(i in names(flist)){
    plist = DimPlot(flist[i], raster = F) + labs(title = i)
}

