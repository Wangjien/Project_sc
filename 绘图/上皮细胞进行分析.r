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
scRNA = NULL

scRNAlist <- SplitObject(scRNA_sub,split.by="sample")
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x))

scRNA <- RunFastMNN(object.list = scRNAlist)
scRNA <- RunUMAP(scRNA, reduction = "mnn", dims = 1:50)
scRNA <- RunTSNE(scRNA, reduction = "mnn", dims = 1:50)
scRNA <- FindNeighbors(scRNA, reduction = "mnn", dims = 1:50)
scRNA <- FindClusters(scRNA,resolution=seq(0.05,1,0.05))

## 使用不同的min.dist
flist = list()
plist = list()
for(index in c(0.001,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)){
    cat('>>>>>>>> 此时的min.dist数值\n',index)
    scRNA <- FindNeighbors(scRNA, reduction = "mnn", dims = 1:50)
    scRNA <- FindClusters(scRNA,resolution=0.5)
    scRNA = RunUMAP(scRNA, reduction = "mnn", min.dist = index, dims = 1:50)
    flist[[as.character(index)]] = scRNA
    tmp = DimPlot(scRNA, raster=F) + labs(title = as.character(index)) + theme(plot.title = element_text(hjust = 0.5))
    plist[[as.character(index)]] = tmp
}

#------------------------------------------------------------------
# liger
#------------------------------------------------------------------
library(liger)



plist = list()
for(i in names(flist)){
    plist = DimPlot(flist[i], raster = F) + labs(title = i)
}

