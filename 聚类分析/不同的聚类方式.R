library(Seurat)
library(SeuratWrapper)
library(dplyr)
library(stringr)
library(rliger)
library(harmony)
library(patchwork)
library(tidyverse)
library(conos)


Run_FastMNN = function(sce, batch, nfeatures){
    if(Reduce('!', is(sce)=='Seurat')){
        scRNAlist = SplitObject(scRNA, split.by = batch)
        scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))
        scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x, nfeatures=nfeatures))
        scRNA <- RunFastMNN(object.list = scRNAlist)
        scRNA <- RunUMAP(scRNA, reduction = "mnn", dims = 1:30)
        scRNA <- FindNeighbors(scRNA, reduction = "mnn", dims = 1:30)
        scRNA <- FindClusters(scRNA)
        return(scRNA)
    }

}

#---------------------多个nfeatures参数结果汇总--------------------------------
Run_plot_DimPlot = function(sce, groupBy=NULL, title = NULL){
    if(is.null(groupBy)){
        tmp = DimPlot(sce, label=T, repel=T, raster=F) + 
            labs(title = title) +
            theme(plot.title = element_text(hjust = 0.5)) 
    }else{
        tmp = DimPlot(sce, label=T, repel=T, raster=F, group.by=groupBy) +
            labs(title = title) +
            theme(plot.title = element_text(hjust = 0.5))
    }
    return(tmp)
}

p0 = Run_plot_DimPlot(sce = NKT_FastMNN_500, title = "HVG:500")
p1 = Run_plot_DimPlot(sce = NKT_FastMNN, title = "HVG:1000")
p2 = Run_plot_DimPlot(sce = NKT_FastMNN_1500, title = "HVG:1500")
p3 = Run_plot_DimPlot(sce = NKT_FastMNN_10000, title = "HVG:10000")

p4 = Run_plot_DimPlot(sce = NKT_FastMNN, title = "HVG:1000", groupBy = 'Response')
p5 = Run_plot_DimPlot(sce = NKT_FastMNN_1500, title = "HVG:1500", groupBy = 'Response')
p6 = Run_plot_DimPlot(sce = NKT_FastMNN_10000, title = "HVG:10000", groupBy = 'Response')

png('./NKT_FastMNN_Diff_HVG.png', height = 2000, width = 8000, res=300)
p0|p1 | p2 | p3
dev.off()

# ------------------ 添加信息后再绘制DimPlot------------------------------------

