library(ggplot2)
library(Seurat)
library(dplyr)
library(dior)
library(tdiyverse)
library(data.table)
library(stringr)
library(futile.logger)


load('/root/wangje/Project/刘老师/大群/Data/大群CCA.Rdata')
scRNA = CreateSeuratObject(scRNA_seurat@assays$RNA@counts)
scRNA_seurat = NULL
gc()
scRNA$sample = str_split_fixed(rownames(scRNA@meta.data), '_[A|T|G|C].*', n = 2)[,1]
scRNA$percent.mt = PercentageFeatureSet(scRNA, pattern = '^MT-')
dior::write_h5(scRNA, file = "/root/wangje/Project/刘老师/大群/Data/bigCluster.h5")

# 每个样本单独聚类
sample_list = list()
for(index in unique(scRNA$sample)){
    print(index)
    suppressMessages({
    tmp = scRNA[,scRNA$sample == index]
    sample_list[[index]] = CreateSeuratObject(tmp@assays$RNA@counts)})
} 

# 进行聚类
flist = list()
plist1 = list()
for(i in 1:length(sample_list)){
    print(names(sample_list[i]))
    suppressMessages({
        sample_list[[i]]$percent.mt = PercentageFeatureSet(sample_list[[i]], pattern = '^MT-')
        flog.info("开始绘制小提琴图！")
        plist1[[names(sample_list[i])]] = VlnPlot(sample_list[[i]], features = c('nCount_RNA',
                                                    'nFeature_RNA', 'percent.mt'), ncol = 3) + 
                                                    labs(subtitle = names(sample_list[i]))
        flog.info("绘制小提琴结束！")
        tmp = NormalizeData(sample_list[[i]])
        tmp = FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
        tmp = ScaleData(tmp)
        tmp = RunPCA(tmp)
        tmp <- FindNeighbors(tmp, dims = 1:30)
        tmp <- FindClusters(tmp, resolution = 0.5)
        tmp <- RunUMAP(tmp, dims = 1:30)          
        flist[[names(sample_list)[i]]] = tmp      
    })
}

# 绘制DimPlot
plist2 = list()
for(i in 1:length(flist)){
    flog.info(names(flist[i]))
    suppressMessages({
        plist2[[names(flist[i])]] = DimPlot(flist[[i]], label=T, repel = T, raster=F)
    })
}