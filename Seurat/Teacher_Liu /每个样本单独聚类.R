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


#----------------------------------------------------------------
# 查看T细胞中的其他亚群marker
Tcells <- c('CD3D','CD3G','CD2')
NK_cells <- c('KLRC1','KLRC3','TRDV9')
Fibroblasts <- c('COL1A1','DCN','LUM')
Myeloids <- c('LYZ','CD68','TYROBP')
Epithelial <- c('CD24','KRT19','EPCAM')
Bcells <- c('CD79A','CD19','MS4A1')
Endothelial <- c('CLDN5','FLT1','RAMP2')
Plasma <- c('IGHG1','JCHAIN','MZB1')
Hepatocytes <- c('ALB','APOB','HP')
Keratinocytes <- c("KRT5","KRT14","FABP5")
DC <- c("LILRA4","CXCR3","IRF7")
Mast <- c("CPA3","TPSABT","TPSB2")

marker.list <- list("NK cell"=NK_cells,"T cell"=Tcells,'B cell'=Bcells,
                    "Plasmas" =Plasma,Myeloids=Myeloids,Fibroblasts=Fibroblasts,
                    Epithelials=Epithelial,Endothelials=Endothelial,Hepatocytes=Hepatocytes,
                    Keratinocytes=Keratinocytes,DC = DC, Mast = Mast)
                    
plotFeature <- function(scRNA_data = scRNA_data,
                        choose = "Feature",
                        col_num = 6, marker.list = marker.list,...) {
    pacman::p_load("Seurat", "ggplot2", "tidyverse")
    DefaultAssay(scRNA_data) <- "RNA"
    plist <- list()
    if (is.null(choose)) {
        message("请选择绘图类型")
    } else if (choose == "Feature") {
        for (i in names(marker.list)) {
            for (j in marker.list[[i]]) {
                #    print(paste0(i,"_",j))
                tmp <- tryCatch(
                    {
                        FeaturePlot(scRNA_data, features = j,raster=F) +
                            theme(legend.position = "right") +
                            labs(title = paste0(i, "_", j))
                    },
                    error = function(e) {
                        message("Error @ ", j)
                        return(NA)
                    },
                    finally = {
                        message(paste0(i, "_", j, "_next..."))
                    }
                )
                plist[[paste0(i, "_", j)]] <- tmp
            }
        }
        p_new <- Filter(Negate(anyNA), plist)
        p <- wrap_plots(p_new, bycol = T, ncol = col_num)
        return(p)
    } else if (choose == "SCpubr") {
        for (i in names(marker.list)) {
            for (j in marker.list[[i]]) {
                #    print(paste0(i,"_",j))
                tmp <- tryCatch(
                    {
                        SCpubr::do_NebulosaPlot(scRNA_data, features = j, viridis_color_map = "H", pt.size = 0.02) +
                            theme(legend.position = "righty") +
                            labs(title = paste0(i, "_", j))
                    },
                    error = function(e) {
                        message("Error @ ", j)
                        return(NA)
                    },
                    finally = {
                        message(paste0(i, "_", j, "_next..."))
                    }
                )
                plist[[paste0(i, "_", j)]] <- tmp
            }
        }
        p_new <- Filter(Negate(anyNA), plist)
        p <- wrap_plots(p_new, bycol = T, ncol = col_num)
        return(p)
    }
    ...
}

# 循环绘图
# 绘制每个样本的大群marker 基因
plist = list()
for(index in 1:length(flist)){
    flog.info(names(flist[index]))
    tmp = plotFeature(flist[[index]],choose = "Feature",col_num = 6, marker.list = marker.list )
    plist[[names(flist[index])]] = tmp
}

for(i in 1:l)