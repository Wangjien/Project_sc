# load
rm(list = ls())
gc()
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(future)
library(SCpubr)
library(SeuratWrappers)

setwd("/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts")
# future
plan(workers=60)
options(future.globals.maxSize= 100*1024^3)

############################################################
# Seurat CCA
############################################################

# load data
# ## -- Endothelials
# load("/root/wangje/Project/刘老师/Endothelials/Liu_Endothelials_seurat.RData")
# Endo <- scRNA_seurat
# rm(scRNA_seurat)
# ## -- CAFs
# load("/root/wangje/Project/刘老师/Fibroblasts/fibroblasts_cca.RData")
# CAF <- scRNA_seurat
# rm(scRNA_seurat)

# merge Endothelials and CAFs
# scRNA <- merge(Endo, CAF)
# rownames(Endo@meta.data) %in% rownames(CAF@meta.data)
load("/root/wangje/wangjien_new/liu-big/scRNA_FastMNN_new.RData")
# select Endothelials and CAFs from big cluster
CAF_Endo <- scRNA[, scRNA$celltype %in% c("Fibroblasts", "Endothelials")]
CAF_Endo <- CreateSeuratObject(CAF_Endo@assays$RNA@counts)

## add infomation
CAF_Endo$percent.mt <- PercentageFeatureSet(CAF_Endo, pattern = "^MT-")
CAF_Endo$Patient <- str_split_fixed(rownames(CAF_Endo@meta.data), "_[A|T|G|C],*", n = 2)[, 1]

## Seurat CCA
scRNAlist <- SplitObject(CAF_Endo, split.by = "Patient")
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x))
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)
scRNA_seurat <- IntegrateData(anchorset = scRNA.anchors)
scRNA_seurat <- ScaleData(scRNA_seurat, vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA")) %>% RunPCA(verbose = FALSE)
scRNA_seurat <- RunUMAP(scRNA_seurat, dims = 1:30)
scRNA_seurat <- FindNeighbors(scRNA_seurat, dims = 1:30) %>% FindClusters(resolution = seq(0.05, 2, 0.05))

## 查看其他细胞类型
# 查看时候混杂了其他的细胞
NKT <- c("CD3D", "CD3G", "CD2")
Fibroblasts <- c("COL1A1", "DCN", "LUM")
Myeloids <- c("LYZ", "CD68", "TYROBP")
Epithelial <- c("CD24", "KRT19", "EPCAM")
Bcells <- c("CD79A", "CD19", "MS4A1")
Endothelial <- c("CLDN5", "FLT1", "RAMP2")
Plasma <- c("IGHG1", "JCHAIN", "MZB1")
Hepatocytes <- c("ALB", "APOB", "HP")
Keratinocytes <- c("KRT5","KRT14","FABP5")
DC <- c('LILRA4','CXCR3','IRF7')
Mast <- c('CPA3','TPSAB1','TPSB2')
osteoblast <- c('COL2A1','SOX9','BMP2')

marker.list <- list(
    "NK&T cell" = NKT, "B cells" = Bcells, Plasmas = Plasma, Myeloids = Myeloids,
    Epithelials = Epithelial, Endothelials = Endothelial, Fibroblasts = Fibroblasts, Hepatocytes = Hepatocytes,
    Keratinocytes = Keratinocytes, DC = DC, Mast = Mast,
    osteoblast = osteoblast
)
# ----------------------- 写成函数 ---------------------------
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
                        FeaturePlot(scRNA_data, features = j) +
                            theme(legend.position = "none") +
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
                            theme(legend.position = "none") +
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
png("./查看是否有其他的细胞_FeaturePlot.png", height = 5000, width = 6000, res = 300)
plotFeature(scRNA_data = scRNA_seurat, choose = "Feature", col_num = 6, marker.list = marker.list)
dev.off()

png("./查看是否有其他的细胞_SCpubr.png", height = 4000, width = 6000, res = 300)
plotFeature(scRNA_data = scRNA_seurat, choose = "SCpubr", col_num = 6, marker.list = marker.list)
dev.off()

## 绘制DimPlot分群去除其他类型的细胞
png("./合并Fibroblasts和Endothelials_cca_umap.png", height = 2000, width = 2000, res = 300)
DimPlot(scRNA_seurat, label = T, repel = T)
dev.off()

## 去除其他类型的细胞（res为2）
scRNA_seurat <- scRNA_seurat[, !scRNA_seurat$integrated_snn_res.2 %in% c(26, 22, 23, 19, 25, 26, 12, 27)]

## 查看筛选后的数据的marker基因
png("/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/查看是否有其他的细胞_FeaturePlot_Sub.png", height = 4000, width = 6000, res = 300)
plotFeature(scRNA_data = scRNA_seurat, choose = "Feature", col_num = 6, marker.list = marker.list)
dev.off()

png("/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/查看是否有其他的细胞_SCpubr_Sub.png", height = 4000, width = 6000, res = 300)
plotFeature(scRNA_data = scRNA_seurat, choose = "SCpubr", col_num = 6, marker.list = marker.list)
dev.off()

## 再次使用Seurat CCA 进行分析
load("/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/CCA/合并Endo和CAFs_CCA_去除部分细胞.RData")
scRNA <- CreateSeuratObject(scRNA_seurat@assays$RNA@counts)
scRNA$Patient <- str_split_fixed(rownames(scRNA@meta.data),"_[A|T|G|C]",n=2)[,1]
scRNA$percent.mt <- PercentageFeatureSet(scRNA, pattern="MT-")
scRNAlist <- SplitObject(scRNA, split.by = "Patient")
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x))

scRNAlist <- SplitObject(scRNA_seurat, split.by = "Patient")
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)
scRNA_seurat <- IntegrateData(anchorset = scRNA.anchors)
scRNA_seurat <- ScaleData(scRNA_seurat, vars.to.regress = c("percent.mt")) %>% RunPCA(verbose = FALSE)
scRNA_seurat <- RunUMAP(scRNA_seurat, dims = 1:30)
scRNA_seurat <- FindNeighbors(scRNA_seurat, dims = 1:30) %>% FindClusters(resolution = seq(0.05, 2, 0.05))
save(scRNA_seurat,file="/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/CCA/合并Endo和CAFs_CCA_去除部分细胞_重新聚类.RData")

## FeaturePlot
setwd("/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/CCA")
png("/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/CCA/重新聚类后查看是否有其他的细胞_FeaturePlot.png", height = 6000, width = 6000, res = 300)
plotFeature(scRNA_data = scRNA_seurat, choose = "Feature", col_num = 6, marker.list = marker.list)
dev.off()

png("/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/CCA/重新聚类后查看是否有其他的细胞_SCpubr.png", height = 6000, width = 6000, res = 300)
plotFeature(scRNA_data = scRNA_seurat, choose = "SCpubr", col_num = 6, marker.list = marker.list)
dev.off()

##############################################################
# FastMNN
##############################################################
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts)
scRNA$Patient <- str_split_fixed(rownames(scRNA@meta.data),"_[A|T|G|C]",n=2)[,1]
scRNA$percent.mt <- PercentageFeatureSet(scRNA, pattern="MT-")
scRNAlist <- SplitObject(scRNA, split.by = "Patient")
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x))


scRNA <- RunFastMNN(object.list = scRNAlist)
scRNA <- RunUMAP(scRNA, reduction = "mnn", dims = 1:30)
scRNA <- RunTSNE(scRNA, reduction = "mnn", dims = 1:30)
scRNA <- FindNeighbors(scRNA, reduction = "mnn", dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = seq(0.05, 1, 0.05))
save(scRNA, file = "/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/FastMNN/去除部分细胞_合并endo和cafs_Fastmnn.RData")
# FeaturePlot
png("/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/FastMNN/重新聚类后查看是否有其他的细胞_FeaturePlot.png", height = 6000, width = 6000, res = 300)
plotFeature(scRNA_data = scRNA, choose = "Feature", col_num = 6, marker.list = marker.list)
dev.off()

png("/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/FastMNN/重新聚类后查看是否有其他的细胞_SCpubr.png", height = 6000, width = 6000, res = 300)
plotFeature(scRNA_data = scRNA, choose = "SCpubr", col_num = 6, marker.list = marker.list)
dev.off()




#################################################################
# 绘图
#################################################################
Features <-c("COL1A1","ACTA2","PDGFRB","MCAM","PDGFRA","DCN") # nolint
layout <-"ABBBBBB"
png("FeaturePlot.png",height=1500,width=10500,res=300)
(DimPlot(scRNA_sub,label = T,repel = T)+
    theme(legend.position = "right"))+
    FeaturePlot(scRNA_sub,features=Features,ncol=6)+
plot_layout(design=layout)
dev.off()

## 使用FastMNN中的大群，提取Endothelials和Fibroblasts亚群，使用Seurat CCA进行聚类，去除doublet cells
scRNA_sub <- scRNA_seurat[,!scRNA_seurat$integrated_snn_res.1 %in% c(16,19,11,10)]


## 绘图
plotFeature <- function(scRNA_seurat = scRNA_seurat) {
    features <- c("PDGFRA", "COL1A1", "ACTA2", "PDGFRB", "MCAM", "PECAM1", "CD34", "VWF","TOP2A","MKI67","STMN1","TUBA1B")
    p <- DimPlot(scRNA_seurat, label = T, repel = T) +
        theme(legend.position = "bottom") +
        FeaturePlot(scRNA_seurat, features = features, ncol = 6) +
        plot_layout(
            design = "AABBBBBB
                      AABBBBBB"
        )
    return(p)
}
png('harmony_O.3DimPlot结果.png', height = 2000, width = 8000, res= 300 )
plotFeature(scRNA_seurat = scRNA.harmony.res)
dev.off()

png('Seurat_CCA_O.3DimPlot结果.png', height = 2000, width = 8000, res = 300)
plotFeature(scRNA_seurat = scRNA.seurat)
dev.off()
