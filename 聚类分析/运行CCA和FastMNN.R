# library
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(SeuratWrappers)
library(SCpubr)
library(future)
# make function
RunCluster <- function(scRNA_data = scRNA_seurat,
                       choose = "CCA",
                       useFuture = FALSE) {
    if (useFuture) {
        plan(multisession, workers = 40)
        options(future.globals.maxSize = 300 * 1024^3)
    }

    scRNA_data <- Seurat::CreateSeuratObject(scRNA_data@assays$RNA@counts)
    scRNA_data$percent.mt <- Seurat::PercentageFeatureSet(scRNA_data, pattern = "^MT")
    scRNA_data$sampleInfo <- stringr::str_split_fixed(rownames(scRNA_data@meta.data), "_[A|T|G|C]", n = 2)[, 1]
    sprintf("数据集包含%s个细胞,%s个基因", dim(scRNA_data)[2], dim(scRNA_data)[1])
    sprintf("数据集包含了%s个样本", length(scRNA_data$sampleInfo))
    scRNAlist <- Seurat::SplitObject(scRNA_data, split.by = "sampleInfo")
    sprintf("样本的细胞量的范围是%s---->%s", range(table(scRNA_data$sampleInfo))[1], range(table(scRNA_data$sampleInfo))[2])
    scRNAlist <- suppressMessages(scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x)))
    scRNAlist <- suppressMessages(lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x)))
    if (is.null(choose)) {
        message("还没有选择聚类方式,请选择聚类方式")
    } else if (choose == "CCA") {
        scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)
        scRNA_CCA <- IntegrateData(anchorset = scRNA.anchors)
        scRNA_CCA <- ScaleData(scRNA_CCA, vars.to.regress = c("percent.mt", "nFeatureRNA", "nCount_RNA")) %>% RunPCA(verbose = FALSE)
        scRNA_CCA <- RunUMAP(scRNA_CCA, dims = 1:30)
        scRNA_CCA <- FindNeighbors(scRNA_CCA, dims = 1:30) %>% FindClusters(resolution = seq(0.05, 1, 0.05))
        scRNA.anchors.res <<- scRNA.anchors
        scRNA.CCA.res <<- scRNA_CCA
    } else if (choose == "FastMNN") {
        scRNA.FastMNN <- RunFastMNN(object.list = scRNAlist)
        scRNA.FastMNN <- RunUMAP(scRNA.FastMNN, reduction = "mnn", dims = 1:30)
        scRNA.FastMNN <- RunTSNE(scRNA.FastMNN, reduction = "mnn", dims = 1:30)
        scRNA.FastMNN <- FindNeighbors(scRNA.FastMNN, reduction = "mnn", dims = 1:30)
        scRNA.FastMNN <- FindClusters(scRNA.FastMNN, resolution = seq(0.05, 1, 0.05))
        scRNA.FastMNN.res <<- scRNA.FastMNN
    }
}

