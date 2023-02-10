# load
library(Seurat)
library(patchwork)
library(tidyverse)
library(harmony)
# Clusters
RunCluster2 <- function(scRNA_data = scRNA_seurat,
                        choose = "harmony") {
    scRNA_data <- Seurat::CreateSeuratObject(scRNA_data@assays$RNA@counts)
    scRNA_data$percent.mt <- Seurat::PercentageFeatureSet(scRNA_data, pattern = "^MT")
    scRNA_data$sampleInfo <- stringr::str_split_fixed(rownames(scRNA_data@meta.data), "_[A|T|G|C]", n = 2)[, 1]
    sprintf("数据集包含%s个细胞,%s个基因", dim(scRNA_data)[2], dim(scRNA_data)[1])
    sprintf("数据集包含了%s个样本", length(scRNA_data$sampleInfo))
    if (choose == "harmony") {
        scRNA.harmony <- scRNA_data %>%
            NormalizeData() %>%
            FindVariableFeatures(nfeatures = 2000) %>%
            ScaleData(vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA")) %>%
            RunPCA(verbose = F) %>%
            RunHarmony(group.by.vars = "sampleInfo") %>%
            RunUMAP(reduction = "harmony", dims = 1:30) %>%
            FindNeighbors(reduction = "harmony", dims = 1:30) %>%
            FindClusters(resolution = seq(0.05, 1, 0.05))
        scRNA.harmony.res <<- scRNA.harmony
    } else if (choose == "SCT") {
        scRNAlist <- SplitObject(scRNA_data)
        scRNAlist <- lapply(X = scRNAlist, FUN = SCTransform)
        features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
        scRNAlist <- PrepSCTIntegration(object.list = scRNAlist, anchor.features = features)
        scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, normalization.method = "SCT", anchor.features = features)
        scRNA.sct <- IntegrateData(anchorset = scRNA.anchors, normalization.method = "SCT")
        scRNA.sct <- RunPCA(scRNA.sct, verbose = FALSE)
        scRNA.sct <- RunUMAP(scRNA.sct, reduction = "pca", dims = 1:30)
        scRNA.sct <- FindNeighbors(scRNA.sct, dims = 1:30, verbose = FALSE)
        scRNA.sct <- FindClusters(scRNA.sct, resolution = seq(0.05, 1, 0.05))
        scRNA.SCT.anchor <<- scRNA.anchors
        scRNA.sct.res <<- scRNA.sct
    }
}
