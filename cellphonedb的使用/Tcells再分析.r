#--------------------------------------------------------------
# Tcells 重新分群
# 2023-06-18
#--------------------------------------------------------------
library(Seurat)
library(harmony)
library(dplyr)
library(patchwork)

# ----------------------------------------------------------------
# Harmony
#-----------------------------------------------------------------
big.fastmnn = '/root/wangje/wangjien_new/liu-big/scRNA_FastMNN_new.RData'
load(big.fastmnn)











# ----------------------------------------------------------------
# rPCA
#-----------------------------------------------------------------
scRNA = scRNA[,scRNA$celltype == 'NK & T cells'] # 101962 cells
scRNA = CreateSeuratObject(scRNA@assays$RNA@counts)
scRNA$sample = stringr::str_split_fixed(rownames(scRNA@meta.data), '_[A|T|G|C].*', n = 2)[,1]
ifnb.list <- SplitObject(scRNA, split.by = "sample")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
immune.combined <- IntegrateData(anchorset = immune.anchors)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 50, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:50)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:50)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# // 保存图片
p1 = DimPlot(scRN)

# // 选择不同的min.dist
plist = list()
flist = list()
for(dist in c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001)){
    print(as.character(dist))
    test.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:50,min.dist=dist)
    test.combined <- FindNeighbors(test.combined, reduction = "pca", dims = 1:50)
    test.combined <- FindClusters(test.combined, resolution = 0.2)
    # plot
    tmp = DimPlot(test.combined,reduction = 'umap') + labs(title = paste0('dist_',dist))
    plist[[as.character(dist)]] = tmp
    flist[[as.character(dist)]] = test.combined
} 
