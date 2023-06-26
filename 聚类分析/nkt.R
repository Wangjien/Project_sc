library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(DoubletFinder)
# load data 
big.fastmnn = '/root/wangje/wangjien_new/liu-big/scRNA_FastMNN_new.RData'
load(big.fastmnn)
# 去除双细胞
scRNA_sub = scRNA[,scRNA$celltype == 'NK & T cells']
for(i in unique(scRNA_sub$sample)){
    tmp = scRNA_sub[,scRNA_sub$sample == i]
    tmp = CreateSeuratObject(test@assays$RNA@counts)
    tmp$percent.mt = PercentageFeatureSet(tmp, pattern='^MT-')
    tmp <- NormalizeData(tmp)
    tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
    tmp <- ScaleData(tmp)
    tmp <- RunPCA(tmp)
    tmp <- RunUMAP(tmp, dims = 1:50)
    ## pK Identification (NO ground-truth) ------------------------------------------------------------------------------------------
    sweep.res.list_kidney <- paramSweep_v3(tmp, PCs = 1:50, sct = FALSE)
    sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
    bcmvn_kidney <- find.pK(sweep.stats_kidney)

    ## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
    # sweep.res.list_kidney <- paramSweep_v3(tmp, PCs = 1:50, sct = FALSE)
    # gt.calls <- tmp@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"]  ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
    # sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
    # bcmvn_kidney <- find.pK(sweep.stats_kidney)

    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- tmp@meta.data$ClusteringResults
    nExp_poi <- round(0.075*nrow(tmp@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    tmp <- doubletFinder_v3(tmp, PCs = 1:50, pN = 0.25, pK = 0.003, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    tmp <- doubletFinder_v3(tmp, PCs = 1:50, pN = 0.25, pK = 0.003, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)


}

## 使用Seurat CCA 进行分析
scRNAlist = SplitObject(scRNA_sub, split.by = 'sample')
scRNAlist <- suppressMessages(scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x)))
scRNAlist <- suppressMessages(lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x)))
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)
scRNA_CCA <- IntegrateData(anchorset = scRNA.anchors)
scRNA_CCA <- ScaleData(scRNA_CCA, vars.to.regress = c("percent.mt", "nFeatureRNA", "nCount_RNA")) %>% RunPCA(verbose = FALSE)
scRNA_CCA <- RunUMAP(scRNA_CCA, dims = 1:50)
scRNA_CCA <- FindNeighbors(scRNA_CCA, dims = 1:50) %>% FindClusters(resolution = seq(0.05, 1, 0.05))
scRNA.anchors.res <<- scRNA.anchors
