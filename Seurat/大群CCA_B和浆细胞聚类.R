library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(SeuratWrappers)
library(harmony)


#######################################################################
################# 读入文件，使用CCA聚类
rm(list = ls())
gc()

load('/root/wangje/Project/刘老师/大群/大群CCA.Rdata')
scRNA <- scRNA_seurat[, scRNA_seurat$celltype %in% 'B & Plasma cells']
# 查看样本细胞数目，从小到大排序
sort(table(scRNA$Patient),decreasing = F)
# 样本的细胞数目小于30个，在进行锚点检测的时候会发生错误，可以将小于30个样本细胞结合
metaData <- scRNA@meta.data %>% 
                select(Rename, Tissue, Rseponse, Patient)
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts)
# meta.data中添加info
scRNA$percent.mt <- PercentageFeatureSet(scRNA, pattern = '^MT-')
scRNA@meta.data <- scRNA@meta.data %>%
                        mutate(
                            Rename = metaData[match(rownames(scRNA@meta.data),rownames(metaData)),'Rename'],
                            Tissue = metaData[match(rownames(scRNA@meta.data),rownames(metaData)),'Tissue'],
                            Rseponse = metaData[match(rownames(scRNA@meta.data),rownames(metaData)),'Rseponse'],
                            Patient = metaData[match(rownames(scRNA@meta.data),rownames(metaData)),'Patient']
                        )
# 合并细胞数目不足30个的样本
scRNA$new_Patient <- scRNA$Patient
scRNA$new_Patient = ifelse(
    scRNA$Patient %in% c('HDYU', 'WYJU', 'CJME_0707'), 'merge_Patient', scRNA$Patient
)
scRNAlist <- SplitObject(scRNA, split.by = 'new_Patient') %>% 
            lapply(function(x){NormalizeData(x)}) %>% 
            lapply(function(x)FindVariableFeatures(x))
scRNA <- scRNAlist %>% FindIntegrationAnchors()
scRNA <- IntegrateData(anchorset = scRNA, k.weight = 48)
scRNA <- scRNA %>% ScaleData(vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percengt.mt'))
scRNA <- scRNA %>% RunPCA(verbose =FALSE)
scRNA <- scRNA %>% RunUMAP(dims = 1:30) %>% FindNeighbors(dims = 1:30) %>% FindClusters(resolution = seq(0.05, 0.8, 0.05))
# 保存文件
saveRDS(scRNA, file = '/root/wangje/Project/刘老师/Bcells/Data/Bcells_CCA.RDS')
######################################################################
######### 使用SeuratWrappers中的的RunFastMNN
scRNAlist <- SplitObject(scRNA, split.by = 'new_Patient') %>% 
            lapply(function(x){NormalizeData(x)}) %>% 
            lapply(function(x)FindVariableFeatures(x))

scRNA <- scRNAlist %>% 
            RunFastMNN() %>% 
            RunUMAP(reduction = 'mnn', dims = 1:30) %>% 
            FindNeighbors(reduction = 'mnn', dims = 1:30) %>%
            FindClusters(reduction = 'mnn', dims = 1:30)

#######################################################################
########## 使用Harmony进行聚类分析
scRNA <- scRNA %>% 
            NormalizeData() %>% 
            FindVariableFeatures(nfeatures = 2000) %>%
            ScaleData(vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA")) %>%
            RunPCA(verbose = F) %>%
            RunHarmony(group.by.vars = "new_Patient") %>%
            RunUMAP(reduction = "harmony", dims = 1:30) %>%
            FindNeighbors(reduction = "harmony", dims = 1:30) %>%
            FindClusters(resolution = seq(0.05, 1, 0.05))

#######################################################################
######### 绘图
DefaultAssay(scRNA) <- 'RNA'
# 查看细胞中是否还有其他的细胞类型
NKT <- c("CD3D", "CD3G", "CD2")
Fibroblasts <- c("COL1A1", "DCN", "LUM")
Myeloids <- c("LYZ", "CD68", "TYROBP")
Epithelial <- c("CD24", "KRT19", "EPCAM")
Bcells <- c("CD79A", "CD19", "MS4A1")
Endothelial <- c("CLDN5", "FLT1", "RAMP2")
Plasma <- c("IGHG1", "JCHAIN", "MZB1")
Hepatocytes <- c("ALB", "APOB", "HP")
Mast <- c('CPA3','TPSAB1','TPSB2')
DC <- c('LILRA4','CXCR3','IRF7')
levels <- c("NK&T cells", "B cells", "Plasmas", "Myeloids", "Epithelials", "Endothelials", "Fibroblasts", "Hepatocytes")
marker.list <- list(
    "NK&T cell" = NKT, "B cells" = Bcells, Plasmas = Plasma, Myeloids = Myeloids,
    Epithelials = Epithelial, Endothelials = Endothelial, Fibroblasts = Fibroblasts, Hepatocytes = Hepatocytes
)
# 绘制并保存图片
## CCA
png("/root/wangje/Project/刘老师/Bcells/Fig/CCA_查看是否有其他类型细胞_FueaturePlot.png",height =5000,width = 6000,res=300)
plotFeature(scRNA_data=scRNA,
            choose="Feature",
            col_num=6,
            marker.list=marker.list)
dev.off()

png("/root/wangje/Project/刘老师/Bcells/Fig/CCA_查看是否有其他类型细胞_Scpubr.png",height =5000,width = 6000,res=300)
plotFeature(scRNA_data=scRNA,
            choose="SCpubr",
            col_num=6,
            marker.list=marker.list)
dev.off() 

## FastMNN
png("/root/wangje/Project/刘老师/Bcells/Fig/FastMNN_查看是否有其他类型细胞_FueaturePlot.png",height =5000,width = 6000,res=300)
plotFeature(scRNA_data=scRNA,
            choose="Feature",
            col_num=6,
            marker.list=marker.list)
dev.off()

png("/root/wangje/Project/刘老师/Bcells/Fig/FastMNN.png",height =5000,width = 6000,res=300)
plotFeature(scRNA_data=scRNA,
            choose="SCpubr",
            col_num=6,
            marker.list=marker.list)
dev.off() 

## Harmony
png("/root/wangje/Project/刘老师/Bcells/Fig/Harmony_查看是否有其他类型细胞_FueaturePlot.png",height =5000,width = 6000,res=300)
plotFeature(scRNA_data=scRNA,
            choose="Feature",
            col_num=6,
            marker.list=marker.list)
dev.off()

png("/root/wangje/Project/刘老师/Bcells/Fig/harmony_SCpubr.png",height =5000,width = 6000,res=300)
plotFeature(scRNA_data=scRNA,
            choose="SCpubr",
            col_num=6,
            marker.list=marker.list)
dev.off() 

# 绘制DimPlot和DimPlot以及其他信息
## DimPlot
Features <- c('CD79A', 'MS4A1', 'MZB1', 'CD38', 'IGHD', 'CD27', 'NEIL1', 'MKI67', 'IGHA1', 'IGHG1',
            'IGHM', 'FCER2', '')
p_Feature <- FeaturePlot(scRNA, features = Features, ncol = 5) + 
    theme(legend.position = 'bottom')

p_dimplot <- DimPlot(scRNA, 
                    group.by = 'integrated_snn_res.0.5',
                    label = T, 
                    repel =T) + 
                    theme(legend.position = 'bottom')

## 小提琴图
DefaultAssay(scRNA) <- 'RNA'
Idents(scRNA) <- 'integrated_snn_res.0.1'
p_vlinplot <- VlnPlot(scRNA,
                    features = Features, 
                    stack = T, 
                    group.by ='integrated_snn_res.0.1',
                    fill = 'ident') +
                    ggsci::scale_fill_igv()
## 绘制热图
library(Seurat)
library(patchwork)
library(ggplot2)
library(BiocParallel)
library(dplyr)
library(tidyverse)

Idents(scRNA) <- scRNA$integrated_snn_res.0.1
scRNA <- scRNA %>% ScaleData()
 scRNA.Findmarkers <- FindAllMarkers(scRNA, 
                                    only.pos = TRUE, 
                                    BPPARAM = MulticoreParam(20))
Top10.genes <- scRNA.Findmarkers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
p_Heatmap <- DoHeatmap(scRNA, 
                    features = Top10.genes$gene, 
                    size = 4) + NoLegend() + 
                    theme(text = element_text(size = 6))

