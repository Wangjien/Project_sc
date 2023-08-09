
suppressMessages({
    library(Seurat)
    library(qs)
    library(dplyr)
    library(ggplot2)
    library(SCP)
    library(stringr)
})

#@ 添加TCR信息
tcr <- read.csv("/root/wangje/Project/刘老师/tcr/new_contig.merge.csv",sep = ',', header = F)
scRNA$TCR <- ifelse(rownames(scRNA@meta.data) %in% tcr$V1,'TCR','No TCR')



p1 = SCP::CellDimPlot(scRNA,group.by = c('celltype'), ncol = 1, raster = F, label = T) + 
    labs(x = "UMAP_1", y = "UMAP_2")

p2 = SCP::CellDimPlot(scRNA,group.by = c('Treat_assess'), ncol = 1, raster = F) + 
    labs(x = "UMAP_1", y = "UMAP_2", subtitle = " ")

p3 = SCP::CellDimPlot(scRNA,group.by = c('Tissue'), ncol = 1, raster = F) + 
    labs(x = "UMAP_1", y = "UMAP_2", subtitle = " ")

p4 = SCP::CellDimPlot(scRNA,group.by = c('TCR'), ncol = 1,raster = F,palcolor = list('grey','red')) + 
    labs(x = "UMAP_1", y = "UMAP_2", subtitle = " ")

ggsave(filename = "./03_DimPlot.png",height = 10, width =15 ,limitsize = F, plot = patchwork::wrap_plots(list(p1,p2,p3,p4),ncol = 2), bg = 'white')

levels <- c("NK & T cells","Fibroblasts","Myeloid cells","Epithelial cells","B & Plasma cells","Endothelial cells","Hepatocytes","Keratinocytes")
scRNA$celltype <- factor(scRNA$celltype, levels = levels)


p5 = SCP::CellDimPlot(scRNA,group.by = c('new_Rename'), ncol = 1,raster = F) + 
    labs(x = "UMAP_1", y = "UMAP_2", subtitle = " ")

p6 = SCP::CellDimPlot(scRNA,group.by = c('Response'), ncol = 1,raster = F) + 
    labs(x = "UMAP_1", y = "UMAP_2", subtitle = " ")
ggsave(filename = "./04_DimPlot.png",height = 5, width =15 ,limitsize = F, plot = patchwork::wrap_plots(list(p5,p6),ncol = 2), bg = 'white')


#@ 查看不用样本的细胞数量
scRNA@meta.data %>% 
	group_by(sample) %>% 
	dplyr::summarize(cells = n(),
		median_nCount_RNA = median(nCount_RNA), 
		median_nFeature_RNA = median(nFeature_RNA),
		median_percent.mt = median(percent.mt),
		range_nCount_RNA = paste0(min(nCount_RNA), '->' ,max(nCount_RNA)),
		range_nFeature_RNA = paste0(min(nFeature_RNA), '->', max(nFeature_RNA)),
		range_percent.mt = paste0(min(percent.mt),'->', max(percent.mt))) %>% arrange(cells) -> data
ggpubr::as_ggplot(gridExtra::tableGrob(data)) -> p3

ggsave(filename = "./01_cellNum.png",height =15, width =30 ,limitsize = F, plot = p3, bg = 'white')

#@ 绘制小提琴图
Idents(scRNA) = scRNA$Rename
p2 <- VlnPlot(scRNA, features = c('nCount_RNA','nFeature_RNA','percent.mt'), ncol = 3, raster =F)
ggsave(filename = "./02_VlnPlot.png",height = 8, width =30 ,limitsize = F, plot = p2, bg = 'white')

p2 <- VlnPlot(scRNA, features = c('nCount_RNA','nFeature_RNA','percent.mt'), ncol = 3, raster =F, pt.size = 0)
ggsave(filename = "./02_VlnPlot_2.png",height = 8, width =30 ,limitsize = F, plot = p2, bg = 'white')

#@ 查看不同分组的细胞变化
p7 = SCP::CellDimPlot(scRNA,group.by = c('celltype'), split.by = 'Treat_assess', raster = F, ncol = 4) + 
    labs(x = "UMAP_1", y = "UMAP_2", subtitle = " ")

p8 = SCP::CellDimPlot(scRNA,group.by = c('celltype'), split.by = 'celltype', raster = F, ncol = 4) + 
    labs(x = "UMAP_1", y = "UMAP_2", subtitle = " ")



#@ 绘制热图
NKT <- c("CD3D", "CD3G", "CD2",'KLRC1','KLRC2')
Fibroblasts <- c("COL1A1", "DCN", "LUM")
Myeloids <- c("LYZ", "CD68", "TYROBP")
Epithelial <- c("CD24", "KRT19", "EPCAM")
Bcells <- c("CD79A", "CD19", "MS4A1","IGHG1", "JCHAIN", "MZB1")
Endothelial <- c("CLDN5", "FLT1", "RAMP2")
Hepatocytes <- c("ALB", "APOB", "HP")
Keratinocytes <- c("KRT5", "KRT14", "FABP5")
DC <- c('LILRA4','CXCR3','IRF7')
Mast <- c('CPA3','TPSAB1','TPSB2')
marker.list <- list("NK & T cell" = NKT, "B cell" = Bcells, Myeloids = Myeloids, 
  Fibroblasts = Fibroblasts, Epithelials = Epithelial,
  Endothelials = Endothelial, Hepatocytes = Hepatocytes, 
  Keratinocytes = Keratinocytes, DC =DC, Mast = Mast)
