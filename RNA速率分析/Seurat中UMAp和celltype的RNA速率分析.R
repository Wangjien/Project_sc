library(Seurat)
library(dplyr)
library(ggplot2)
library(SCP)
library(qs)
library(destiny)
library(velocyto.R)
library(SeuratWrappers)
setwd('/root/wangje/Project/刘老师/NK_T/Data/CCA新的分群')
scRNA<- qread('CCA新的分群结果.qs')
# 对CD4细胞类型进行拟时间分析
scRNA_cd4 <- scRNA[,scRNA$celltype %in% c('Treg','SOX4+ CD4','TFH cells','Naive T cells','TH1-like')]




# 读入loom文件
# loom文件的地址在/root/wangje/Project/刘老师/Loom文件
loom_file <- read.loom.matrices('./merged.loom')
# 读入Seurat对象
scRNA <- qread('/root/wangje/Project/刘老师/NK_T/Data/CCA新的分群/CCA新的分群结果.qs')
test <- as.Seurat(loom_file) # 将loom文件转换成Seurat格式
rownames(test@meta.data) <- gsub(':','_',rownames(test@meta.data))



de_filter <- filter(scRNA@tools$DEtest_celltype$AllMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 1)
colors_use = colorRampPalette(c("#3288BD", "white", "red"))(50)
colors_use = colorRampPalette(c("#20528e", "white", "#790f1e"))(50)
ht2 <- GroupHeatmap(
  srt = scRNA, 
  features = de_filter$gene, 
  group.by = "celltype",
#   split.by = "Phase", 
  cell_split_palette = "Dark2",
  nlabel = 20, show_row_names = FALSE,
  heatmap_palcolor = colors_use
)
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
ht2$plot

ht3 <- GroupHeatmap(
  srt = scRNA, 
  features = de_filter$gene, 
  feature_split = de_filter$group1, 
  group.by = "celltype",
  nlabel = 50, show_row_names = FALSE,
  species = "Homo_sapiens", db = "GO_BP", 
  anno_terms = TRUE, 
  anno_keys = TRUE, 
  heatmap_palcolor = colors_use,
  anno_features = TRUE,
  GO_simplify = TRUE
)


scRNA$celltype <- ifelse(scRNA$celltype %in% c('Keratinocytes','Hepatocytes'),'Non-malignant epithelial cells',scRNA$celltype)
scRNA$celltype <- factor(scRNA$celltype, levels = c('Myeloid cells' ,'Epithelial cells' ,'Endothelial cells' ,'Fibroblasts' ,'NK & T cells' ,'B & Plasma cells','Non-malignant epithelial cells'))

ht <- FeatureHeatmap(
  srt = scRNA, 
  group.by = "celltype", 
  features = DEGs$gene, 
  feature_split = DEGs$group1,
  species = "Homo_sapiens", 
  db = c("GO_BP"), anno_terms = TRUE,
#   feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  heatmap_palcolor = colors_use,
  height = 8, width = 6
)


 ht <- FeatureHeatmap(
      srt = scRNA,
      group.by = "celltype",
      features = DEGs$gene,
      feature_split = DEGs$group1,
      species = "Homo_sapiens",
      db = c("GO_BP"), anno_terms = TRUE,
    #   feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
      heatmap_palcolor = colors_use,
      height = 8, width = 6,topTerm = 5, 
      terms_fontsize = 13 # 调整GO Term字体大小
    )



mycolors = c('#73c8c8','#d3c2dd','#e3852c','#babb3f','#a87cba','#5289ad','#cb5b52','#4da846')
# CellDimPlot(scRNA, group.by = 'celltype', xlab = 'UMAP1', ylab = 'UMAP2') + Seurat::NoLegend() + 
#     theme(axis.title = element_text(size = 20, colour = 'black'), axis.text = element_text(size = 18, colour = 'black'),
#     ,axis.line = element_line(linewidth = 1,colour = 'black'))-> p1
# ggsave(filename = './DimPlot01.pdf', height = 5, width= 5, plot = p1, bg = 'white')
# ggsave(filename = './DimPlot01.png', height = 5, width= 5, plot = p1, bg = 'white')


CellDimPlot(scRNA, group.by = 'celltype', xlab = 'UMAP1', ylab = 'UMAP2') + Seurat::NoLegend() + 
    theme(axis.title = element_text(size = 20, colour = 'black'), axis.text = element_text(size = 18, colour = 'black'),
    ,panel.border = element_rect(fill=NA,color="black", size=1,linetype="solid"))-> p1
ggsave(filename = './DimPlot01.pdf', height = 5, width= 5, plot = p1, bg = 'white')
ggsave(filename = './DimPlot01.png', height = 5, width= 5, plot = p1, bg = 'white')

## Group
mycolors <- c('#a9cce1','#2077b3','#b0e08b','#2ea327')
CellDimPlot(scRNA, group.by = 'Group', xlab = 'UMAP1', ylab = 'UMAP2', show_stat = FALSE,palcolor=mycolors ) + Seurat::NoLegend() + 
    theme(axis.title = element_text(size = 20, colour = 'black'), axis.text = element_text(size = 18, colour = 'black'),
    panel.border = element_blank()) + Seurat::NoAxes() -> p2

ggsave(filename = './DimPlot02.pdf', height = 5, width= 5, plot = p2, bg = 'white')
ggsave(filename = './DimPlot02.png', height = 5, width= 5, plot = p2, bg = 'white')

## Response
CellDimPlot(scRNA, group.by = '', xlab = 'UMAP1', ylab = 'UMAP2', show_stat = FALSE) + Seurat::NoLegend() + 
    theme(axis.title = element_text(size = 20, colour = 'black'), axis.text = element_text(size = 18, colour = 'black'),
    panel.border = element_blank()) + Seurat::NoAxes() -> p2

ggsave(filename = './DimPlot02.pdf', height = 5, width= 5, plot = p2, bg = 'white')
ggsave(filename = './DimPlot02.png', height = 5, width= 5, plot = p2, bg = 'white')

## sample
source('/root/wangje/Project/刘老师/script/AddInformation.R')
scRNA <- AddInfo(scRNA)
scRNA$sample <- stringr::str_split_fixed(scRNA$new_Rename,'-', n = 2)[,1]
scRNA$sample <- factor(scRNA$sample, levels=paste0('P',1:17))

CellDimPlot(scRNA, group.by = 'sample', xlab = 'UMAP1', ylab = 'UMAP2', show_stat = FALSE) + Seurat::NoLegend() + 
    theme(axis.title = element_text(size = 20, colour = 'black'), axis.text = element_text(size = 18, colour = 'black'),
    panel.border = element_blank()) + Seurat::NoAxes() -> p3
ggsave(filename = './DimPlot_sample.pdf', height = 5, width= 5, plot = p3, bg = 'white')
ggsave(filename = './DimPlot_sample.png', height = 5, width= 5, plot = p3, bg = 'white')

# TCR
scRNA$TCR <- factor(scRNA$TCR, levels = c('TCR','No TCR'))