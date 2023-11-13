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
  species = "Homo_sapiens", 
  db = "KEGG", # 富集类型 
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

colors_use = colorRampPalette(c("#a54894", "black", "#e5d952"))(50)
ht <- FeatureHeatmap(
    srt = scRNA,
    group.by = "celltype",
    split.by = 'Group',
    features = new_DEGs$gene,
    feature_split = new_DEGs$group1,
    species = "Homo_sapiens",
    db = c("KEGG"), # 富集类型 
    anno_terms = FALSE,
    keys_width = unit(10, "in"),
    terms_width = unit(10, "in"),
  #   feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
    heatmap_palcolor = colors_use,
    show_row_names = FALSE,
    nlabel = 0,
    height = 16, width = 6,topTerm = 5, 
    # ht_params = list(width = unit(4, "cm")),
    terms_fontsize = 0 # 调整GO Term字体大小
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
mycolors <- c('#2c91a4','#a04745','#346799','#8a9d4e')
CellDimPlot(scRNA, group.by = 'Group', xlab = 'UMAP1', ylab = 'UMAP2', show_stat = FALSE,palcolor=mycolors ) + Seurat::NoLegend() + 
    theme(axis.title = element_text(size = 20, colour = 'black'), axis.text = element_text(size = 18, colour = 'black'),
    panel.border = element_blank()) + Seurat::NoAxes() -> p2

ggsave(filename = './DimPlot02.pdf', height = 5, width= 5, plot = p2, bg = 'white')
ggsave(filename = './DimPlot02.png', height = 5, width= 5, plot = p2, bg = 'white')

# 加上legend
mycolors <- c('#2c91a4','#a04745','#346799','#8a9d4e')
CellDimPlot(scRNA, group.by = 'Group', xlab = 'UMAP1', ylab = 'UMAP2', show_stat = FALSE,palcolor=mycolors )+ 
    theme(axis.title = element_text(size = 20, colour = 'black'), axis.text = element_text(size = 18, colour = 'black'),
    panel.border = element_blank()) + Seurat::NoAxes() -> p2
p2  <- p2 + guides(color=guide_legend(override.aes=list(size=5)))    

ggsave(filename = './DimPlot02_2.pdf', height = 5, width= 5, plot = p2, bg = 'white')
ggsave(filename = './DimPlot02_2.png', height = 5, width= 5, plot = p2, bg = 'white')


## Response
CellDimPlot(scRNA, group.by = 'Response', xlab = 'UMAP1', ylab = 'UMAP2', show_stat = FALSE) + Seurat::NoLegend() + 
    theme(axis.title = element_text(size = 20, colour = 'black'), axis.text = element_text(size = 18, colour = 'black'),
    panel.border = element_blank()) + Seurat::NoAxes() -> p2

ggsave(filename = './DimPlot_Responde.pdf', height = 5, width= 5, plot = p2, bg = 'white')
ggsave(filename = './DimPlot_Response.png', height = 5, width= 5, plot = p2, bg = 'white')

## sample
source('/root/wangje/Project/刘老师/script/AddInformation.R')
scRNA <- AddInfo(scRNA)
scRNA$sample <- stringr::str_split_fixed(scRNA$new_Rename,'-', n = 2)[,1]
scRNA$sample <- factor(scRNA$sample, levels=paste0('P',1:17))

CellDimPlot(scRNA, group.by = 'sample', xlab = 'UMAP1', ylab = 'UMAP2', show_stat = FALSE, legend.direction= 'horizontal', legend.position = 'bottom') + 
    theme(axis.title = element_text(size = 20, colour = 'black'), axis.text = element_text(size = 18, colour = 'black'),
    panel.border = element_blank()) + Seurat::NoAxes() -> p3
p3 <- p3 + guides(color=guide_legend(override.aes=list(size=6)))    
ggsave(filename = './DimPlot_sample.pdf', height = 5, width= 8, plot = p3, bg = 'white')
ggsave(filename = './DimPlot_sample.png', height = 5, width= 8, plot = p3, bg = 'white')

# 没有legend
CellDimPlot(scRNA, group.by = 'sample', xlab = 'UMAP1', ylab = 'UMAP2', show_stat = FALSE, legend.direction= 'horizontal', legend.position = 'bottom') + 
    theme(axis.title = element_text(size = 20, colour = 'black'), axis.text = element_text(size = 18, colour = 'black'),
    panel.border = element_blank()) + Seurat::NoAxes() -> p3
p3 <- p3 + Seurat::NoLegend()  
ggsave(filename = './DimPlot_sample.pdf', height = 5, width= 5, plot = p3, bg = 'white')
ggsave(filename = './DimPlot_sample.png', height = 5, width= 5, plot = p3, bg = 'white')


# TCR
scRNA$TCR <- factor(scRNA$TCR, levels = c('TCR','No TCR'))
mycolors <- c('#81c1d7','#b9b9b5')
CellDimPlot(scRNA, group.by = 'TCR', xlab = 'UMAP1', ylab = 'UMAP2', palcolor = mycolors,show_stat = FALSE, legend.direction= 'horizontal', legend.position = 'bottom') + 
    theme(axis.title = element_text(size = 20, colour = 'black'), axis.text = element_text(size = 18, colour = 'black'),
    panel.border = element_blank()) + Seurat::NoAxes() -> p3
# p3 <- p3 + guides(color=guide_legend(override.aes=list(size=6)))
p3 <- p3 + Seurat::NoLegend()

ggsave(filename = './DimPlot_TCR.pdf', height = 5, width= 5, plot = p3, bg = 'white')
ggsave(filename = './DimPlot_TCR.png', height = 5, width= 5, plot = p3, bg = 'white')

ids5 <- c(
    "P1-Pre", "P1-Post", "P2-Pre", "P2-Post", "P3-Pre", "P3-Post", "P4-Pre", "P4-Post", "P5-Pre", "P5-Post", "P6-Pre", "P6-Post",
    "P7-Pre", "P7-Post", "P8-Pre", "P8-Post", "P9-Pre","P10-Pre", "P10-Post", "P11-Pre", "P11-Post", "P12-Pre", "P12-Post", "P13-Pre", "P13-Post",
    "P14-Pre", "P14-Post", "P15-Pre", "P16-Pre", "P17-Pre"
)
scRNA$new_Rename <- factor(scRNA$new_Rename , levels= ids5)
scRNA$sample  <- stringr::str_split_fixed(scRNA$new_Rename, '-', n=2)[,1]
scRNA$sample  <- factor(scRNA$sample, levels = paste0('P',1:17))

mycolors  <- c('#da7100','#a976bc','#bd2918','#2c71a2')
CellDimPlot(scRNA, group.by = 'Tissue', xlab = 'UMAP1', ylab = 'UMAP2',palcolor = mycolors,show_stat = FALSE, legend.direction= 'horizontal', legend.position = 'bottom') + 
    theme(axis.title = element_text(size = 20, colour = 'black'), axis.text = element_text(size = 18, colour = 'black'),
    panel.border = element_blank()) + Seurat::NoAxes() -> p3
# p3 <- p3 + guides(color=guide_legend(override.aes=list(size=6)))   
p3 <- p3 + Seurat::NoLegend() 
ggsave(filename = './DimPlot_Tissue.pdf', height = 5, width= 5, plot = p3, bg = 'white')
ggsave(filename = './DimPlot_Tissue.png', height = 5, width= 5, plot = p3, bg = 'white')




p1 <- CellStatPlot(scRNA, stat.by = "celltype", group.by = "new_Rename", label = FALSE, plot_type = 'trend', 
  legend.direction = 'horizontal',legend.position = 'top')
p1 + guides(color=guide_legend(override.aes=list(size=6))) + 
    theme(legend.text = element_text(size = 10, color = 'black'),axis.title.x = element_blank())


