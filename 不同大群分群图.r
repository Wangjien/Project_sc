###########################################
####### 大群
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse)

# 添加字体
setwd("/root/wangje/Project/刘老师/大群/Data/")
load("/root/wangje/Project/刘老师/大群/Data/大群CCA.Rdata")

scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(1, 11, 16, 4)), "new_celltype"] <- "NK & T cells"
scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(12, 14)), "new_celltype"] <- "B & Plasma cells"
# scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(12)), "new_celltype"] <- "Plasmas" # nolint
scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(8, 9, 17, 15, 0, 20)), "new_celltype"] <- "Myeloid cells"
scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(6, 21, 10)), "new_celltype"] <- "Fibroblasts"
scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(5, 3, 2, 19, 18)), "new_celltype"] <- "Cancer cells"
scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(18)), "new_celltype"] <- "Non-malignant epithelial cells"
scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(7, 13)), "new_celltype"] <- "Endothelial cells"

## 图片style
mytheme <- theme(
    plot.title = element_text(size = 20, colour = "black", hjust = 0.5, family = "Arial"),
    axis.text = element_text(size = 15, color = "black", family = "Arial"),
    axis.title = element_text(size = 17, color = "black", family = "Arial"),
    legend.text = element_text(size = 15, color = "black", family = "Arial"),
    legend.title = element_text(size = 17, color = "black", family = "Arial")
)

ids4 <- c(
    "CJME", "CJME_0707", "CMDI", "CMDI_0624", "HDYU", "HDYU_0720", "HXZH", "HXZH_0220", "LCLI", "LCLI_0623", "WYJU", "WYJU_0122",
    "ZXME", "ZXME_0223", "ZJLI_0116", "ZJLI_0312", "CZYI", "CZYI_0702", "FHXHBS1", "FHXHBS2", "HEJI", "HEJX", "LAWE", "LAWE_0309",
    "ZEYI", "ZEYI_0204", "WZLA", "FYYI", "ZFXI", "LIPE"
)
ids5 <- c(
    "P1-Pre", "P1-Post", "P2-Pre", "P2-Post", "P3-Pre", "P3-Post", "P4-Pre", "P4-Post", "P5-Pre", "P5-Post", "P6-Pre", "P6-Post",
    "P7-Pre", "P7-Post", "P8-Pre", "P8-Post", "P10-Pre", "P10-Post", "P11-Pre", "P11-Post", "P12-Pre", "P12-Post", "P13-Pre", "P13-Post",
    "P14-Pre", "P14-Post", "P9-Pre", "P15-Pre", "P16-Pre", "P17-Pre"
)

for (i in 1:length(ids4)) {
    print(ids4[i])
    print(ids5[i])
    scRNA_seurat@meta.data[which(scRNA_seurat$Patient == ids4[i]), "new_Rename"] <- ids5[i]
}
levels <- ids5
scRNA_seurat$new_Rename <- factor(scRNA_seurat$new_Rename, levels = levels)

add_umap_axes <- function(p) {
    umap_1_min <- min(p$data$UMAP_1)
    umap_2_min <- min(p$data$UMAP_2)

    p +
        geom_segment(aes(
            x = umap_1_min, y = umap_2_min,
            xend = umap_1_min + 4, yend = umap_2_min
        ), colour = "black", size = 0.5, arrow = arrow(length = unit(0.5, "cm"), type = "closed")) +

        geom_segment(aes(
            x = umap_1_min, y = umap_2_min,
            xend = umap_1_min, yend = umap_2_min + 4
        ), colour = "black", size = 0.5, arrow = arrow(length = unit(0.5, "cm"), type = "closed")) +

        annotate("text",
            x = umap_1_min + 2, y = umap_2_min - 1, label = "UMAP1",
            color = "black", size = 9, family = 'Arial'
        ) +

        annotate("text",
            x = umap_1_min - 1, y = umap_2_min + 2, label = "UMAP2",
            color = "black", size = 9,family='Arial', angle = 90
        )
}

p1 <- DimPlot(scRNA_seurat, group.by = "Tissue", label = F, raster = F, repel = T, label.size = 11, pt.size = 1e-28) +
    mytheme+ # 调整legend间的距离
    labs(title = "") +
    guides(colour = guide_legend(override.aes = list(size = 10)))

p1 <- add_umap_axes(p1)
png("./01.png", height = 3000, width = 3000, res = 300)
p1
dev.off()

p2 <- DimPlot(scRNA_seurat, group.by = "new_celltype", label = F, raster = F, repel = T, label.size = 11, pt.size = 1e-28) +
    mytheme + 
   labs(title = "") +
    guides(colour = guide_legend(override.aes = list(size = 10)))

p2 <- add_umap_axes(p2)
png("./02.png", height = 3000, width = 3000, res = 300)
p2
dev.off()


#####################################################
############## 样本信息
scRNA_seurat$sample <- stringr::str_split_fixed(scRNA_seurat$new_Rename, '-', n= 2)[,1]

p3 <- DimPlot(scRNA_seurat, group.by = "sample", label = F, raster = F, repel = T, label.size = 11, pt.size = 1e-28) +
    mytheme+
    labs(title = "") +
    guides(colour = guide_legend(override.aes = list(size = 10)))

p3 <- add_umap_axes(p3)
png("./03.png", height = 3000, width = 3000, res = 300)
p3
dev.off()

##################################################
######### TCR信息

