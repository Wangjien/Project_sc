# library
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(cowplot)
library(stringr)
# load data
setwd("/root/wangje/Project/刘老师/大群/")
load("/root/wangje/Project/刘老师/大群/大群CCA.Rdata")
# add celltype
scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(1, 11, 16, 4)), "celltype"] <- "NK & T cells"
scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(12, 14)), "celltype"] <- "B & Plasma cells"
# scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(12)), "celltype"] <- "Plasmas" # nolint
scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(8, 9, 17, 15, 0, 20)), "celltype"] <- "Myeloid cells"
scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(6, 21, 10)), "celltype"] <- "Fibroblasts"
scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(5, 3, 2, 19, 18)), "celltype"] <- "Epithelial cells"
# scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(18)), "celltype"] <- "Non-malignant epithelial cells"
scRNA_seurat@meta.data[which(scRNA_seurat@meta.data$integrated_snn_res.0.45 %in% c(7, 13)), "celltype"] <- "Endothelial cells"

# 气泡图颜色
heatmap_color <- RColorBrewer::brewer.pal(name = "RdBu", n = 11)
pal <- rev(colorRampPalette(heatmap_color)(100))

# 循环
NKT <- c("CD3D", "CD3G", "CD2")
Fibroblasts <- c("COL1A1", "DCN", "LUM")
Myeloids <- c("LYZ", "CD68", "TYROBP")
Epithelial <- c("CD24", "KRT19", "EPCAM")
Bcells <- c("CD79A", "CD19", "MS4A1")
Endothelial <- c("CLDN5", "FLT1", "RAMP2")
Plasma <- c("IGHG1", "JCHAIN", "MZB1")
Hepatocytes <- c("ALB", "APOB", "HP")
Keratinocytes <- c("KRT5", "KRT14", "FABP5")
marker.list <- list("NK&T cell" = NKT, "B cell" = Bcells, "Plasmas" = Plasma, Myeloids = Myeloids, Fibroblasts = Fibroblasts, Epithelials = Epithelial, Endothelials = Endothelial, Hepatocytes = Hepatocytes, Keratinocytes = Keratinocytes)


# 循环绘制不同分辨率下的DimPlot和DotPlot
for (index in grep("^integrated", colnames(scRNA@meta.data))) {
    pdim <- scRNA %>% DimPlot(reduction = "umap", label = T, label.box = F, label.size = 6, repel = T, group.by = colnames(scRNA@meta.data)[index], label.color = "black", raster = F) +
        annotate(geom = "segment", y = Inf, yend = Inf, color = "black", x = -Inf, xend = Inf, size = 1) +
        annotate(geom = "segment", x = Inf, xend = Inf, color = "black", y = -Inf, yend = Inf, size = 1) +
        theme(legend.position = "bottom") + labs(title = colnames(scRNA@meta.data)[index])

    print(index)
    Idents(scRNA) <- scRNA@meta.data[, index]
    pDot <- DotPlot(scRNA, assay = "RNA", dot.scale = 4, features = marker.list) +
        scale_color_gradientn(colours = pal) +
        annotate(geom = "segment", y = Inf, yend = Inf, color = "black", x = -Inf, xend = Inf, size = 1) +
        annotate(geom = "segment", x = Inf, xend = Inf, color = "black", y = -Inf, yend = Inf, size = 1) +
        annotate(geom = "segment", x = Inf, xend = Inf, color = "black", y = -Inf, yend = Inf, size = 1) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 8), strip.text.x = element_text(size = 8, angle = 360)) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
        labs(x = "Gene Marker", y = "Cluster") + theme(legend.position = "bottom")
    pwd <- "./"
    outfile <- paste0(pwd, "Image_大群CCA_", colnames(scRNA@meta.data)[index], ".png")

    layout <- "
  AABBBBBB"
    png(filename = outfile, height = 2000, width = 5000, res = 300)
    print(pdim + pDot + plot_layout(design = layout) + plot_annotation(title = "", tag_levels = "A"))
    dev.off()
}

# 调整legend顺序
# scRNA_seurat$celltype <- factor(scRNA_seurat$celltype,levels = c("NK & T cells","B & Plasma cells","Myeloid cells",
# "Epithelial cells","Endothelial cells","Fibroblasts","Non-malignant epithelial cells"))
scRNA_seurat$celltype <- factor(scRNA_seurat$celltype, levels = c(
    "NK & T cells", "B & Plasma cells", "Myeloid cells",
    "Epithelial cells", "Endothelial cells", "Fibroblasts"
))

p <- DimPlot(scRNA_seurat, group.by = "celltype")
p1 <- DimPlot(scRNA_seurat, group.by = "celltype", label = F, raster = F, repel = T, label.size = 11, pt.size = 0.0001) +
    theme(
        legend.position = "bottom", axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.line.x = element_blank(),
        axis.line.y = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.text = element_text(size = 20), legend.key.height = unit(2, "line")
    ) + # 调整legend间的距离
    labs(title = "") +
    guides(colour = guide_legend(override.aes = list(size = 8))) +
    geom_segment(
        aes(
            x = min(p$data$UMAP_1), y = min(p$data$UMAP_2),
            xend = min(p$data$UMAP_1) + 4, yend = min(p$data$UMAP_2)
        ),
        colour = "black", size = 0.15, arrow = arrow(length = unit(0.3, "cm"), type = "closed")
    ) +
    geom_segment(
        aes(
            x = min(p$data$UMAP_1), y = min(p$data$UMAP_2),
            xend = min(p$data$UMAP_1), yend = min(p$data$UMAP_2) + 4
        ),
        colour = "black", size = 0.15, arrow = arrow(length = unit(0.3, "cm"), type = "closed")
    ) +
    annotate("text",
        x = min(p$data$UMAP_1) + 2, y = min(p$data$UMAP_2) - 1, label = "UMAP1",
        color = "black", size = 6, fontface = "bold"
    ) +
    annotate("text",
        x = min(p$data$UMAP_1) - 1, y = min(p$data$UMAP_2) + 2, label = "UMAP2",
        color = "black", size = 6, fontface = "bold", angle = 90
    )
png("./Image_Seurat_DimPlot_celltype_2.png", height = 2000, width = 2000, res = 300)
p1
dev.off()

png("./Image_Seurat_DimPlot_celltype_NoLegend.png", height = 2000, width = 2000, res = 300)
p1 + theme(legend.position = "none")
dev.off()



# ====== DimPlot by orgian ====================
p <- DimPlot(scRNA_seurat, group.by = "Tissue")
p3 <- DimPlot(scRNA_seurat, group.by = "Tissue", label = F, raster = F, repel = T, label.size = 11, pt.size = 1e-28) +
    theme(
        legend.position = "bottom", axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.line.x = element_blank(),
        axis.line.y = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.text = element_text(size = 20), legend.key.height = unit(2, "line")
    ) + # 调整legend间的距离
    labs(title = "") +
    guides(colour = guide_legend(override.aes = list(size = 8))) +
    geom_segment(
        aes(
            x = min(p$data$UMAP_1), y = min(p$data$UMAP_2),
            xend = min(p$data$UMAP_1) + 4, yend = min(p$data$UMAP_2)
        ),
        colour = "black", size = 0.15, arrow = arrow(length = unit(0.3, "cm"), type = "closed")
    ) +
    geom_segment(
        aes(
            x = min(p$data$UMAP_1), y = min(p$data$UMAP_2),
            xend = min(p$data$UMAP_1), yend = min(p$data$UMAP_2) + 4
        ),
        colour = "black", size = 0.15, arrow = arrow(length = unit(0.3, "cm"), type = "closed")
    ) +
    annotate("text",
        x = min(p$data$UMAP_1) + 2, y = min(p$data$UMAP_2) - 1, label = "UMAP1",
        color = "black", size = 6, fontface = "bold"
    ) +
    annotate("text",
        x = min(p$data$UMAP_1) - 1, y = min(p$data$UMAP_2) + 2, label = "UMAP2",
        color = "black", size = 6, fontface = "bold", angle = 90
    )
# 保存图片
png(filename = "Image_Seurat_DimPlot_Organ.png", height = 2000, width = 2000, res = 300)
p3
dev.off()
ggsave(filename = "Image_Seurat_DimPlot_Organ.pdf", height = 12, width = 12, dpi = 600, plot = p3)

# 保存图片没有legengd
png(filename = "Image_Seurat_DimPlot_Organ_Nolegend.png", height = 2000, width = 2000, res = 300)
p3 + theme(legend.position = "none")
dev.off()
ggsave(filename = "Image_Seurat_DimPlot_Organ_Nolegend.pdf", height = 12, width = 12, dpi = 600, plot = p3 + theme(legend.position = "none"))

# ========= DimPlot by Rename============
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
    scRNA_seurat@meta.data[which(scRNA_seurat$Patient == ids4[i]), "Rename"] <- ids5[i]
}
levels <- ids5
scRNA_seurat$new_Rename <- factor(scRNA_seurat$new_Rename, levels = levels)

## 修改名称
scRNA_seurat$new_Patient <- str_split_fixed(scRNA_seurat$Patient, "_", n = 2)[, 1]
scRNA_seurat$new_Patient <- gsub("FHXHBS1|FHXHBS2", "FHXHBS", scRNA_seurat$new_Patient)
scRNA_seurat$new_Patient <- gsub("HEJX", "HEJI", scRNA_seurat$new_Patient)

scRNA_seurat$Rename<- str_split_fixed(scRNA_seurat$Rename,"-",n=2)[,1]

# ---绘图----
p <- DimPlot(scRNA_seurat, group.by = "new_patient")
p4 <- DimPlot(scRNA_seurat, group.by = "new_patient", label = F, raster = F, repel = T, label.size = 11, pt.size = 1e-28) +
    theme(
        legend.position = "bottom", axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.line.x = element_blank(),
        axis.line.y = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.text = element_text(size = 20), legend.key.height = unit(2, "line")
    ) + # 调整legend间的距离
    labs(title = "") +
    guides(colour = guide_legend(override.aes = list(size = 8))) +
    geom_segment(
        aes(
            x = min(p$data$UMAP_1), y = min(p$data$UMAP_2),
            xend = min(p$data$UMAP_1) + 4, yend = min(p$data$UMAP_2)
        ),
        colour = "black", size = 0.15, arrow = arrow(length = unit(0.3, "cm"), type = "closed")
    ) +
    geom_segment(
        aes(
            x = min(p$data$UMAP_1), y = min(p$data$UMAP_2),
            xend = min(p$data$UMAP_1), yend = min(p$data$UMAP_2) + 4
        ),
        colour = "black", size = 0.15, arrow = arrow(length = unit(0.3, "cm"), type = "closed")
    ) +
    annotate("text",
        x = min(p$data$UMAP_1) + 2, y = min(p$data$UMAP_2) - 1, label = "UMAP1",
        color = "black", size = 6, fontface = "bold"
    ) +
    annotate("text",
        x = min(p$data$UMAP_1) - 1, y = min(p$data$UMAP_2) + 2, label = "UMAP2",
        color = "black", size = 6, fontface = "bold", angle = 90
    )
png("./Image_Seurat_DimPlot_Patients.png", height = 2000, width = 2000, res = 300)
p4
dev.off()

png("./Image_Seurat_DimPlot_Patients_NoLegend.png", height = 2000, width = 2000, res = 300)
p4 + theme(legend.position = "none")
dev.off()

## -------------------------------- Treatment ----------------------------------------------------
unique(scRNA_seurat$Treat_assess)
scRNA_seurat$Treat_assess <- factor(scRNA_seurat$Treat_assess, levels = c("R_Pre", "R_Post", "NR_Pre", "NR_Post"))
p <- DimPlot(scRNA_seurat, group.by = "Treat_assess")
p5 <- DimPlot(scRNA_seurat, group.by = "Treat_assess", label = F, raster = F, repel = T, label.size = 11, pt.size = 1e-28) +
    theme(
        legend.position = "bottom", axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.line.x = element_blank(),
        axis.line.y = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.text = element_text(size = 20), legend.key.height = unit(2, "line")
    ) + # 调整legend间的距离
    labs(title = "") +
    guides(colour = guide_legend(override.aes = list(size = 8))) +
    geom_segment(
        aes(
            x = min(p$data$UMAP_1), y = min(p$data$UMAP_2),
            xend = min(p$data$UMAP_1) + 4, yend = min(p$data$UMAP_2)
        ),
        colour = "black", size = 0.15, arrow = arrow(length = unit(0.3, "cm"), type = "closed")
    ) +
    geom_segment(
        aes(
            x = min(p$data$UMAP_1), y = min(p$data$UMAP_2),
            xend = min(p$data$UMAP_1), yend = min(p$data$UMAP_2) + 4
        ),
        colour = "black", size = 0.15, arrow = arrow(length = unit(0.3, "cm"), type = "closed")
    ) +
    annotate("text",
        x = min(p$data$UMAP_1) + 2, y = min(p$data$UMAP_2) - 1, label = "UMAP1",
        color = "black", size = 6, fontface = "bold"
    ) +
    annotate("text",
        x = min(p$data$UMAP_1) - 1, y = min(p$data$UMAP_2) + 2, label = "UMAP2",
        color = "black", size = 6, fontface = "bold", angle = 90
    )

png("./Image_Seurat_DimPlot_Treatment.png", height = 2000, width = 2000, res = 300)
p5
dev.off()

png("./Image_Seurat_DimPlot_Treatment_NoLegend.png", height = 2000, width = 2000, res = 300)
p5 + theme(legend.position = "none")
dev.off()

## --------------------------------------------- TCR ------------------------------------------------
TCR <- read.csv("/root/wangje/Project/刘老师/大群/new_contig.merge.csv", header = F, fill = T)
scRNA_seurat$TCR <- "No TCR"
TCR_rownames <- na.omit(match(TCR$V1, rownames(scRNA_seurat@meta.data)))
scRNA_seurat@meta.data[TCR_rownames, ]$TCR <- "TCR"

scRNA_seurat$TCR <- factor(scRNA_seurat$TCR, levels = c("TCR", "No TCR"))
p <- DimPlot(scRNA_seurat, group.by = "TCR")
p6 <- DimPlot(scRNA_seurat, group.by = "TCR", label = F, raster = F, repel = T, label.size = 11, pt.size = 1e-28, cols = c("red", "grey")) +
    theme(
        legend.position = "bottom", axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.line.x = element_blank(),
        axis.line.y = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.text = element_text(size = 20), legend.key.height = unit(2, "line")
    ) + # 调整legend间的距离
    labs(title = "") +
    guides(colour = guide_legend(override.aes = list(size = 10))) +
    geom_segment(
        aes(
            x = min(p$data$UMAP_1), y = min(p$data$UMAP_2),
            xend = min(p$data$UMAP_1) + 4, yend = min(p$data$UMAP_2)
        ),
        colour = "black", size = 0.15, arrow = arrow(length = unit(0.3, "cm"), type = "closed")
    ) +
    geom_segment(
        aes(
            x = min(p$data$UMAP_1), y = min(p$data$UMAP_2),
            xend = min(p$data$UMAP_1), yend = min(p$data$UMAP_2) + 4
        ),
        colour = "black", size = 0.15, arrow = arrow(length = unit(0.3, "cm"), type = "closed")
    ) +
    annotate("text",
        x = min(p$data$UMAP_1) + 2, y = min(p$data$UMAP_2) - 1, label = "UMAP1",
        color = "black", size = 6, fontface = "bold"
    ) +
    annotate("text",
        x = min(p$data$UMAP_1) - 1, y = min(p$data$UMAP_2) + 2, label = "UMAP2",
        color = "black", size = 6, fontface = "bold", angle = 90
    )

png("./Image_Seurat_DimPlot_TCR.png", height = 2000, width = 2000, res = 300)
p6
dev.off()

png("./Image_Seurat_DimPlot_TCR_NoLegend.png", height = 2000, width = 2000, res = 300)
p6 + theme(legend.position = "none")
dev.off()

# 分开每个病人进行绘图
png("./大群UMAP_splitByPatient.png", height =5000,width = 6000,res=300)
DimPlot(scRNA_seurat,group.by = "celltype",split.by = "Rename",raster=F,label = F,ncol = 6)+
    theme(legend.position = "bottom")
dev.off()



## 计算比例
scRNA_seurat$patients <- paste0(scRNA_seurat$Patient, "--", scRNA_seurat$Treat_assess, "--", scRNA_seurat$Tissue)
Proportion_analys <- function(data = scRNA_sub, IDs = "patients", celltype = "celltype") { # nolint
    ids <- table(data@meta.data[, IDs]) %>% as.data.frame() # 计算每个样本的全部细胞数目
    # ids <- table(scRNA_seurat[,celltype=="NK & T cells"]@meta.data[,IDs]) %>% as.data.frame()
    df.plot <- table(data@meta.data[, c(IDs, celltype)]) %>% as.data.frame()
    for (i in 1:dim(ids)[1]) {
        prop.temp <- df.plot[which(df.plot$patients %in% ids[i, 1]), "Freq"] / ids[i, 2]
        df.plot[which(df.plot$patients %in% ids[i, 1]), "SampleCellNum"] <- ids[i, 2]
        df.plot[which(df.plot$patients %in% ids[i, 1]), "prop"] <- prop.temp
        df.plot <<- df.plot
    }
}
Proportion_analys(data = scRNA_seurat)

# 添加信息
df.plot$group <- str_split_fixed(df.plot$patients, "--", n = 3)[, 2]
df.plot$Tissue <- str_split_fixed(df.plot$patients, "--", n = 3)[, 3]
df.plot$patients <- str_split_fixed(df.plot$patients, "--", n = 3)[, 1]
df.plot$Patient <- df.plot$patients
df.plot$Patient <- str_split_fixed(df.plot$Patient, "_", n = 2)[, 1]
df.plot$Patient <- gsub("FHXHBS2", "FHXHBS1", df.plot$Patient)
df.plot$Patient <- gsub("HEJX", "HEJI", df.plot$Patient)
df.plot$group <- factor(df.plot$group, levels = c("R_Pre", "R_Post", "NR_Pre", "NR_Post"))

# 绘图
compaired <- list(c("R_Pre", "R_Post"), c("R_Pre", "NR_Pre"), c("R_Post", "NR_Post"), c("NR_Pre", "NR_Post"))
df.plot$text <- paste0(df.plot$Patient, ":", round(df.plot$prop, digits = 4))
symnum.args <- list(cutpoints = c(0, 0.005, 0.01, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns"))
plist <- list()
for (celltype in unique(df.plot$celltype)) {
    tmp <- df.plot[df.plot$celltype == celltype, ]
    p <- ggplot(tmp, aes(x = group, y = prop)) +
        geom_boxplot(aes(fill = group), show.legend = F, width = 0.6, outlier.colour = NA) +
        geom_jitter(width = 0.1, size = 2, fill = "black", color = "black") + # 绘制散点
        geom_line(aes(group = Patient), color = "darkgrey", lwd = 0.5) + # 配对样本间连线
        ggrepel::geom_text_repel(aes(label = patients), label.size = 0.05) +
        theme(
            panel.grid = element_blank(),
            legend.position = "none",
            axis.line = element_line(colour = "black", size = 1.2),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 23, face = "bold"),
            plot.subtitle = element_text(size = 23, hjust = 0.5),
            axis.text = element_text(size = 23, color = "black"),
            axis.title = element_text(size = 23, color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        ) +
        labs(x = "", y = "Fraction", title = celltype) +
        # ggpubr::stat_compare_means(method = "wilcox.test", hide.ns = F, comparisons = compaired, label = "p.format", symnum.args = symnum.args, size = 6) +
        ggpubr::stat_compare_means(method = "wilcox.test", hide.ns = F, comparisons = compaired, label = "p.format") +
        scale_fill_manual(values = c("#3b9aac", "#a74947", "#4470a1", "#8da24f", "#6d5a87"))
    plist[[celltype]] <- p
}

# 保存图片
png("./刘老师大群细胞比例_数字.png", height = 3500, width = 4500, res = 300)
p_freq <- wrap_plots(plist, bycol = T, ncol = 3)
p_freq
dev.off()

# -------------------------------- 20230206 大群标肝细胞和内皮细胞 -----------------------------------
marker.list <- list(
    Hepatocytes = Hepatocytes,
    Epithelial = Epithelial,
    Keratinocytes = Keratinocytes
)

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
png("/root/wangje/Project/刘老师/大群/大群CCA_肝细胞_角质细胞和内皮marker.png", height = 2000, width = 9000, res = 300)
p1 <- plotFeature(scRNA_data = scRNA_seurat, choose = "Feature", col_num = 9, marker.list = marker.list)
p2 <- plotFeature(scRNA_data = scRNA_seurat, choose = "SCpubr", col_num = 9, marker.list = marker.list)
p1 / p2
dev.off()
# 查看细胞比例
scRNA_seurat$Treat_assessment <- factor(scRNA_seurat$Treat_assess,levels = c('R_Pre','R_Post','NR_Pre','NR_Post'))    
cellnum <- scRNA_seurat@meta.data %>%
    group_by(celltype) %>%
    count(Treat_assess)
write.table(cellnum,file="/root/wangje/Project/刘老师/大群/大群细胞量.txt",sep="\t",row.names=F,quote=F)

## ---------------- 大群计算比例，在原有基础上增加3Her2数据 ---------------
scRNA_sub <- scRNA_seurat
scRNA_sub$patients <- paste0(scRNA_sub$Patient,"--",scRNA_sub$Treat_assess,"--",scRNA_sub$Tissue)
Proportion_analys <- function(scRNA_seurat, IDs = "patients", celltype = "celltype") { # nolint
    ids <- table(scRNA_seurat@meta.data[, IDs]) %>% as.data.frame()
    df.plot <- table(scRNA_seurat@meta.data[, c(IDs, celltype)]) %>% as.data.frame()
    for (i in 1:dim(ids)[1]) {
        prop.temp <- df.plot[which(df.plot$patients %in% ids[i, 1]), "Freq"] / ids[i,2]
        df.plot[which(df.plot$patients %in% ids[i, 1]),"SampleCellNum"] <- ids[i,2]
        df.plot[which(df.plot$patients %in% ids[i, 1]), "prop"] <- prop.temp
        df.plot <<- df.plot
    }
}
# 计算刘老师大群的比例
scRNA_sub %>% Proportion_analys()
df.plot <- df.plot %>%
    as_tibble() %>%
    mutate(
        sample = str_split_fixed(df.plot$patients, "--", n = 3)[, 1],
        group = str_split_fixed(df.plot$patients, "--", n = 3)[, 2],
        tissue <- str_split_fixed(df.plot$patients, "--", n = 3)[, 3]
    )
Liu_big <- df.plot
df.plot <- NULL

# ---------------------- 计算3HER2的大群比例 -------------------------
# 读入大群数据
scRNA_seurat@meta.data[which(scRNA_seurat$integrated_snn_res.0.4 %in% c(8)), "celltype"] <- "B & Plasma cells"
scRNA_seurat@meta.data[which(scRNA_seurat$integrated_snn_res.0.4 %in% c(10)), "celltype"] <- "Endothelial cells"
scRNA_seurat@meta.data[which(scRNA_seurat$integrated_snn_res.0.4 %in% c(3, 5, 11)), "celltype"] <- "Epithelial cells"
scRNA_seurat@meta.data[which(scRNA_seurat$integrated_snn_res.0.4 %in% c(9)), "celltype"] <- "Fibroblasts"
scRNA_seurat@meta.data[which(scRNA_seurat$integrated_snn_res.0.4 %in% c(6, 12)), "celltype"] <- "Myeloid cells"
scRNA_seurat@meta.data[which(scRNA_seurat$integrated_snn_res.0.4 %in% c(0, 1, 2, 4, 7)), "celltype"] <- "NK & T cells"
save(scRNA_seurat,file="/root/wangje/文献数据/5HER2/Wu_etal_2021_BRCA_scRNASeq/Wu/重新分析3HER2/bigcluster_3her2_cca.rdata")
load('/root/wangje/文献数据/5HER2/Wu_etal_2021_BRCA_scRNASeq/Wu/重新分析3HER2/bigcluster_3her2_cca.rdata')
## 计算3HER2的比例
scRNA_seurat$Patient <- "3HER2"
scRNA_seurat$Treat_assess <- "TRNaive"
scRNA_seurat$Tissue <- 'Breast'
scRNA_seurat$patients <- paste0(scRNA_seurat$Patient,"--",scRNA_seurat$Treat_assess,"--",scRNA_seurat$Tissue)
scRNA_seurat %>% Proportion_analys()
df.plot <- df.plot %>%
    as_tibble() %>%
    mutate(
        sample = str_split_fixed(df.plot$patients, "--", n = 3)[, 1],
        group = str_split_fixed(df.plot$patients, "--", n = 3)[, 2],
        tissue <- str_split_fixed(df.plot$patients, "--", n = 3)[, 3]
    )
# 3HER2数据
write.table(df.plot,file="/root/wangje/Project/刘老师/大群/结合3HER2数据进行分析/3HER2_大群CCA以每个样本的不同celltype数目计算比例.txt",quote=F,sep="\t",row.names = F)
# Liu 大群数据
write.table(Liu_big,file="/root/wangje/Project/刘老师/大群/结合3HER2数据进行分析/Liu_大群CCA以每个样本的不同celltype数目计算比例.txt",quote = F,sep = "\t",row.names = F)
# 合并两个数据的比例
df.merge <- rbind(Liu_big,df.plot)

