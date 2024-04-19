#### 对结果进行绘图
# 调用需要的R包
# ------------
load_packages_silently <- function(package_names) {
  not_installed <- character(0)

  for (pkg in package_names) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste0("Package '", pkg, "' is not installed or loaded properly. Skipping..."))
      not_installed <- c(not_installed, pkg)
    } else {
      suppressPackageStartupMessages({
        library(pkg, character.only = TRUE)
      })
      message(paste0("Package '", pkg, "' successfully loaded."))
    }
  }

  # 返回未成功加载的包列表
  return(not_installed)
}

package_list <- c("dplyr", "Seurat", "SeuratObject", "ggplot2", "getopt", "patchwork", "ggpubr", "qs", "rstatix", "SCP", "data.table")
load_packages_silently(package_list)

# 读入数据
setwd("/root/wangje/Project/刘老师/大群/new_Result/Data")
scRNA <- qread("cca_大群_new.qs")

## 颜色
library(RColorBrewer)
cols <- RColorBrewer::brewer.pal(n = 12, name = "Paired")

# 1) 绘制UMAP结果
# -------------
## celltype
reduction <- "umap"
file_prefix <- "大群"
p1 <- DimPlot(scRNA, group.by = "celltype", cols = cols, reduction = reduction, raster = TRUE) +
  ggtitle(label = paste0("nCells:", dim(scRNA)[2])) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  SCP::theme_blank(xlab = "UMAP1", ylab = "UMAP2")
ggsave(filename = paste0(file_prefix, "_celltype_DimPlot01.png"), height = 4, width = 5, plot = p1, bg = "white", limitsize = F)
ggsave(filename = paste0(file_prefix, "_celltype_DimPlot01_NoLegend.png"), height = 4, width = 4, plot = p1 + Seurat::NoLegend(), bg = "white", limitsize = F)


source("/root/wangje/Project/刘老师/script/AddInformation.R")
scRNA <- AddInfo(scRNA)
scRNA$Treatment <- scRNA$Treat_assess
scRNA$sample <- stringr::str_split_fixed(scRNA$new_Rename, "-", n = 2)[, 1]
scRNA$sample <- factor(scRNA$sample, levels = paste0("P", 1:17))
scRNA$sample <- forcats::fct_drop(scRNA$sample)
scRNA$Treatment1 <- stringr::str_split_fixed(scRNA$new_Rename, "-", n = 2)[, 2]
scRNA$Treatment1 <- factor(scRNA$Treatment1, levels = c("Pre", "Post"))
scRNA$Treatment
scRNA$Treatment <- factor(scRNA$Treatment, levels = c("R_Pre", "R_Post", "NR_Pre", "NR_Post"))

## celltype
# ----------------------
# file_prefix <- "大群"
reduction <- 'umap'
cols <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
file_prefix <- "Epithelial_FastMNN"
p1 <- DimPlot(scRNA, group.by = "celltype", cols = cols, reduction = reduction, raster = TRUE) +
  ggtitle(label = paste0("nCells:", dim(scRNA)[2])) +
  guides(color = guide_legend(override.aes = list(size = 6), ncol = 1)) +
  SCP::theme_blank(xlab = "UMAP1", ylab = "UMAP2") + theme(legend.title = element_blank())
ggsave(filename = paste0(file_prefix, "_celltype_DimPlot01.png"), height = 4, width = 5, plot = p1, bg = "white", limitsize = F)
ggsave(filename = paste0(file_prefix, "_celltype_DimPlot01_NoLegend.png"), height = 4, width = 4, plot = p1 + Seurat::NoLegend(), bg = "white", limitsize = F)
ggsave(filename = paste0(file_prefix, "_celltype_DimPlot01_NoLegend.pdf"), height = 4, width = 4, plot = p1 + Seurat::NoLegend(), bg = "white", limitsize = F)
ggsave(filename = paste0(file_prefix, "_celltype_DimPlot01_NoLegend.pdf"), height = 4, width = 5, plot = p1 , bg = "white", limitsize = F)

p1 <- DimPlot(scRNA, group.by = "celltype", reduction = reduction, raster = TRUE, label = T, repel = T, cols = cols) +
  ggtitle(label = paste0("nCells:", dim(scRNA)[2])) +
  guides(color = guide_legend(override.aes = list(size = 6), ncol = 1)) +
  SCP::theme_blank(xlab = "UMAP1", ylab = "UMAP2") + theme(legend.title = element_blank())
ggsave(filename = paste0(file_prefix, "_celltype_DimPlot02.png"), height = 4, width = 5, plot = p1, bg = "white", limitsize = F)

## sample
# ---------------------
p1 <- SCP::CellDimPlot(scRNA, group.by = c("sample"), reduction = reduction, raster = TRUE, ncol = 3, show_stat = FALSE) +
  # ggtitle(label = paste0("nCells:", dim(scRNA)[2])) +
  guides(color = guide_legend(override.aes = list(size = 6)), col = 5) +
  SCP::theme_blank(xlab = "UMAP1", ylab = "UMAP2") + theme(legend.title = element_blank())
ggsave(filename = paste0(file_prefix, "_sample_DimPlot02.png"), height = 4, width = 7, plot = p1, bg = "white", limitsize = F)
ggsave(filename = paste0(file_prefix, "_sample_DimPlot02_NoLegend.png"), height = 4, width = 4, plot = p1 & Seurat::NoLegend(), bg = "white", limitsize = F)
ggsave(filename = paste0(file_prefix, "_sample_DimPlot02_NoLegend.pdf"), height = 4, width = 4, plot = p1 & Seurat::NoLegend(), bg = "white", limitsize = F)

p1 <- SCP::CellDimPlot(scRNA, group.by = "sample", reduction = reduction, raster = TRUE, show_stat = FALSE) +
  ggtitle(label = paste0("nCells:", dim(scRNA)[2])) +
  guides(color = guide_legend(override.aes = list(size = 6), ncol = 1)) +
  SCP::theme_blank(xlab = "UMAP1", ylab = "UMAP2") + theme(legend.title = element_blank())
ggsave(filename = paste0(file_prefix, "_sample_DimPlot02.png"), height = 4, width = 5, plot = p1, bg = "white", limitsize = F)

## Treat_assess
# ----------------------
scRNA$Treatment <- scRNA$Treat_assess
p1 <- SCP::CellDimPlot(scRNA, group.by = c("Treatment"), reduction = reduction, raster = TRUE, ncol = 3, show_stat = FALSE) +
  # ggtitle(label = paste0("nCells:", dim(scRNA)[2])) +
  guides(color = guide_legend(override.aes = list(size = 6)), ncol = 5) +
  SCP::theme_blank(xlab = "UMAP1", ylab = "UMAP2") + theme(legend.title = element_blank())
ggsave(filename = paste0(file_prefix, "_treatment_DimPlot02.png"), height = 4, width = 7, plot = p1, bg = "white", limitsize = F)
ggsave(filename = paste0(file_prefix, "_treatment_DimPlot02_NoLegend.png"), height = 4, width = 4, plot = p1 & Seurat::NoLegend(), bg = "white", limitsize = F)
ggsave(filename = paste0(file_prefix, "_treatment_DimPlot02_NoLegend.pdf"), height = 4, width = 4, plot = p1 & Seurat::NoLegend(), bg = "white", limitsize = F)

p1 <- SCP::CellDimPlot(scRNA, group.by = "Treatment", reduction = reduction, raster = TRUE, show_stat = FALSE) +
  ggtitle(label = paste0("nCells:", dim(scRNA)[2])) +
  guides(color = guide_legend(override.aes = list(size = 6), ncol = 1)) +
  SCP::theme_blank(xlab = "UMAP1", ylab = "UMAP2") + theme(legend.title = element_blank())
ggsave(filename = paste0(file_prefix, "_treatment_DimPlot02.png"), height = 4, width = 5, plot = p1, bg = "white", limitsize = F)

## Tissue
# -----------------------
p1 <- SCP::CellDimPlot(scRNA, group.by = c("Tissue"), reduction = reduction, raster = TRUE, ncol = 3, show_stat = FALSE) +
  # ggtitle(label = paste0("nCells:", dim(scRNA)[2])) +
  guides(color = guide_legend(override.aes = list(size = 6)), ncol = 5) +
  SCP::theme_blank(xlab = "UMAP1", ylab = "UMAP2") + theme(legend.title = element_blank())
ggsave(filename = paste0(file_prefix, "_Tissue_DimPlot02.png"), height = 4, width = 7, plot = p1, bg = "white", limitsize = F)
ggsave(filename = paste0(file_prefix, "_Tissue_DimPlot02_NoLegend.png"), height = 4, width = 4, plot = p1 & Seurat::NoLegend(), bg = "white", limitsize = F)
ggsave(filename = paste0(file_prefix, "_Tissue_DimPlot02_NoLegend.pdf"), height = 4, width = 4, plot = p1 & Seurat::NoLegend(), bg = "white", limitsize = F)

p1 <- SCP::CellDimPlot(scRNA, group.by = "Tissue", reduction = reduction, raster = TRUE, show_stat = FALSE) +
  ggtitle(label = paste0("nCells:", dim(scRNA)[2])) +
  guides(color = guide_legend(override.aes = list(size = 6), ncol = 1)) +
  SCP::theme_blank(xlab = "UMAP1", ylab = "UMAP2") + theme(legend.title = element_blank())
ggsave(filename = paste0(file_prefix, "_Tissue_DimPlot02.png"), height = 4, width = 5, plot = p1, bg = "white", limitsize = F)

## split Treatment
# -----------------------
p1 <- SCP::CellDimPlot(scRNA, group.by = c("celltype"), split.by = "Treatment", reduction = reduction, raster = FALSE, show_stat = FALSE, ncol = 4, bg_color = "white") &
  # ggtitle(label = paste0("nCells:", dim(scRNA)[2])) +
  guides(color = guide_legend(override.aes = list(size = 6)), ncol = 4) &
  SCP::theme_blank(xlab = "UMAP1", ylab = "UMAP2") & theme(legend.title = element_blank())
ggsave(filename = paste0(file_prefix, "_treatment_DimPlot03.png"), height = 4, width = 17, plot = p1, bg = "white", limitsize = F)
ggsave(filename = paste0(file_prefix, "_treatment_DimPlot03_NoLegend.png"), height = 4, width = 16, plot = p1 & Seurat::NoLegend(), bg = "white", limitsize = F)
ggsave(filename = paste0(file_prefix, "_treatment_DimPlot03_NoLegend.pdf"), height = 4, width = 16, plot = p1 & Seurat::NoLegend(), bg = "white", limitsize = F)

p1 <- SCP::CellDimPlot(scRNA, group.by = "celltype", split.by = "Treatment", reduction = reduction, raster = FALSE, show_stat = FALSE, ncol = 4, bg_color = "white") &
  ggtitle(label = paste0("nCells:", dim(scRNA)[2])) &
  guides(color = guide_legend(override.aes = list(size = 6), ncol = 4)) &
  SCP::theme_blank(xlab = "UMAP1", ylab = "UMAP2") + theme(legend.title = element_blank())
ggsave(filename = paste0(file_prefix, "_treatment_DimPlot03.png"), height = 4, width = 17, plot = p1, bg = "white", limitsize = F)

## split sample
# -----------------------
scRNA$sample <- forcats::fct_drop(scRNA$sample)
scRNA$new_Rename <- forcats::fct_drop(scRNA$new_Rename)
p1 <- SCP::CellDimPlot(scRNA, group.by = c("celltype"), split.by = "new_Rename", reduction = reduction, raster = FALSE, show_stat = FALSE, ncol = 5, bg_color = "white") &
  # ggtitle(label = paste0("nCells:", dim(scRNA)[2])) +
  guides(color = guide_legend(override.aes = list(size = 6)), ncol = 4) &
  SCP::theme_blank(xlab = "UMAP1", ylab = "UMAP2") & theme(legend.title = element_blank())
ggsave(filename = paste0(file_prefix, "_sample_DimPlot03.png"), height = 20, width = 20, plot = p1, bg = "white", limitsize = F)
ggsave(filename = paste0(file_prefix, "_sample_DimPlot03_NoLegend.png"), height = 20, width = 20, plot = p1 & Seurat::NoLegend(), bg = "white", limitsize = F)
ggsave(filename = paste0(file_prefix, "_sample_DimPlot03_NoLegend.pdf"), height = 20, width = 20, plot = p1 & Seurat::NoLegend(), bg = "white", limitsize = F)

p1 <- SCP::CellDimPlot(scRNA, group.by = "celltype", split.by = "sample", reduction = reduction, raster = FALSE, show_stat = FALSE, ncol = 5, bg_color = "white") &
  ggtitle(label = paste0("nCells:", dim(scRNA)[2])) &
  guides(color = guide_legend(override.aes = list(size = 6), ncol = 1)) &
  SCP::theme_blank(xlab = "UMAP1", ylab = "UMAP2") + theme(legend.title = element_blank())
ggsave(filename = paste0(file_prefix, "_sample_DimPlot03.png"), height = 20, width = 20, plot = p1, bg = "white", limitsize = F)

## 热图分析
DefaultAssay(scRNA) <- "RNA"
Seurat::Idents(scRNA) <- "celltype"
# scRNA <- Seurat::ScaleData(scRNA, vars.to.regress='percent.mt')
scRNA <- Seurat::ScaleData(scRNA, features = rownames(scRNA))
scRNA.Findmarkers <- FindAllMarkers(
  scRNA,
  only.pos = TRUE,
  min.pct = 0.15,
  logfc.threshold = 0.15,
  BPPARAM = MulticoreParam(20)
)

Top10.genes <- scRNA.Findmarkers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC)
write.csv(Top10.genes, file = "FastMNNO2_Top30基因.csv", quote = F, row.names = T)
write.csv(scRNA.Findmarkers, file = "FastMNNO2_Top基因.csv", quote = F, row.names = T)


## 绘制热图
# --------------------
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
# 提取符合条件的scale.data
test <- scRNA
Top10.genes %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice(1:20) %>%
  dplyr::select(cluster, gene) -> test.gene

# test$cluster <- varhandle::unfactor(test$cluster)
# test <- split(test$gene,test$cluster)

# 提取细胞
test@meta.data %>%
  dplyr::select(celltype) %>%
  tibble::rownames_to_column(var = "cell") %>%
  dplyr::group_by(celltype) %>%
  dplyr::sample_n(500) -> select.cells

# 提取scale.data
# data.all <- scRNA@assays$RNA@scale.data[test$gene,select.cells$cell]
data.all <- test[["RNA"]]$scale.data[test.gene$gene, select.cells$cell]

## ---------------------
# 获得样本对应的sample名
source("/root/wangje/Project/刘老师/script/AddInformation.R")
test <- AddInfo(test)
test$sample <- stringr::str_split_fixed(test$new_Rename, "-", n = 2)[, 1]

# 获得sample和Patient的对应名称
test@meta.data %>%
  dplyr::select(Patient, sample) %>%
  dplyr::distinct() -> result1
sample_patient <- split(result1$sample, result1$Patient)
# 获得Patient对应得sample名
cells <- colnames(data.all)
cells.sample <- stringr::str_split_fixed(cells, "_[A|T|G|C].*", n = 2)[, 1]
cells.sample <- sapply(cells.sample, function(x) {
  sample_patient[[x]]
})
# 设置factor
cells.sample <- factor(cells.sample, levels = paste0("P", 1:17))

# 添加分组
test@meta.data %>%
  dplyr::select(Patient, Treat_assess) %>%
  dplyr::distinct() %>%
  varhandle::unfactor() -> result2
group_patient <- split(result2$Treat_assess, result2$Patient)

cells <- colnames(data.all)
cells_anno <- stringr::str_split_fixed(cells, "_[A|T|G|C].*", n = 2)[, 1]
cells.group <- sapply(cells_anno, function(x) {
  group_patient[[x]]
})
# 设置cells.group得factor
cells.group <- factor(cells.group, levels = c("R_Pre", "R_Post", "NR_Pre", "NR_Post"))

# 注释颜色
cols1 <- c("#33a02c", "#1f78b4", "#a6cee3", "#b2df8a", "#a2c9dc", "#1a71a9", "#acd388", "#329838", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00")[1:length(unique(test$celltype))]
names(cols1) <- levels(test$celltype)

# sample
cols2 <- c(
  "#33a02c", "#1f78b4", "#a6cee3", "#b2df8a", "#0f964d", "#016c4c", "#0d4090", "#5c960a", "#b37c06", "#a80541", "#049964", "#4c0b7b", "#05862a",
  "#2d7003", "#012c86", "#c10d62", "#940dcb", "#b50c44", "#0c8a0d", "#257704", "#a30426"
)[1:length(levels(cells.sample))]
names(cols2) <- levels(cells.sample)

# group
cols3 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")[1:length(levels(test$Treat_assess))]
names(cols3) <- levels(test$Treat_assess)

# col_fun1= colorRamp2(c(0,400, 10000, 40000), c('#ffd4cf',"#ffd4cf", "#ff7463", "#ff1b00"))
# cell num
cell.num <- table(test$celltype)
# colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100)
# values <- seq(-2, 2, length.out = 101)[-101]
# col_fun <- colorRamp2(values, colors)
col_fun1 <- colorRamp2(c(200, 70000), c("white", "red"))
column_ha <- HeatmapAnnotation(foo2 = cell.num, col = list(foo2 = col_fun1))


# 头部注释
top_anno <- HeatmapAnnotation(
  celltype = rep(levels(test$celltype), each = 500), group = cells.group, sample = cells.sample,
  cellNum = rep(cell.num, each = 500),
  col = list(celltype = cols1, cellNum = col_fun1, sample = cols2, group = cols3)
)

# 右侧注释
genenames <- split(test.gene$gene, test.gene$cluster)
# ha <- rowAnnotation(foo= anno_empty(border = TRUE, width = max_text_width(unlist(genenames)) + unit(4, "mm")))
ha <- rowAnnotation(foo = anno_block(gp = gpar(fill = cols1, col = "white")))

col_fun <- colorRamp2(c(-2, 0, 2), c("#f400f4", "black", "#fbfb02"))

h1 <- Heatmap(
  data.all,
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  use_raster = FALSE, # 是否栅格化，
  name = "Scaled Expression",
  column_split = rep(c(LETTERS[1:length(unique(test$celltype))]), each = 500),
  column_gap = unit(0, "mm"), border = TRUE, border_gp = gpar(col = "white"),
  row_split = rep(LETTERS[1:length(unique(test$celltype))], each = 20),
  row_gap = unit(0, "mm"),
  row_title = NULL,
  column_title = levels(test$celltype), # 设置列名
  column_title_side = "bottom",
  column_title_gp = gpar(fontsize = 20, fontface = "bold"),
  column_title_rot = 90,
  top_annotation = top_anno,
  right_annotation = ha,
  width = unit(10, "cm"),
  height = unit(9, "cm")
)
# draw(h1,heatmap_legend_side = "right")
## 修饰右侧注释的颜色
# 设置legend的位置
png("marker热图.png", height = 3000, width = 4000, res = 200)
print(draw(h1, heatmap_legend_side = "right"))
dev.off()

pdf("结合内皮上皮细胞_marker热图.pdf", height = 8, width = 11)
draw(h1, heatmap_legend_side = "right")
dev.off()

### 计算比例
# ------------------------
library(rstatix)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggExtra)

## 计算细胞比例
# scRNA$celltype <- scRNA$RNA_snn_res.0.2
df1 <- read.table("/root/wangje/Project/刘老师/细胞Fraction/CCA大群_sample_nCells.txt", sep = "\t", header = T)
df2 <- scRNA@meta.data %>%
  dplyr::select(Rename, celltype, Treat_assess) %>%
  dplyr::group_by(Rename, celltype, Treat_assess) %>%
  dplyr::count()
df <- df2 %>%
  dplyr::left_join(df1, by = c("Rename")) %>%
  dplyr::mutate(Prop = n / AllCells)
df1 <- df

lapply(unique(df1$celltype), FUN = function(x) {
  tryCatch(
    {
      df_sub <- df1[df1$celltype == x, ]
      # plot_boxplot
      p1 <- ggplot(df_sub) +
        geom_boxplot(
          aes(x = celltype, y = Prop, color = Treat_assess, fill = Treat_assess),
          position = position_dodge(width = 0.8),
          outlier.shape = NA, color = "black"
        ) +
        geom_jitter(aes(x = celltype, y = Prop, fill = Treat_assess),
          color = "black", position = position_dodge(width = 0.8), pch = 21
        )
      # 计算显著性
      stat.test <- df_sub %>%
        group_by(celltype) %>%
        wilcox_test(Prop ~ Treat_assess, comparisons = list(c("R_Pre", "R_Post"), c("NR_Pre", "NR_Post"), c("R_Pre", "NR_Pre"), c("R_Post", "NR_Post"))) %>%
        add_xy_position(x = "celltype", dodge = 0.8)

      # 添加新的显著性
      stat.test$new_signif <- case_when(
        0 <= stat.test$p & stat.test$p < 0.01 ~ "***",
        0.01 <= stat.test$p & stat.test$p < 0.05 ~ "**",
        0.05 <= stat.test$p & stat.test$p < 0.1 ~ "*",
        TRUE ~ "ns"
      )
      stat.test <- stat.test %>% filter(new_signif != "ns")
      p2 <- p1 + stat_pvalue_manual(
        stat.test,
        label = "new_signif", tip.length = 0.00, size = 8,
        hide.ns = FALSE
      ) + labs(y = "Fraction", x = "")
      p2 <- p2 +
        theme_classic(base_size = 20, base_line_size = 1) +
        scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#ff7f01", "#fb9a99")) +
        theme(
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 18, color = "black"),
          axis.text.y = element_text(size = 22, colour = "black"),
          axis.title = element_text(size = 20, colour = "black"),
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5, colour = "black", size = 20, angle = 0)
        ) +
        scale_x_discrete(expand = c(0.45, 0))
      return(p2)
    },
    error = function(e) {
      message(e$message)
      return(NULL)
    }
  )
}) -> plist

p <- ggpubr::ggarrange(plotlist = plist, ncol = 3, nrow = 1, common.legend = TRUE, align = "v", legend = "bottom")

pdf("Epithelials_FastMNN比例箱线图_Arial02_大群.pdf", width = 9, height = 4, family = "ArialMT", bg = "white")
p
dev.off()
ggpubr::ggarrange(plotlist = plist, ncol = 4, nrow = 2, common.legend = TRUE, align = "v", legend = "bottom")

png("./Epithelials_FastMNN_比例箱线图.png", height = 5, width = 13, units = "in", res = 300)
p
dev.off()


## 2) 富集分析
# --------------
library(Scillus)
library(forcats)
library(ggplot2)
library(patchwork)
library(org.Hs.eg.db)

## 寻找差异基因
Idents(scRNA) <- scRNA$celltype
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.1, min.pct = 0, only.pos = T)
# 保存文件
write.csv(markers, file = "CCA大群FindAllMarkers.csv", quote = F)

# GO富集
p.go <- lapply(unique(markers$cluster), function(x) {
  print(x)
  plot_cluster_go(markers, cluster_name = x, org = "human", ont = "BP")
})

png("./大群GO_BP富集.png", height = 18, width = 30, units = "in", res = 300)
plot_all_cluster_go(markers, org = "human", ont = "BP")
dev.off()

## 富集结果
sc_GO <- function(markers, cluster_name, topn = 100, org, prefix = "", ...) {
  pkg_name <- ifelse(org == "human", "org.Hs.eg.db", "org.Mm.eg.db")
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop(paste("Package", pkg_name, "needed for this function to work. Please install it."),
      call. = FALSE
    )
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop(paste("Package \"clusterProfiler\" needed for this function to work. Please install it."),
      call. = FALSE
    )
  }
  gene_list <- markers %>%
    filter(.data$cluster == cluster_name) %>%
    arrange(.data$p_val_adj) %>%
    head(topn) %>%
    pull(.data$gene)
  db <- if (org == "human") {
    org.Hs.eg.db::org.Hs.eg.db
  } else {
    org.Mm.eg.db::org.Mm.eg.db
  }
  res <- clusterProfiler::enrichGO(
    gene = gene_list, OrgDb = db,
    keyType = "SYMBOL", ...
  )
  df <- as_tibble(res@result) %>%
    arrange(.data$p.adjust)
  write.csv(df, file = paste0(prefix, "_GO.csv"), quote = FALSE)
  df <- df %>%
    head(10) %>%
    mutate(cluster = cluster_name) %>%
    mutate(Description = stringr::str_to_title(.data$Description)) %>%
    mutate(Description = fct_reorder(.data$Description, dplyr::desc(.data$p.adjust)))
  ggplot(df, mapping = aes(x = .data$Description, y = -log10(.data$p.adjust))) +
    geom_bar(aes(fill = .data$Count), stat = "identity") +
    scale_fill_gradient2("Gene Count",
      low = "lightgrey",
      mid = "#feb24c", high = "#bd0026"
    ) +
    coord_flip() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    xlab("Gene Ontology") +
    ylab(bquote("-log"[10] ~ " adjusted p-value")) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    ggtitle(cluster_name)
}

p.go <- lapply(unique(markers$cluster), function(x) {
  print(x)
  sc_GO(markers, cluster_name = x, org = "human", ont = "BP", prefix = x)
})






p <- ggplot() +
  geom_jitter(
    data = dt,
    aes(x = cluster, y = avg_log2FC, color = label),
    size = 0.85,
    width = 0.4
  )
p


p <- ggplot() +
  geom_jitter(
    data = dt,
    aes(x = cluster, y = avg_log2FC, color = label),
    size = 0.85,
    width = 0.4
  ) +
  geom_jitter(
    data = df.top30,
    aes(x = cluster, y = avg_log2FC, color = label),
    size = 4,
    width = 0.4
  )

## 绘制细胞比例柱形图 ---------------------------------------------------------------------------------------------
library(Seurat)
library(patchwork)
library(magrittr)
library(ggplot2)
library(dplyr)
library(SCP)
library(ggalluvial)
library(COSG)
library(stringr)
library(qs)

## 大群
# -------
setwd("/root/wangje/Project/刘老师/大群/new_Result/Data/IMG")
scRNA <- qread("../cca_大群_new.qs")
source("/root/wangje/Project/刘老师/script/AddInformation.R")
scRNA %<>% AddInfo()
scRNA@meta.data %<>% dplyr::mutate(sample = fct_drop(factor(str_split_fixed(new_Rename, "-", n = 2)[, 1], levels = paste0("P", 1:17))))
## 计算细胞比例
df1 <- read.table("/root/wangje/Project/刘老师/细胞Fraction/CCA大群_sample_nCells.txt", sep = "\t", header = T)
df2 <- scRNA@meta.data %>%
  dplyr::select(new_Rename, celltype, Treat_assess) %>%
  dplyr::group_by(new_Rename, celltype, Treat_assess) %>%
  dplyr::count() %>%
  dplyr::left_join(df1, by = c("new_Rename")) %>%
  dplyr::mutate(Fraction = n / AllCells)
write.csv(df2, "CCA大群细胞比例.csv", quote = F) # 保存文件

## 绘制柱形图
## 虚线连接
# ----------------------
p1 <- df2 %>%
  ggplot(aes(x = new_Rename, y = Fraction, fill = celltype, stratum = celltype, alluvium = celltype)) +
  geom_flow(width = 0.7, curve_type = "linear", alpha = 0.3) +
  geom_stratum(width = 0.7, alpha = 0.9, color = "black") +
  geom_alluvium(width = 0.7, curve_type = "linear", fill = NA, color = "black", linetype = "dashed") +
  scale_y_continuous(
    expand = c(0, 0),
    name = NULL
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)), ncol = length(unique(df2$celltype))) +
  labs(title = NULL) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(12, name = "Paired")) +
  theme_test(base_line_size = 1.3, base_rect_size = 1.3) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 20, angle = 90, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.ticks.y = element_line(linewidth = 1.3, color = "black"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_line(linewidth = 1.3, color = "black"),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 22, angle = 90),
    legend.text = element_text(size = 20, color = "black")
  )

## 非虚线连接
# -----------------------
p2 <- df2 %>%
  ggplot(aes(x = new_Rename, y = Fraction, fill = celltype, stratum = celltype, alluvium = celltype)) +
  geom_flow(width = 0.7, curve_type = "linear", alpha = 0.3) +
  geom_stratum(width = 0.7, alpha = 0.9, color = "black") +
  geom_alluvium(width = 0.7, curve_type = "linear", fill = NA, color = "black") +
  scale_y_continuous(
    expand = c(0, 0),
    name = NULL
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)), ncol = length(unique(df2$celltype))) +
  labs(title = NULL) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(12, name = "Paired")) +
  theme_test(base_line_size = 1.3, base_rect_size = 1.3) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 20, angle = 90, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.ticks.y = element_line(linewidth = 1.3, color = "black"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_line(linewidth = 1.3, color = "black"),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 22, angle = 90),
    legend.text = element_text(size = 20, color = "black")
  )

## 拼图
# -------------------------
library(dplyr)
library(ggplot2)
p3 <- df2 %>% ggplot() +
  geom_bar(aes(x = new_Rename, y = Fraction, fill = celltype), stat = "identity", width = 0.7, size = 0.5, colour = "#222222") +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(x = " ", y = "Fraction") +
  theme(
    axis.text.y = element_text(size = 15, color = "black"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks = element_line(linewidth = 1.3),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    legend.text = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 18, colour = "black")
  ) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(12, name = "Paired"))

num <- as.data.frame(table(scRNA$new_Rename))
p4 <- ggplot(num, aes(x = Var1, y = Freq)) +
  geom_bar(color = "black", fill = "#1f77b3", width = 0.7, stat = "identity") +
  theme_classic() +
  labs(x = " ", y = "Cell Num") +
  theme(
    axis.text.y = element_text(size = 15, color = "black"),
    axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks = element_line(linewidth = 1.3),
    legend.text = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 18, colour = "black")
  )
# 合并图片
p3 <- p3 + p4 + patchwork::plot_layout(nrow = 2, heights = c(7, 1.5))
# 保存图片
pdf("./大群柱形图.pdf", family = "ArialMT", height = 21, width = 17)
p1 + p2 + p3 + plot_layout(ncol = 1, nrow = 3)
dev.off()

## ----------------------------------------------------------------------------------------------------------------
## 不同分组
## Tisue
# 写出样本数量文件
df3 <- scRNA@meta.data %>%
  dplyr::group_by(new_Rename, sample, Tissue) %>%
  dplyr::count()
write.table(df3, file = "/root/wangje/Project/刘老师/细胞Fraction/CCA大群_Tissue_nCells.txt",sep = '\t', quote = F)
df <- as.data.frame(prop.table(table(scRNA$Tissue, scRNA$celltype), margin = 1)) %>%
  dplyr::rename("Tissue" = "Var1", "celltype" = "Var2", "Freq" = "Freq") %>%
  mutate(Fraction = 100 * .$Freq)

p1 <- df %>%
  ggplot(aes(x = Tissue, y = Fraction, fill = celltype, stratum = celltype, alluvium = celltype)) +
  geom_flow(width = 0.7, curve_type = "linear", alpha = 0.3) +
  geom_stratum(width = 0.7, alpha = 0.9, color = "black") +
  geom_alluvium(width = 0.7, curve_type = "linear", fill = NA, color = "black", linetype = "dashed") +
  scale_y_continuous(
    expand = c(0, 0),
    name = NULL
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)), ncol = length(unique(df2$celltype))) +
  labs(title = NULL, y = "Fraction (%)") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(12, name = "Paired")) +
  theme_test(base_line_size = 1.3, base_rect_size = 1.3) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 20, angle = 90, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.ticks.y = element_line(linewidth = 1.3, color = "black"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_line(linewidth = 1.3, color = "black"),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 22, angle = 90),
    legend.text = element_text(size = 20, color = "black")
  )

## 非虚线连接
# -----------------------
p2 <- df %>%
  ggplot(aes(x = Tissue, y = Fraction, fill = celltype, stratum = celltype, alluvium = celltype)) +
  geom_flow(width = 0.7, curve_type = "linear", alpha = 0.3) +
  geom_stratum(width = 0.7, alpha = 0.9, color = "black") +
  geom_alluvium(width = 0.7, curve_type = "linear", fill = NA, color = "black") +
  scale_y_continuous(
    expand = c(0, 0),
    name = NULL
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)), ncol = length(unique(df2$celltype))) +
  labs(title = NULL) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(12, name = "Paired")) +
  theme_test(base_line_size = 1.3, base_rect_size = 1.3) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 20, angle = 90, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.ticks.y = element_line(linewidth = 1.3, color = "black"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_line(linewidth = 1.3, color = "black"),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 22, angle = 90),
    legend.text = element_text(size = 20, color = "black")
  )

## 拼图
# -------------------------
library(dplyr)
library(ggplot2)
p3 <- df %>% ggplot() +
  geom_bar(aes(x = Tissue, y = Fraction, fill = celltype), stat = "identity", width = 0.7, size = 0.5, colour = "#222222") +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(x = " ", y = "Fraction") +
  theme(
    axis.text.y = element_text(size = 15, color = "black"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks = element_line(linewidth = 1.3),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    legend.text = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 18, colour = "black")
  ) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(12, name = "Paired"))

num <- as.data.frame(table(scRNA$Tissue))
p4 <- ggplot(num, aes(x = Var1, y = Freq)) +
  geom_bar(color = "black", fill = "#1f77b3", width = 0.7, stat = "identity") +
  theme_classic() +
  labs(x = " ", y = "Cell Num") +
  theme(
    axis.text.y = element_text(size = 15, color = "black"),
    axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks = element_line(linewidth = 1.3),
    legend.text = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 18, colour = "black")
  )
# 合并图片
p3 <- p3 + p4 + patchwork::plot_layout(nrow = 2, heights = c(7, 1.5))
# 保存图片
pdf("./大群柱形图_tissue.pdf", family = "ArialMT", height = 18, width = 6)
p1 + p2 + p3 + plot_layout(ncol = 1, nrow = 3)
dev.off()

## ----------------------------------------------------------------------------------------------------------------
## 不同分组
## Treatment
# 写出样本数量文件
df3 <- scRNA@meta.data %>%
  dplyr::group_by(new_Rename, sample, Treat_assess) %>%
  dplyr::count()
write.table(df3, file = "/root/wangje/Project/刘老师/细胞Fraction/CCA大群_Treat_assess_nCells.txt",sep = '\t', quote = F)
df <- as.data.frame(prop.table(table(scRNA$Treat_assess, scRNA$celltype), margin = 1)) %>%
  dplyr::rename("Treat_assess" = "Var1", "celltype" = "Var2", "Freq" = "Freq") %>%
  mutate(Fraction = 100 * .$Freq)

p1 <- df %>%
  ggplot(aes(x = Treat_assess, y = Fraction, fill = celltype, stratum = celltype, alluvium = celltype)) +
  geom_flow(width = 0.7, curve_type = "linear", alpha = 0.3) +
  geom_stratum(width = 0.7, alpha = 0.9, color = "black") +
  geom_alluvium(width = 0.7, curve_type = "linear", fill = NA, color = "black", linetype = "dashed") +
  scale_y_continuous(
    expand = c(0, 0),
    name = NULL
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)), ncol = length(unique(df2$celltype))) +
  labs(title = NULL, y = "Fraction (%)") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(12, name = "Paired")) +
  theme_test(base_line_size = 1.3, base_rect_size = 1.3) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 20, angle = 90, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.ticks.y = element_line(linewidth = 1.3, color = "black"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_line(linewidth = 1.3, color = "black"),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 22, angle = 90),
    legend.text = element_text(size = 20, color = "black")
  )

## 非虚线连接
# -----------------------
p2 <- df %>%
  ggplot(aes(x = Treat_assess, y = Fraction, fill = celltype, stratum = celltype, alluvium = celltype)) +
  geom_flow(width = 0.7, curve_type = "linear", alpha = 0.3) +
  geom_stratum(width = 0.7, alpha = 0.9, color = "black") +
  geom_alluvium(width = 0.7, curve_type = "linear", fill = NA, color = "black") +
  scale_y_continuous(
    expand = c(0, 0),
    name = NULL
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)), ncol = length(unique(df2$celltype))) +
  labs(title = NULL) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(12, name = "Paired")) +
  theme_test(base_line_size = 1.3, base_rect_size = 1.3) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 20, angle = 90, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.ticks.y = element_line(linewidth = 1.3, color = "black"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_line(linewidth = 1.3, color = "black"),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 22, angle = 90),
    legend.text = element_text(size = 20, color = "black")
  )

## 拼图
# -------------------------
library(dplyr)
library(ggplot2)
p3 <- df %>% ggplot() +
  geom_bar(aes(x = Treat_assess, y = Fraction, fill = celltype), stat = "identity", width = 0.7, size = 0.5, colour = "#222222") +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(x = " ", y = "Fraction") +
  theme(
    axis.text.y = element_text(size = 15, color = "black"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks = element_line(linewidth = 1.3),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    legend.text = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 18, colour = "black")
  ) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(12, name = "Paired"))

num <- as.data.frame(table(scRNA$Treat_assess))
p4 <- ggplot(num, aes(x = Var1, y = Freq)) +
  geom_bar(color = "black", fill = "#1f77b3", width = 0.7, stat = "identity") +
  theme_classic() +
  labs(x = " ", y = "Cell Num") +
  theme(
    axis.text.y = element_text(size = 15, color = "black"),
    axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks = element_line(linewidth = 1.3),
    legend.text = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 18, colour = "black")
  )
# 合并图片
p3 <- p3 + p4 + patchwork::plot_layout(nrow = 2, heights = c(7, 1.5))
# 保存图片
pdf("./大群柱形图_Treat_assess.pdf", family = "ArialMT", height = 6, width = 18)
p1 + p2 + p3 + plot_layout(ncol = 3, nrow = 1)
dev.off()

## -------------------------------------------------------------------------------------------------------------------------
## 拟时间分析
## -------------------------------------------------------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(SCP)
library(qs)
library(patchwork)

setwd("/root/wangje/Project/刘老师/NK_T/new_Result/Data")
scRNA <- qread('./new_CCA新分群结果.qs')
# 添加信息
source("/root/wangje/Project/刘老师/script/AddInformation.R")
scRNA  <-  AddInfo(scRNA)
scRNA$Treatment  <-  scRNA$Treat_assess
scRNA$sample  <- stringr::str_split_fixed(scRNA$new_Rename, "-", n  =  2)[, 1]
scRNA$sample  <-  factor(scRNA$sample, levels  =  paste0("P", 1:17))
scRNA$sample  <- forcats::fct_drop(scRNA$sample)
scRNA$Treatment1  <- stringr::str_split_fixed(scRNA$new_Rename, "-", n  =  2)[, 2]
scRNA$Treatment1  <-  factor(scRNA$Treatment1, levels  =  c("Pre", "Post"))
scRNA$Treatment
scRNA$Treatment  <-  factor(scRNA$Treatment, levels  =  c("R_Pre", "R_Post", "NR_Pre", "NR_Post"))
# 绘图
## celltype
cols  <- RColorBrewer::brewer.pal(n  =  12, name  =  "Paired")
reduction <- 'umap'
p1  <-  DimPlot(scRNA, 
    group.by  =  "celltype", 
    cols  =  cols, 
    reduction  =  reduction,
    raster  =  TRUE) +
    ggtitle(label  =  paste0("nCells:", dim(scRNA)[2])) +
    guides(color  =  guide_legend(override.aes  =  list(size  =  6), ncol  =  1)) +
    SCP::theme_blank(xlab  =  "UMAP1", ylab  =  "UMAP2") +
    theme(legend.title  =  element_blank())

## sample
p2  <- SCP::CellDimPlot(scRNA, 
    group.by  =  c("sample"), 
    reduction  =  reduction, 
    raster  =  TRUE, 
    ncol  =  3, 
    show_stat  =  FALSE) +
 # ggtitle(label = paste0("nCells:", dim(scRNA)[2])) +
    guides(color  =  guide_legend(override.aes  =  list(size  =  6)), col  =  5) +
    SCP::theme_blank(xlab  =  "UMAP1", ylab  =  "UMAP2") +
    theme(legend.title  =  element_blank())
    
## Treatment
p3 <- DimPlot(scRNA, 
    group.by  =  "Treatment", 
    cols  =  cols, 
    reduction  =  reduction,
    raster  =  TRUE) +
    ggtitle(label  =  paste0("nCells:", dim(scRNA)[2])) +
    guides(color  =  guide_legend(override.aes  =  list(size  =  6), ncol  =  1)) +
    SCP::theme_blank(xlab  =  "UMAP1", ylab  =  "UMAP2") +
    theme(legend.title  =  element_blank())

## Tissue
p4 <- DimPlot(scRNA, 
    group.by  =  "Tissue", 
    cols  =  cols, 
    reduction  =  reduction,
    raster  =  TRUE) +
    ggtitle(label  =  paste0("nCells:", dim(scRNA)[2])) +
    guides(color  =  guide_legend(override.aes  =  list(size  =  6), ncol  =  1)) +
    SCP::theme_blank(xlab  =  "UMAP1", ylab  =  "UMAP2") +
    theme(legend.title  =  element_blank())
## 输出图片
png('./NK_T细胞DimPlot01.png', height=4, width = 18, unit='in', res=300)
p1 + p2 + p3 +p4 + plot_layout(ncol = 4)
dev.off()


scRNA_sub <- NormalizeData(scRNA_sub)
scRNA_sub <- FindVariableFeatures(scRNA_sub)
scRNA_sub <- ScaleData(scRNA_sub)
scRNA_sub <- RunPCA(scRNA_sub)

scRNA_sub <- RunUMAP(scRNA_sub,dims=1:30)
scRNA_sub <- FindNeighbors(scRNA_sub, reduction = "integrated.cca", dims = 1:30)
scRNA_sub <- FindClusters(scRNA_sub, resolution = seq(0.1,2,0.1))

## 使用SCP进行多种方法聚类
sce <- SCP::Integration_SCP(scRNA_sub, batch = 'sample',integration_method = 'BBKNN',vars_to_regress = 'percent.mt',cluster_resolution = seq(0.1,2,0.1))
qsave(sce, file = '02_SCP_CD4.qs')
# fastMNN
sce <- SCP::Integration_SCP(scRNA_sub, batch = 'sample',integration_method = 'fastMNN',vars_to_regress = 'percent.mt',cluster_resolution = seq(0.1,2,0.1))
qsave(sce, file = '02_SCP_CD4.qs')
# fastMNN
sce <- SCP::Integration_SCP(scRNA_sub, batch = 'sample',integration_method = 'Seurat',vars_to_regress = 'percent.mt',cluster_resolution = seq(0.1,2,0.1))
qsave(sce, file = '02_SCP_CD4.qs')
# Harmony
sce <- SCP::Integration_SCP(scRNA_sub, batch = 'sample',integration_method = 'Harmony',vars_to_regress = 'percent.mt',cluster_resolution = seq(0.1,2,0.1))
qsave(sce, file = '02_SCP_CD4.qs')
# scVI
sce <- SCP::Integration_SCP(scRNA_sub, batch = 'sample',integration_method = 'scVI',vars_to_regress = 'percent.mt',cluster_resolution = seq(0.1,2,0.1))
qsave(sce, file = '03_SCP_cd4.qs')




## ---------------------------------------------------------------------------------------------
## CD8 + Tcells
scRNA_sub <- scRNA[,scRNA$celltype1 %in% c('CD8+ Tn','Texp','Tex')]
sce <- SCP::Integration_SCP(scRNA_sub, batch = 'sample',integration_method = 'BBKNN',vars_to_regress = 'percent.mt',cluster_resolution = seq(0.1,2,0.1))
qsave(sce, file = '02_SCP_CD8.qs')
# fastMNN
sce <- SCP::Integration_SCP(scRNA_sub, batch = 'sample',integration_method = 'fastMNN',vars_to_regress = 'percent.mt',cluster_resolution = seq(0.1,2,0.1))
qsave(sce, file = '02_SCP_CD8.qs')
# fastMNN
sce <- SCP::Integration_SCP(scRNA_sub, batch = 'sample',integration_method = 'Seurat',vars_to_regress = 'percent.mt',cluster_resolution = seq(0.1,2,0.1))
qsave(sce, file = '02_SCP_CD8.qs')
# Harmony
sce <- SCP::Integration_SCP(scRNA_sub, batch = 'sample',integration_method = 'Harmony',vars_to_regress = 'percent.mt',cluster_resolution = seq(0.1,2,0.1))
qsave(sce, file = '02_SCP_CD8.qs')



p2 <- ggplot(data = res.all,
             aes(x = Count, y = Description)) +
  geom_point(aes(size = Count, color = -log10(pvalue))) + #气泡大小及颜色设置
  scale_color_distiller(palette = "Spectral",direction = -1) +
  labs(x = "Gene Number",
       y = "",
       title = "Dotplot of Enriched KEGG Pathways",
       size = "gene number") + #图例名
  theme_bw() 


# 多样本进行双细胞检测，需要分别对每一个样本进行分析
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize") %>% 
  FindVariableFeatures(., selection.method = "vst", nfeatures = 3000)
scRNA <- ScaleData(scRNA, vars.to.regress = c("percent.mt","nCount_RNA"), verbose = T)
scRNA <- RunPCA(scRNA, npcs = 50, verbose = FALSE)

#去批次
library(harmony)
scRNA <- RunHarmony(scRNA,"sample")
Seurat::ElbowPlot(scRNA, ndims = 50) # 查看不同PC数的分布，选择平台期的pc
scRNA <- RunUMAP(scRNA,  dims = 1:30, reduction = "harmony") # 选择前30个PC
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:40) 
scRNA  <- FindClusters(object = scRNA , resolution = seq(from = 0.1, to = 1.0,  by = 0.1))
# 保存文件
qsve(scRNA, file = "04_Harmony.qs")



df_sig  <- df[df$p_val_adj < 0.05, ]
library(clusterProfiler)
library(ggplot2)

group <- data.frame(gene=df_sig$gene,
                    group=df_sig$cluster)

Gene_ID <- bitr(df_sig$gene, fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")

#构建文件并分析
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')

data_GO <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

data_GO_sim <- simplify(data_GO, 
                        cutoff=0.7, 
                        by="p.adjust", 
                        select_fun=min)
