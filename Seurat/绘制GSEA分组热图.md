# 绘制GSEA分组热图

## 01 进行差异分析,获得差异基因

### 文章中的图片

Tumor and immune reprogramming duringimmunotherapy in advanced renal cell carcinoma（Fig.3E）

```r
# 差异分析
library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)

# 读入数据
load("/root/wangje/Project/刘老师/Myeloids/Data/Myeloids_CCA.RData")

Find_DEGs <- function(data, ident = 'new_celltype',ident.1 = "R_Post", ident.2 = "R_Pre", group = "Treat_assess") {
  flist <- list()
  if (DefaultAssay(data) != "RNA") {
    DefaultAssay(data) <- "RNA"
    Idents(data) <- ident
  }
  for (celltype in unique(Idents(data))) {
    print(celltype)
    res <- FindMarkers(
      object = data,
      ident.1 = ident.1,
      ident.2 = ident.2,
      group.by = group,
      subset.ident = celltype,
      # thresh.use = 0.99,
      logfc.threshold = 0,
      min.pct = 0,
      test.use = "MAST",
      only.pos = F
    )
    print(dim(res))
    flist[[celltype]] <- res
  }
  return(flist)
}

# 分组分析
R_Post_R_Pre <- Find_DEGs(data = scRNA_seurat, group = "Treat_assess", ident.1 = "R_Post", ident.2 = "R_Pre")
NR_Post_NR_Pre <- Find_DEGs(data = scRNA_seurat, ident.1 = "NR_Post", ident.2 = "NR_Pre", group = "Treat_assess")
NR_Post_R_Post <- Find_DEGs(data = scRNA_seurat, ident.1 = "NR_Post", ident.2 = "R_Post", group = "Treat_assess")
NR_Pre_R_Pre <- Find_DEGs(data = scRNA_seurat, ident.1 = "NR_Pre", ident.2 = "R_Pre", group = "Treat_assess")

```

## 02 进行GSEA分析

```R
# 读取GeneSet文件，并且生成geneSet列表
geneSet <- read.table('/root/wangje/Project/刘老师/Myeloids/Data/GSEA/FIG3基因集合.txt', header = T,sep = '\t')
geneSet <- geneSet %>% split(x = .$gene, f = .$term)

# 进行GSEA分析
fgsea_analysis <- function(DEGs,geneSet = geneSet, nperm = 1000){
    result <- list()
    for (name in names(DEGs)) {
        print(name)
        DEGs[[name]]$gene <- rownames(DEGs[[name]])
        geneList = DEGs[[name]] %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene, avg_log2FC)
        ranks<- deframe(geneList)
        res <- fgsea::fgsea(geneSet, stats = ranks, nperm = nperm)
        # 根据p值和q value对分组进行筛选
        res$group <- case_when(res$pval < 0.05 & res$padj < 0.25 ~ 'postivate',
                                       TRUE ~ 'ns')
        # 添加celltype
        res$celltype <- name
        # 重新设置NES值
        res$NES <- ifelse(res$group == 'ns', NA, res$NES)
        result[[name]] <- res
    }
    return(result)
}

R_Post_R_Pre_fgsea <- fgsea_analysis(DEGs = R_Post_R_Pre, geneSet = geneSet, nperm = 1000)
NR_Post_NR_Pre_fgsea <- fgsea_analysis(DEGs = NR_Post_NR_Pre, geneSet = geneSet, nperm = 1000)
NR_Post_R_Post_fgsea <- fgsea_analysis(DEGs = NR_Post_R_Post, geneSet = geneSet, nperm = 1000)
NR_Pre_R_Pre_fgsea <- fgsea_analysis(DEGs = NR_Pre_R_Pre, geneSet = geneSet, nperm = 1000)
```

## 03 整理获得的GSEA数据

```R
tidy_data <- function(fgsea_list){
    tmp_data = Reduce(rbind,fgsea_list, accumulate = F) %>%
        as.data.frame() %>% 
        select(pathway, NES, celltype) %>% 
        tidyr::pivot_wider(names_from = celltype, values_from = NES) %>% 
        as.data.frame() %>% 
        column_to_rownames(var = 'pathway')
    return(tmp_data)
}
R_Post_R_Pre_fgsea_tidy <- tidy_data(fgsea_list = R_Post_R_Pre_fgsea)
NR_Post_NR_Pre_fgsea_tidy <- tidy_data(fgsea_list = NR_Post_NR_Pre_fgsea)
NR_Post_R_Post_fgsea_tidy <- tidy_data(fgsea_list = NR_Post_R_Post_fgsea)
NR_Pre_R_Pre_fgsea_tidy <- tidy_data(fgsea_list = NR_Pre_R_Pre_fgsea)
```

## 04 绘制热图

```R
plot_heatmap <- function(mat_wider, column_title = column_title){
    ht <- Heatmap(
        as.matrix(mat_wider),
        col = col_fun,
        cluster_rows = F,
        cluster_columns = F,
        # row_names_gp = gpar(col),
        row_split = rep(1:nrow(mat_wider),each = 1),
        column_split = rep(1:ncol(mat_wider),each = 1),
        show_column_names = T,
        column_names_rot = 30,
        show_row_names = T,
        border=F,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 1:ncol(mat_wider)))),
        column_names_side = "top",
        column_title_gp = gpar(fonsize = 20, y = 0.5),
        row_names_side = "left",
        name = "NES (GSEA Preranked,p <0.05)",
        row_title = NULL,
        column_title = column_title,
        row_names_gp = gpar(fontsize = 8,col = 'black'),
        width = 0.6,
        height = 1,
        # cell_fun = function(j, i, x, y, width, height, fill) {
            # grid.text(sprintf("%.1f", mat_wider[i, j]), x, y, gp = gpar(fontsize = 10))
    # },
        na_col = '#d3d2d2',
        heatmap_legend_param = list(legend_direction = 'horizontal',legend_position ="left") # 将legend放到图像的横放
    )
    return(ht)
}

# 分群绘制
R_Post_R_Pre_fgsea_tidy_ht <- plot_heatmap(mat_wider = R_Post_R_Pre_fgsea_tidy, column_title = 'R_Post_R_Pre')
NR_Post_NR_Pre_fgsea_tidy_ht <- plot_heatmap(mat_wider = NR_Post_NR_Pre_fgsea_tidy, column_title = 'NR_Post_NR_Pre')
NR_Post_R_Post_fgsea_tidy_ht <- plot_heatmap(mat_wider = NR_Post_R_Post_fgsea_tidy, column_title = 'NR_Post_R_Post')
NR_Pre_R_Pre_fgsea_tidy_ht <- plot_heatmap(mat_wider = NR_Pre_R_Pre_fgsea_tidy, column_title = 'NR_Pre_R_Pre')
```

## 05 拼接热图

```R
# 横向拼接
h_all <- R_Post_R_Pre_fgsea_tidy_ht + NR_Post_NR_Pre_fgsea_tidy_ht + NR_Post_R_Post_fgsea_tidy_ht + NR_Pre_R_Pre_fgsea_tidy_ht
h_all <- draw(h_all, heatmap_legend_side = "bottom")

png('/root/wangje/Project/刘老师/Myeloids/Fig/GSEA/GSEA如图/髓系大群_热图.png', height = 2000,width = 3000,res=300)
h_all
dev.off()
```

[![ppRDNOx.md.png](https://s1.ax1x.com/2023/03/31/ppRDNOx.md.png)](https://imgse.com/i/ppRDNOx)

