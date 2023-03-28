library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(gghalves)
library(ggunchained)
library(tidyverse)

# 读入数据
load('/root/wangje/Project/刘老师/Myeloids/Data/Myeloids_CCA.RData')




# ------------------------封装成函数-------------------------------
calculate_IDO1 <- function(gene, scRNA_seurat) {
  IDO1 <- data.frame(matrix(1:24, ncol = 2))
  rownames(IDO1) <- unique(scRNA_seurat$celltype)
  colnames(IDO1) <- c("R", "NR")
  # 循环计算
  for (i in unique(scRNA_seurat$celltype)) {
    tmp_R_Pre <- scRNA_seurat[, scRNA_seurat$celltype == i & scRNA_seurat$Treat_assess == "R_Pre"] %>%
      FetchData(vars = gene) %>%
      rename_with(~"value", 1)
    tmp_R_Post <- scRNA_seurat[, scRNA_seurat$celltype == i & scRNA_seurat$Treat_assess == "R_Post"] %>%
      FetchData(vars = gene) %>%
      rename_with(~"value", 1)
    tmp_NR_Pre <- scRNA_seurat[, scRNA_seurat$celltype == i & scRNA_seurat$Treat_assess == "NR_Pre"] %>%
      FetchData(vars = gene) %>%
      rename_with(~"value", 1)
    tmp_NR_Post <- scRNA_seurat[, scRNA_seurat$celltype == i & scRNA_seurat$Treat_assess == "NR_Post"] %>%
      FetchData(vars = gene) %>%
      rename_with(~"value", 1)

    R_Pre_Post <- wilcox.test(tmp_R_Post$value, tmp_R_Pre$value)$p.value
    NR_Pre_Post <- wilcox.test(tmp_NR_Post$value, tmp_NR_Pre$value)$p.value

    IDO1[rownames(IDO1) == i, colnames(IDO1) == "R"] <- (-log10(R_Pre_Post))
    print(paste0(i, "---->",(- log10(R_Pre_Post))))
    IDO1[rownames(IDO1) == i, colnames(IDO1) == "NR"] <- (-log10(NR_Pre_Post))
    print(paste0(i, "---->",(- log10(NR_Pre_Post))))
  }
  return(IDO1)
}

IDO1 <- calculate_IDO1(gene = "IDO1", scRNA_seurat)
DefaultAssay(scRNA_seurat) <- 'RNA'
PDCD1LG2 <- calculate_IDO1(gene = 'PDCD1LG2', scRNA_seurat)
CD274 <- calculate_IDO1(gene = 'CD274',scRNA_seurat)

###############################################################
######### 使用FindMarkers函数进行计算
Idents(scRNA_seurat) <- 'celltype'
DefaultAssay(scRNA_seurat) <- 'RNA'
FindMarkers(scRNA_seurat,
            ident.1 = 'R_Post',
            ident.2 = 'R_Pre',
            group.by = 'Treat_assess',
            only.pos = FALSE)
flist <- list()
for(i in unique(scRNA_seurat$celltype)){
  print(i)
  flist[[i]] <- FindMarkers(scRNA_seurat,
            ident.1 = 'R_Post',
            ident.2 = 'R_Pre',
            group.by = 'Treat_assess',
            subset.ident = i,
            only.pos = FALSE)
  print(dim(flist[[i]]))
}

flist2 <- list()
Idents(scRNA_seurat) <- 'celltype'
for(i in unique(scRNA_seurat$celltype)){
  print(i)
  flist2[[i]] <- FindMarkers(scRNA_seurat,
            ident.1 = 'NR_Post',
            ident.2 = 'NR_Pre',
            group.by = 'Treat_assess',
            subset.ident = i,
            only.pos = FALSE)
  print(dim(flist2[[i]]))
}

flist3 <- list()
Idents(scRNA_seurat) <- 'celltype'
for(i in unique(scRNA_seurat$celltype)){
  print(i)
  flist3[[i]] <- FindMarkers(scRNA_seurat,
            ident.1 = 'NR_Post',
            ident.2 = 'R_Post',
            group.by = 'Treat_assess',
            subset.ident = i,
            only.pos = FALSE)
  print(dim(flist3[[i]]))
}

flist4 <- list()
Idents(scRNA_seurat) <- 'celltype'
for(i in unique(scRNA_seurat$celltype)){
  print(i)
  flist4[[i]] <- FindMarkers(scRNA_seurat,
            ident.1 = 'NR_Pre',
            ident.2 = 'R_Pre',
            group.by = 'Treat_assess',
            subset.ident = i,
            only.pos = FALSE)
  print(dim(flist4[[i]]))
}


plist1 <- list()
for (i in names(flist4)) {
  print(i)
  p <- EnhancedVolcano::EnhancedVolcano(
    flist4[[i]],
    lab = rownames(flist4[[i]]),
    x = "avg_log2FC",
    y = "p_val_adj",
    xlab = bquote(~ log[2] ~ "fold change"),
    ylab = bquote(~ -log[10] ~ "adjusted p-value"),
    cutoffLineType = "longdash",
    cutoffLineCol = "red",
    pCutoff = 0.05,
    FCcutoff = 0.2,
    axisLabSize = 25,
    title = paste0(i, "_NR_Pre vs R_Pre"),
    subtitle = NULL,
    titleLabSize = 25,
    labSize = 5,
    gridlines.major = F,
    gridlines.minor = F,
    border = "full",
    borderWidth = 1,
    pointSize = 3
  ) + theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25)
  )
  plist1[[i]] <- p
}

# dir.create('/root/wangje/Project/刘老师/Myeloids/Fig/火山图')
png('/root/wangje/Project/刘老师/Myeloids/Fig/火山图/髓系NR_Pre_vs_R_Pre.png', height = 3000,width = 4000)
wrap_plots(plist1, ncol = 4)
dev.off()

################################################################
######## 绘制热图
col_fun = colorRamp2(c(0, 1.3, 8), c("blue", "white", "red"))
#-----> CD274
h1 <- Heatmap(as.matrix(t(CD274)), 
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        column_names_gp = gpar(fontsize = 20),
        column_names_rot = 45,
        column_split = rep(LETTERS[1:12],1),
        row_split = rep(LETTERS[1:2]),
        row_names_gp = gpar(col = 'black',fontsize=20, fill='grey'),
        row_title = 'CD274(PD-L1)',
        row_title_gp = gpar(col = 'black', fontface = 'bold', fontsize = 12),
        row_title_rot = 0,
        column_title = NULL,
        width = unit(16, "cm"), height = unit(1.5, "cm"),
        name = '-log10(p-value)',
        col = col_fun,
        border_gp = gpar(col = '#cbcbcb')
) 

# -----> IDO1
h2 <- Heatmap(as.matrix(t(IDO1)), 
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        column_names_gp = gpar(fontsize = 20),
        column_names_rot = 45,
        column_split = rep(LETTERS[1:12],1),
        row_split = rep(LETTERS[1:2]),
        row_names_gp = gpar(col = 'black',fontsize=20, fill='grey'),
        row_title = 'IDO1',
        row_title_gp = gpar(col = 'black', fontface = 'bold', fontsize = 12),
        row_title_rot = 0,
        column_title = NULL,
        width = unit(16, "cm"), height = unit(1.5, "cm"),
        name = '-log10(p-value)',
        col = col_fun,
        border_gp = gpar(col = '#cbcbcb')
)

# ------> PD-L2
h3 <- Heatmap(as.matrix(t(PDCD1LG2)), 
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        column_names_gp = gpar(fontsize = 20),
        column_names_rot = 45,
        column_split = rep(LETTERS[1:12],1),
        row_split = rep(LETTERS[1:2]),
        row_names_gp = gpar(col = 'black',fontsize=20, fill='grey'),
        row_title = 'PDCD1LG2(PD-L2)',
        row_title_gp = gpar(col = 'black', fontface = 'bold', fontsize = 12),
        row_title_rot = 0,
        column_title = NULL,
        width = unit(16, "cm"), height = unit(1.5, "cm"),
        name = '-log10(p-value)',
        col = col_fun,
        border_gp = gpar(col = '#cbcbcb'),
        na_col = 'black'
)
png('/root/wangje/Project/刘老师/Myeloids/complexheatmap/PD-L1_PD-L2_IDO--Post-Pre.png',height = 2000,width = 4000,res=300)
draw(h1 %v% h2 %v% h3)
dev.off()

head(flist[[1]])
for(i in names(flist)){
  print(i)
  print(flist[[i]][rownames(flist[[i]]) == 'IDO1',])
}

#############################################################################
##################### 绘制髓系大群火山图
library(Seurat)
library(patchwork)
library(dplyr)
library(tidyverse)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(EnhancedVolcano)
# 读如数据
load("/root/wangje/Project/刘老师/Myeloids/Data/Myeloids_CCA.RData")
scRNA_seurat$new_celltype <- case_when(
  scRNA_seurat$celltype %in% grep('^Macro',unique(scRNA_seurat$celltype), value = T) ~ 'Macrophage',
  scRNA_seurat$celltype %in% c('Mono_CD14') ~ 'Mono_CD14',
  scRNA_seurat$celltype %in% c('migDC','cDC1_CLEC9A','cDC2_CD1C') ~ 'DC',
  scRNA_seurat$celltype %in% c('Mast cells') ~'Mast cells',
  TRUE ~ 'Proliferating'
)
# save(scRNA_seurat, file = '/root/wangje/Project/刘老师/Myeloids/Data/Myeloids_CCA.RData')
# flist <- list()
# for(i in unique(scRNA_seurat$new_celltype)){
#   DefaultAssay(scRNA_seurat) <- 'RNA'
#   Idents(scRNA_seurat) <- 'new_celltype'
#   print(i)
#   flist[[i]] <- Seurat::FindMarkers(
#     object = scRNA_seurat,
#     ident.1 = 'NR_Pre',
#     ident.2 = 'NR_Post',
#     subset.indent = i,
#     group.by = 'Treat_assess',
#     only.pos = FALSE
#   )
# }
# Myeloids_NR_Pre_NR_Post <- flist

find_markers_by_celltype <- function(scRNA_seurat, ident.1, ident.2){
  flist <- list()
  DefaultAssay(scRNA_seurat) <- 'RNA'
  Idents(scRNA_seurat) <- 'new_celltype'
  for(i in unique(scRNA_seurat$new_celltype)){
    print(i)
    flist[[i]] <- Seurat::FindMarkers(
      object = scRNA_seurat,
      ident.1 = ident.1,
      ident.2 = ident.2,
      group.by = 'Treat_assess',
      subet.ident = i,
      only.pos = FALSE
    )
    print(dim(flist[[i]]))
  }
  return(flist)
}
R_Post_R_Pre <- find_markers_by_celltype(scRNA_seurat, ident.1 = 'R_Post', ident.2 = 'R_Pre')
NR_Post_NR_Pre <- find_markers_by_celltype(scRNA_seurat, ident.1 = 'NR_Post', ident.2 = 'NR_Pre')

find_markers_by_celltype <- function(scRNA_seurat, ident.1, ident.2){
  flist <- list()
  DefaultAssay(scRNA_seurat) <- 'RNA'
  Idents(scRNA_seurat) <- 'celltype'
  for(i in unique(scRNA_seurat$celltype)){
    print(i)
    test <- Seurat::FindMarkers(
      object = scRNA_seurat,
      ident.1 = ident.1,
      ident.2 = ident.2,
      subset.indent = i,
      group.by = 'Treat_assess',
      only.pos = FALSE
    )
  }
  return(dim(test))
}

######################################################################
########### 寻找差异基因

DefaultAssay(scRNA_seurat) <- 'RNA'
Idents(scRNA_seurat) <- scRNA_seurat$new_celltype

Find_DEGs <- function(data,ident.1 = 'R_Post', ident.2 = 'R_Pre', group='Treat_assess'){
  flist = list()
  if(DefaultAssay(data) != 'RNA'){
    DefaultAssay(data) = 'RNA'
    Idents(data) = 'new_celltype'
  } 
  for(celltype in unique(data$new_celltype)){
    print(celltype)
    res = FindMarkers(
      object = data,
      ident.1 = ident.1,
      ident.2 = ident.2,
      group.by = group,
      subset.ident = celltype,
      only.pos = F
    )
    print(dim(res))
    flist[[celltype]] <- res
  }
  return(flist)
}

R_Post_R_Pre <- Find_DEGs(data = scRNA_seurat, group = 'Treat_assess',ident.1 = 'R_Post', ident.2 = 'R_Pre')
NR_Post_NR_Pre <- Find_DEGs(data = scRNA_seurat, ident.1 = 'NR_Post', ident.2 = 'NR_Pre',group = 'Treat_assess')
NR_Post_R_Post <- Find_DEGs(data = scRNA_seurat, ident.1 = 'NR_Post', ident.2 = 'R_Post',group = 'Treat_assess')
NR_Pre_R_Pre <- Find_DEGs(data = scRNA_seurat, ident.1 = 'NR_Pre', ident.2 = 'R_Pre',group = 'Treat_assess')

##############################################################################
######### 进行GO富集分析
library(GOSemSim)
library(clusterProfiler)
library(org.Hs.eg.db)

# data = R_Post_R_Pre
# plist_up = list()
# plist_down = list()
# for(i in names(data)){
# print(i)
# 区分出上调和下调的基因 avg_log2FC
# data_up = data[[i]] %>% dplyr::filter(avg_log2FC >0)
# data_down = data[[i]] %>% dplyr::filter(avg_log2FC <0)
# tmp_up <- enrichGO(
#                    gene = rownames(data_up),
#                    OrgDb = "org.Hs.eg.db",
#                    keyType = "SYMBOL",
#                    ont = "BP",
#                    pAdjustMethod = "BH"
# )
# tmp_down <- enrichGO(
#                    gene = rownames(data_down),
#                    OrgDb = "org.Hs.eg.db",
#                    keyType = "SYMBOL",
#                    ont = "BP",
#                    pAdjustMethod = "BH",
# )
# tmp_up <- enrichGO(
#                    gene = rownames(data_up),
#                    OrgDb = "org.Hs.eg.db",
#                    keyType = "SYMBOL",
#                    ont = "BP",
#                    pAdjustMethod = "BH",
# )
# print(dim(tmp_up))
# print(dim(tmp_down))
# plist_up[[i]] = barplot(tmp_up, showCategory = 15) + labs(title= paste0(i,'R_Post vs R_Pre up'))
# plist_down[[i]] = barplot(tmp_down, showCategory = 15) + labs(title= paste0(i,'R_Post vs R_Pre down'))
# 
# }


go_enrichment_analysis <- function(data,
                                   title1 = '',
                                   title2 = '') {
    plist_up <- list()
    plist_down <- list()
    
    for (i in names(data)) {
        print(i)
        # 区分出上调和下调的基因 avg_log2FC
        data_up = data[[i]] %>% dplyr::filter(avg_log2FC > 0)
        data_down = data[[i]] %>% dplyr::filter(avg_log2FC < 0)
        
        tmp_up <- enrichGO(
            gene = rownames(data_up),
            OrgDb = "org.Hs.eg.db",
            keyType = "SYMBOL",
            ont = "BP",
            pAdjustMethod = "BH"
        )
        
        tmp_down <- enrichGO(
            gene = rownames(data_down),
            OrgDb = "org.Hs.eg.db",
            keyType = "SYMBOL",
            ont = "BP",
            pAdjustMethod = "BH"
        )
        
        plist_up[[i]] = barplot(tmp_up, showCategory = 15) + labs(title = paste0(i, '-->', title1))
        plist_down[[i]] = barplot(tmp_down, showCategory = 15) + labs(title = paste0(i, '-->', title2))
    }
    
    result_list <- list(plist_up, plist_down)
    return(result_list)
}

R_Post_R_Pre_GO <- go_enrichment_analysis(data = R_Post_R_Pre, title1 = 'R_Post vs R_Pre up', title2 = 'R_Post vs R_Pre down')
NR_Post_NR_Pre_GO <- go_enrichment_analysis(data = NR_Post_NR_Pre, title1 = 'NR_Post vs NR_Pre up', title2 = 'NR_Post vs NR_Pre down')
NR_Post_R_Post_GO <- go_enrichment_analysis(data = NR_Post_R_Post, title1 = 'NR_Post vs R_Post up', title2 = 'NR_Post vs R_Post down')
NR_Pre_R_Pre_GO <- go_enrichment_analysis(data = NR_Pre_R_Pre, title1 = 'NR_Pre vs R_Pre up', title2 = 'NR_Pre vs R_Pre down')

# save plot
# R_Post vs R_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/GO富集/R_Post_R_Pre_GO.png', height = 6000, width = 12000, res = 300)
p1 <- wrap_plots(R_Post_R_Pre_GO[[1]], ncol = 5)
p2 <- wrap_plots(R_Post_R_Pre_GO[[2]], ncol = 5)
p1/p2
dev.off()

# NR_Post vs NR_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/GO富集/NR_Post_NR_Pre_GO.png', height = 6000, width = 12000, res = 300)
p1 <- wrap_plots(NR_Post_NR_Pre_GO[[1]], ncol = 5)
p2 <- wrap_plots(NR_Post_NR_Pre_GO[[2]], ncol = 5)
p1/p2
dev.off()

# NR_Post vs R_Post
png('/root/wangje/Project/刘老师/Myeloids/Fig/GO富集/NR_Post_R_Post_GO.png', height = 6000, width = 12000, res = 300)
p1 <- wrap_plots(NR_Post_R_Post_GO[[1]], ncol = 5)
p2 <- wrap_plots(NR_Post_R_Post_GO[[2]], ncol = 5)
p1/p2
dev.off()

# NR_Pre vs R_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/GO富集/NR_Pre_R_Pre_GO.png', height = 6000, width = 12000, res = 300)
p1 <- wrap_plots(NR_Pre_R_Pre_GO[[1]], ncol = 5)
p2 <- wrap_plots(NR_Pre_R_Pre_GO[[2]], ncol = 5)
p1/p2
dev.off()


####################
## 髓系小群数据
# save plot
# R_Post vs R_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/GO富集/celltype_R_Post_R_Pre_GO.png', height = 12000, width = 13000, res = 300)
p1 <- wrap_plots(R_Post_R_Pre_GO[[1]], ncol = 6)
p2 <- wrap_plots(R_Post_R_Pre_GO[[2]], ncol = 6)
p1/p2
dev.off()

# NR_Post vs NR_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/GO富集/celltype_NR_Post_NR_Pre_GO.png', height = 12000, width = 13000, res = 300)
p1 <- wrap_plots(NR_Post_NR_Pre_GO[[1]], ncol = 6)
p2 <- wrap_plots(NR_Post_NR_Pre_GO[[2]], ncol = 6)
p1/p2
dev.off()

# NR_Post vs R_Post
png('/root/wangje/Project/刘老师/Myeloids/Fig/GO富集/celltype_NR_Post_R_Post_GO.png', height = 12000, width = 12000, res = 300)
p1 <- wrap_plots(NR_Post_R_Post_GO[[1]], ncol = 6)
p2 <- wrap_plots(NR_Post_R_Post_GO[[2]], ncol = 6)
p1/p2
dev.off()

# NR_Pre vs R_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/GO富集/celltype_NR_Pre_R_Pre_GO.png', height = 12000, width = 13000, res = 300)
p1 <- wrap_plots(NR_Pre_R_Pre_GO[[1]], ncol = 6)
p2 <- wrap_plots(NR_Pre_R_Pre_GO[[2]], ncol = 6)
p1/p2
dev.off()

save.image(file = '/root/wangje/Project/刘老师/Myeloids/Fig/GO富集/髓系小群GO.RData')
################################################################################
######### KEGG 分析

kegg_enrichmernt_analysis <- function(data, title1 = '', title2 = ''){
    if(!typeof(data) == 'list'){
        cat('输入数据不是列表结构')
    }
    plist_up <- list()
    plist_down <- list()
    for (i in names(data)) {
        print(i)
        # 提取出上下调的基因
        tmp_up <- data[[i]] %>% dplyr::filter(avg_log2FC > 0 )
        tmp_down <- data[[i]] %>% dplyr::filter(avg_log2FC < 0 )
        # 讲gene symbol转换为Enterz id
        genelist_up <- clusterProfiler::bitr(
            rownames(tmp_up),
            fromType = 'SYMBOL',
            toType = 'ENTREZID',
            OrgDb = 'org.Hs.eg.db'
            ) %>% pull(ENTREZID)
        
        genelist_down <- clusterProfiler::bitr(
            rownames(tmp_down),
            fromType = 'SYMBOL',
            toType = 'ENTREZID',
            OrgDb = 'org.Hs.eg.db'
        ) %>% pull(ENTREZID)
        print(length(genelist_up))
        print(length(genelist_down))
        
        kegg_up <- clusterProfiler::enrichKEGG(
            gene = genelist_up,
            organism = 'hsa'
        )
        kegg_down <- clusterProfiler::enrichKEGG(
            gene = genelist_down,
            organism = 'hsa'
        )
        plist_up[[i]] <- barplot(kegg_up,showCategory = 15) + labs(title = paste0(i, '--->',title1))
        plist_down[[i]] <- barplot(kegg_down,showCategory = 15) + labs(title = paste0(i, '--->',title2))
    }
    result_list <- list(plist_up, plist_down)
    return(result_list)
}
R_Post_R_Pre_kegg <- kegg_enrichmernt_analysis(R_Post_R_Pre,title1 = 'R_Post vs R_Pre Up', title2 = 'R_Post vs R_Pre down')
NR_Post_NR_Pre_kegg <- kegg_enrichmernt_analysis(NR_Post_NR_Pre,title1 = 'NR_Post vs NR_Pre Up', title2 = 'NR_Post vs NR_Pre down')
NR_Post_R_Post_kegg <- kegg_enrichmernt_analysis(NR_Post_R_Post,title1 = 'NR_Post vs R_Post Up', title2 = 'NR_Post vs R_Post down')
NR_Pre_R_Pre_kegg <- kegg_enrichmernt_analysis(NR_Pre_R_Pre,title1 = 'NR_Pre vs R_Pre Up', title2 = 'NR_Pre vs R_Pre down')




# save plot
# R_Post vs R_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/R_Post_R_Pre_kegg.png', height = 6000, width = 12000, res = 300)
p1 <- wrap_plots(R_Post_R_Pre_kegg[[1]], ncol = 5)
p2 <- wrap_plots(R_Post_R_Pre_kegg[[2]], ncol = 5)
p1/p2
dev.off()

# NR_Post vs NR_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/NR_Post_NR_Pre_kegg.png', height = 6000, width = 12000, res = 300)
p1 <- wrap_plots(NR_Post_NR_Pre_kegg[[1]], ncol = 5)
p2 <- wrap_plots(NR_Post_NR_Pre_kegg[[2]], ncol = 5)
p1/p2
dev.off()

# NR_Post vs R_Post
png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/NR_Post_R_Post_kegg.png', height = 6000, width = 12000, res = 300)
p1 <- wrap_plots(NR_Post_R_Post_kegg[[1]], ncol = 5)
p2 <- wrap_plots(NR_Post_R_Post_kegg[[2]][c(1,3,4,5)], ncol = 4) 
p1/p2
dev.off()

NR_Post_R_Post_kegg[[2]][c(1,3,4,5)]

# NR_Pre vs R_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/NR_Pre_R_Pre_kegg.png', height = 6000, width = 10000, res = 300)
p1 <- wrap_plots(NR_Pre_R_Pre_kegg[[1]], ncol = 5)
p2 <- wrap_plots(NR_Pre_R_Pre_kegg[[1]], ncol = 5) 
p1/p2
dev.off()
save.image(file = "/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/KEGG和GO富集.RData")

## --------------------------------------------------------------------------------
## kegg celltype
## --------------------------------------------------------------------------------
# R_Post vs R_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/R_Post_R_Pre_celltype_kegg.png', height = 12000, width = 15000, res = 300)
p1 <- wrap_plots(R_Post_R_Pre_kegg[[1]], ncol = 6)
p2 <- wrap_plots(R_Post_R_Pre_kegg[[2]], ncol = 6)
p1/p2
dev.off()

# NR_Post vs NR_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/NR_Post_NR_Pre_celltype_kegg.png', height = 12000, width = 15000, res = 300)
p1 <- wrap_plots(NR_Post_NR_Pre_kegg[[1]], ncol = 6)
p2 <- wrap_plots(NR_Post_NR_Pre_kegg[[2]], ncol = 6)
p1/p2
dev.off()

# NR_Post vs R_Post
png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/NR_Post_R_Post_celltype_kegg.png', height = 12000, width = 12000, res = 300)
p1 <- wrap_plots(NR_Post_R_Post_kegg[[1]], ncol = 6)
p2 <- wrap_plots(NR_Post_R_Post_kegg[[2]][c(1,2,4:12)], ncol = 6) 
p1/p2
dev.off()


# NR_Pre vs R_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/NR_Pre_R_Pre_celltype_kegg.png', height = 12000, width = 15000, res = 300)
p1 <- wrap_plots(NR_Pre_R_Pre_kegg[[1]], ncol = 6)
p2 <- wrap_plots(NR_Pre_R_Pre_kegg[[2]], ncol = 6) 
p1/p2
dev.off()



###############################################################################
########### Reactome
# install.packages('ReactomePA')
# BiocManager::install('ReactomePA')
# dir.create('/root/wangje/Project/刘老师/Myeloids/Fig/Reactome富集')
library(ReactomePA)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(org.Hs.eg.db)

reactome_enrichmernt_analysis <- function(data, title1 = '', title2 = ''){
    if(!typeof(data) == 'list'){
        cat('输入数据不是列表结构')
    }
    plist_up <- list()
    plist_down <- list()
    for (i in names(data)) {
        print(i)
        # 提取出上下调的基因
        tmp_up <- data[[i]] %>% dplyr::filter(avg_log2FC > 0 )
        tmp_down <- data[[i]] %>% dplyr::filter(avg_log2FC < 0 )
        # 讲gene symbol转换为Enterz id
        genelist_up <- clusterProfiler::bitr(
            rownames(tmp_up),
            fromType = 'SYMBOL',
            toType = 'ENTREZID',
            OrgDb = 'org.Hs.eg.db'
        ) %>% pull(ENTREZID)
        
        genelist_down <- clusterProfiler::bitr(
            rownames(tmp_down),
            fromType = 'SYMBOL',
            toType = 'ENTREZID',
            OrgDb = 'org.Hs.eg.db'
        ) %>% pull(ENTREZID)
        
        #进行reactome富集分析
        reactome_up = enrichPathway(
            gene = genelist_up,
            organism = 'human',
            pAdjustMethod = 'BH'
        )
        
        reactome_down = enrichPathway(
            gene = genelist_down,
            organism = 'human',
            pAdjustMethod = 'BH'
        )
        
        # flist_up[[i]] <- reactome_up %>% as.data.frame()
        # flist_down[[i]] <- reactome_down %>% as.data.frame()
        
        
        plist_up[[i]] <- barplot(reactome_up, showCategory = 15) + labs(title = paste0(i, '--->',title1))
        plist_down[[i]] <- barplot(reactome_down, showCategory = 15) + labs(title = paste0(i, '--->',title2))
    }
    result_list <- list(plist_up, plist_down)
    return(result_list)
}

R_Post_R_Pre_Reactome <- reactome_enrichmernt_analysis(R_Post_R_Pre,title1 = 'R_Post vs R_Pre Up', title2 = 'R_Post vs R_Pre down')
NR_Post_NR_Pre_Reactome <- reactome_enrichmernt_analysis(NR_Post_NR_Pre,title1 = 'NR_Post vs NR_Pre Up', title2 = 'NR_Post vs NR_Pre down')
NR_Post_R_Post_Reactome <- reactome_enrichmernt_analysis(NR_Post_R_Post,title1 = 'NR_Post vs R_Post Up', title2 = 'NR_Post vs R_Post down')
NR_Pre_R_Pre_Reactome <- reactome_enrichmernt_analysis(NR_Pre_R_Pre,title1 = 'NR_Pre vs R_Pre Up', title2 = 'NR_Pre vs R_Pre down')


## save plot
# R_Post vs R_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/Reactome富集/R_Post_R_Pre_Reactome.png', height = 6000, width = 12000, res = 300)
p1 <- wrap_plots(R_Post_R_Pre_Reactome[[1]], ncol = 5)
p2 <- wrap_plots(R_Post_R_Pre_Reactome[[2]], ncol = 5)
p1/p2
dev.off()

# NR_Post vs NR_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/Reactome富集/NR_Post_NR_Pre_Reactome.png', height = 6000, width = 12000, res = 300)
p1 <- wrap_plots(NR_Post_NR_Pre_Reactome[[1]], ncol = 5)
p2 <- wrap_plots(NR_Post_NR_Pre_Reactome[[2]], ncol = 5)
p1/p2
dev.off()

# NR_Post vs R_Post
png('/root/wangje/Project/刘老师/Myeloids/Fig/Reactome富集/NR_Post_R_Post_Reactome_02.png', height = 6000, width = 12000, res = 300)
p1 <- wrap_plots(NR_Post_R_Post_Reactome[[1]], ncol = 5)
p2 <- wrap_plots(NR_Post_R_Post_Reactome[[2]][c(1,3,4,5)], ncol = 4) 
p1/p2
dev.off()


# NR_Pre vs R_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/Reactome富集/NR_Pre_R_Pre_Reactome.png', height = 6000, width = 10000, res = 300)
p1 <- wrap_plots(NR_Pre_R_Pre_Reactome[[1]], ncol = 5)
p2 <- wrap_plots(NR_Pre_R_Pre_Reactome[[2]], ncol = 5) 
p1/p2
dev.off()

save.image(file = '/root/wangje/Project/刘老师/Myeloids/Fig/Reactome富集/髓系大群Reactome.RData')


## -----------------------------------------------------------------------------
## 髓系celltype
## -----------------------------------------------------------------------------
# R_Post vs R_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/Reactome富集/R_Post_R_Pre_celltype_Reactome.png', height = 12000, width = 15000, res = 300)
p1 <- wrap_plots(R_Post_R_Pre_Reactome[[1]], ncol = 6)
p2 <- wrap_plots(R_Post_R_Pre_Reactome[[2]], ncol = 6)
p1/p2
dev.off()

# NR_Post vs NR_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/Reactome富集/NR_Post_NR_Pre_Reactome.png', height = 12000, width = 15000, res = 300)
p1 <- wrap_plots(NR_Post_NR_Pre_Reactome[[1]], ncol = 6)
p2 <- wrap_plots(NR_Post_NR_Pre_Reactome[[2]], ncol = 6)
p1/p2
dev.off()

# NR_Post vs R_Post
png('/root/wangje/Project/刘老师/Myeloids/Fig/Reactome富集/NR_Post_R_Post_celltype_Reactome.png', height = 12000, width = 15000, res = 300)
p1 <- wrap_plots(NR_Post_R_Post_Reactome[[1]], ncol = 6)
p2 <- wrap_plots(NR_Post_R_Post_Reactome[[2]][c(1:2,4:12)], ncol = 6) 
p1/p2
dev.off()


# NR_Pre vs R_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/Reactome富集/NR_Pre_R_Pre_celltype_Reactome.png', height = 12000, width = 15000, res = 300)
p1 <- wrap_plots(NR_Pre_R_Pre_Reactome[[1]], ncol = 6)
p2 <- wrap_plots(NR_Pre_R_Pre_Reactome[[2]], ncol = 6) 
p1/p2
dev.off()


################################################################################
########### 使用AddmoudleScore 计算score
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(stringr)


setwd('/root/wangje/Project/刘老师/Myeloids/VISION_score')
gmt_file <- list.files('/root/wangje/Project/刘老师/Myeloids/VISION_score', pattern = 'gmt.txt$')

Run_AddMoudleScore <- function(object,gmt_path='/root/wangje/Project/刘老师/Myeloids/VISION_score'){
    DefaultAssay(object) <- 'RNA'
    gmt_file <-
        list.files(gmt_path, pattern = 'gmt.txt', full.names = F)
    result <- list()
    
    for (file in gmt_file) {
        print(file)
        name = str_split_fixed(file, '\\.', n = 2)[, 1]
        tmp <- read.table(file, header = F)
        gene <- as.vector(t(tmp)[3:dim(tmp)[2]])
        genelist <- list()
        genelist[[name]] <- gene
        print(genelist)
        
        res <- AddModuleScore(object = object,
                              features = genelist,
                              name = names(genelist))
        result[[name]] <- res@meta.data %>% dplyr::select(last_col())
    }
    return(result)
    
}

test <- Run_AddMoudleScore(object = scRNA_seurat,gmt_path ='/root/wangje/Project/刘老师/Myeloids/VISION_score')


#--------> 接上面的结果进行绘图
library(gghalves)
library(ggunchained)
library(ggpubr)

id1 <- c("CJME","CMDI","HDYU","HXZH","LCLI","WYJU","WZLA","ZXME","ZJLI_0116")
id2 <- c("CZYI","FHXHBS1","FYYI","HEJI","LAWE","ZEYI","ZFXI","LIPE")
id3 <- c("CJME_0707","CMDI_0624","HDYU_0720","HXZH_0220","LCLI_0623","WYJU_0122","ZXME_0223","ZJLI_0312")

Plot_AddmoudleScore <- function(object){
    plist <- list()
    for (i in seq_along(object)) {
        print(names(object[i]))
        # 添加新的信息
        object[[i]]$Patient <- str_split_fixed(rownames(object[[i]]),'_[A|T|G|C].*', n = 2)[,1]
        object[[i]]$group <- case_when(
            object[[i]]$Patient %in% id1 ~'R_Pre',
            object[[i]]$Patient %in% id2 ~ 'NR_Pre',
            object[[i]]$Patient %in% id3 ~ 'R_Post',
            TRUE ~'NR_Post'
        )
        object[[i]]$group1 <- str_split_fixed(object[[i]]$group,'_',n = 2)[,1]
        object[[i]]$group2 <- str_split_fixed(object[[i]]$group,'_',n=2)[,2]
        object[[i]]$group1 <- factor(object[[i]]$group1, levels = c('R','NR'))
        object[[i]]$group2 <- factor(object[[i]]$group2, levels = c('Pre','Post'))
        
        p <- ggplot(object[[i]], aes_string(x = "group1", y = colnames(object[[i]])[1], fill="group2"))+
            geom_split_violin(trim = TRUE) + 
            geom_half_boxplot(
                data = object[[i]] %>% dplyr::filter(group2 == 'Pre'), aes_string(x = 'group1', y = colnames(object[[i]])[1]),
                width = 0.15, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef = 0, color = "black", linetype = 1, fill = "white"
            )+
            geom_half_boxplot(
                data = object[[i]] %>% dplyr::filter(group2 == 'Post'), aes_string(x = 'group1', y = colnames(object[[i]])[1]),
                width = 0.15, side = 'r',notch = FALSE, notchwidth = .4, outlier.shape = NA, coef = 0, color = "black", linetype = 1, fill = "white"
            )+
            labs(x = NULL, y = "Signature score", fill = "Group", title = colnames(object[[i]])[1]) +
            theme_classic() +
            theme(
                text = element_text(size = 20), 
                axis.title = element_text(size = 23, color = "black"),
                axis.text = element_text(size = 20, colour = "black"),
                plot.title = element_text(hjust = 0.5)
            ) +
            scale_fill_manual(values = c("#3b9aac", "#a74947", "#4470a1", "#8da24f", "#6d5a87")) +
            ggpubr::stat_compare_means(method = "wilcox.test")
        plist[[i]] <- p
    }
    return(plist)
}
Plot_AddmoudleScore(object = test)

###############################################################################
########## GSEA 分析
library(ggplot2)
library(GSEABase)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)

# 使用msigdbr生成需要的gmt文件
# msigdbr_species()
# human <- msigdbr(species = "Homo sapiens")
# human[1:4,1:4]

# 读入基因集文件
geneSet = read.table('/root/wangje/Reference/GSEA_gmt/Human/GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I.gmt', header = F)
geneSet <- as.data.frame(t(geneSet))
geneSet$term <- geneSet[1,1]
colnames(geneSet) <- c('term','gene')
geneSet <- geneSet[,c(2,1)]
geneSet <- geneSet[4:nrow(geneSet),1:2]
rownames(geneSet) <- 1:nrow(geneSet)

head(geneSet)
# term  gene
# 1 GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I  TAP1
# 2 GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I  NCF2
# 3 GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I  TAP2
# 4 GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I   MR1
# 5 GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I PSMC5
# 6 GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I PSME2

# 进行Gase分析
GSEA_analysis <- function(DEGs, geneSet=geneSet){
    flist = list()
    if(typeof(DEGs) == 'list'){
        for(name in names(DEGs)){
            print(name)
            geneList = DEGs[[name]]$avg_log2FC
            names(geneList)= rownames(DEGs[[name]])
            geneList=sort(geneList,decreasing = T)
            print(head(geneList))
            egmt <- GSEA(geneList, TERM2GENE=geneSet, 
                         minGSSize = 1,
                         pvalueCutoff = 0.99,
                         verbose=FALSE)
            gsea_results_df <- egmt@result
            outFile = paste0('/root/wangje/Project/刘老师/Myeloids/Data/GSEA/',name,'_',unique(geneSet$term),'txt')
            write.table(gsea_results_df,file = outFile, sep = '\t',quote = F)
            flist[[name]] <- gseaplot2(egmt,geneSetID = unique(geneSet$term),pvalue_table=T)
            
        }
        return(flist)
    }
    
}

R_Post_R_Pre_GSEA <- GSEA_analysis(DEGs = R_Post_R_Pre, geneSet = geneSet)



#### 使用fgsea进行分析
ifng <- read.table('/root/wangje/Reference/GSEA_gmt/Human/IFNG.txt', header = T , sep = '\t')

fgsea_set = rbind(geneSet, ifng)
fgsea_analysis <- function(DEGs, fgsea_set = fgsea_set){
    flist <- list()
    if(typeof(DEGs) == 'list'){
        for(name in names(DEGs)){
            print(name)
            DEGs[[name]]$gene <- rownames(DEGs[[name]])
            geneList = DEGs[[name]] %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene, avg_log2FC)
            ranks<- deframe(geneList)
            test <- fgsea::fgsea(fgsea_sets, stats = ranks,nperm = 1000,nproc=1)
            flist[[name]] <- plotEnrichment(test[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
                                            ranks) + labs(title="R_Post_R_Pre_HALLMARK_INTERFERON_GAMMA_RESPONSE")
        }
        return(flist)
    }
}

fgsea_set <- list()
fgsea_set[['GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I']] <- 
    as.vector(geneSet$V1[4:94])

test <- fgsea_analysis(DEGs = R_Post_R_Pre, fgsea_set = fgsea_set)

plotEnrichment(test[["GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I"]],
               ranks) + labs(title="GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I")

# ---------------------------------------------------------------------
name = 1
DEGs <- NR_Post_NR_Pre
DEGs[[name]]$gene <- rownames(DEGs[[name]])
geneList = DEGs[[name]] %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene, avg_log2FC)
ranks<- deframe(geneList)
test <- fgsea::fgsea(fgsea_set, stats = ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000,nproc=1)

fgsea_set



NR_Post_R_Post[[1]]$gene <- rownames(NR_Post_R_Post[[1]])
geneList = NR_Post_R_Post[[1]] %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene, avg_log2FC)
ranks<- deframe(geneList)


test1 <- fgsea(fgsea_sets,stats = ranks, nperm=500,nproc=1)

plotEnrichment(test1[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
               ranks) + labs(title="HALLMARK_INTERFERON_GAMMA_RESPONSE")















################################################################################
########### 绘制PD-L1, PD_L2以及IDO基因的表达分数
DefaultAssay(scRNA_seurat) <- 'RNA'
violins_data <- VlnPlot(
    scRNA_seurat,
    group.by = 'celltype',
    # cols = c("limegreen", "navy"),
    features = c('CD274', 'IDO1', 'PDCD1LG2'),
    ncol = 3,
    log = FALSE,
    combine = FALSE
)

id1 <- c("CJME","CMDI","HDYU","HXZH","LCLI","WYJU","WZLA","ZXME","ZJLI_0116")
id2 <- c("CZYI","FHXHBS1","FYYI","HEJI","LAWE","ZEYI","ZFXI","LIPE")
id3 <- c("CJME_0707","CMDI_0624","HDYU_0720","HXZH_0220","LCLI_0623","WYJU_0122","ZXME_0223","ZJLI_0312")

plot_violins <- function(data, id1, id2, id3) {
    
    if(!is.list(data)){
        stop('输入的文件不是列表类型.')
    }
    
    plist <- list()
    
    for (i in seq_along(data)) {
        tmp_data <- data[[i]]$data
        tmp_data$sample <- str_extract(rownames(tmp_data), ".*(?=_[A|T|G|C])")
        tmp_data <- as.data.frame(tmp_data)
        tmp_data$group <- case_when(
            tmp_data$sample %in% id1 ~ "R_Pre",
            tmp_data$sample %in% id2 ~ "NR_Pre",
            tmp_data$sample %in% id3 ~ "R_Post",
            TRUE ~ "NR_Post"
        )
        
        unique_ident <- unique(tmp_data$ident)
        
        for (j in unique_ident) {
            tmp <- subset(tmp_data, ident == j)
            tmp$group <- factor(tmp$group, levels = c('R_Pre','R_Post','NR_Pre','NR_Post'))
            p <- ggplot(data = tmp, aes_string(x = 'group', y = colnames(tmp)[1], fill = 'group')) +
                geom_violin(trim = FALSE, width = 0.5) +
                geom_jitter(size = 0.05, alpha = 0.3, width = 0.3) +
                theme_classic() +
                labs(x = '', title = j) +
                theme(
                    axis.title = element_text(size = 20, face = 'bold', colour = 'black'),
                    axis.text.x = element_text(size = 20, angle = 45, colour = 'black', hjust = 1, vjust = 1),
                    axis.text.y = element_text(size = 20, colour = 'black'),
                    plot.title = element_text(size = 20, face = 'bold', color = 'black', hjust = 0.5)
                ) +
                stat_compare_means(
                    comparisons = list(c("R_Pre", "R_Post"), c("R_Pre", "NR_Pre"), c("R_Post", "NR_Pre"), c("NR_Pre", "NR_Post")),
                    method = 'wilcox.test'
                )
            
            plist[[paste0(i,'--',j)]] <- p
        }
    }
    
    return(plist)
}

plist2 <- plot_violins(data = violins_data, id1 = id1, id2 = id2, id3 = id3)
# dir.create('/root/wangje/Project/刘老师/Myeloids/Fig/PD-L1_PD-L2_IDO1')
png('/root/wangje/Project/刘老师/Myeloids/Fig/PD-L1_PD-L2_IDO1/PD-L1_PD-L2_IDO1_小提琴图.png', height = 6000 , width = 14000, res=300)
wrap_plots(plist2, ncol = 12)
dev.off()

##############
# 髓系大群
DefaultAssay(scRNA_seurat) <- 'RNA'
violins_data <- VlnPlot(
    scRNA_seurat,
    group.by = 'new_celltype',
    # cols = c("limegreen", "navy"),
    features = c('CD274', 'IDO1', 'PDCD1LG2'),
    ncol = 3,
    log = FALSE,
    combine = FALSE
)

plist1 <- plot_violins(data = violins_data, id1 = id1, id2 = id2, id3 = id3)


png('/root/wangje/Project/刘老师/Myeloids/Fig/PD-L1_PD-L2_IDO1/PD-L1_PD-L2_IDO1_大群_小提琴图.png', height = 6000 , width = 8000, res=300)
wrap_plots(plist1,ncol = 5)
dev.off()


DefaultAssay(scRNA_seurat)

VlnPlot(scRNA_seurat, features = 'IDO1', group.by = 'Treat_assess', split.by = 'celltype')

scRNA_seurat$Treat_assess <- factor(scRNA_seurat$Treat_assess, levels = c('R_Pre','R_Post','NR_Pre','NR_Post'))

plist <- list()
for (i in unique(scRNA_seurat$celltype)) {
    print(i)
    plist[[i]] <- VlnPlot(scRNA_seurat[,scRNA_seurat$celltype == i],features = c('PDCD1LG2'),
                          group.by = 'Treat_assess')+
        stat_compare_means(
            comparisons = list(c("R_Pre", "R_Post"), c("R_Pre", "NR_Pre"), c("R_Post", "NR_Pre"), c("NR_Pre", "NR_Post")),
            method = 'wilcox.test'
        )
}

p <- wrap_plots(plist, ncol = 12)
p1 <- wrap_plots(plist, ncol = 12)
p2 <- wrap_plots(plist, ncol = 12)

png('/root/wangje/Project/刘老师/Myeloids/Fig/PD-L1_PD-L2_IDO1/PD-L1_PD-L2_IDO1_小提琴图02.png', height = 6000 , width = 14000, res=300)
p
dev.off()














##############################################################################
########### 差异基因火山图
compaired <- data.frame(
    C1 = c('R_Post','NR_Pre','NR_Post', 'NR_Post'),
    C2 = c('R_Pre','R_Pre','NR_Pre','R_Post')
)
markers_results_P <- map2(compaired$C1, compaired$C2, function(x1, y1) {
    print(paste0(x1,"--",y1))
    DefaultAssay(scRNA_seurat) <- 'RNA'
    Seurat::FindMarkers(object = scRNA_seurat,
                      ident.1 = x1,
                      ident.2 = y1,
                      group.by = "Treat_assess",
                      subset.ident = 'Macrophage',
                      only.pos = FALSE)})

names(markers_results_P) <- paste(compaired$C1,compaired$C2,sep = "-vs-")


########################
plist1 <- list()
for (i in names(flist)) {
  print(i)
  p <- EnhancedVolcano::EnhancedVolcano(
    flist[[i]],
    lab = rownames(flist[[i]]),
    x = "avg_log2FC",
    y = "p_val_adj",
    xlab = bquote(~ log[2] ~ "fold change"),
    ylab = bquote(~ -log[10] ~ "adjusted p-value"),
    cutoffLineType = "longdash",
    cutoffLineCol = "red",
    pCutoff = 0.05,
    FCcutoff = 0.2,
    axisLabSize = 25,
    title = paste0(i, "_R_Pre vs NR_Pre"),
    subtitle = paste0('Up in R_Post <------   ------> Up in NR_Post'),
    titleLabSize = 25,
    labSize = 5,
    gridlines.major = F,
    gridlines.minor = F,
    border = "full",
    borderWidth = 1,
    pointSize = 3
  ) + theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25),
    plot.title = element_text(size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 15, hjust = 0.5, vjust = -1)
  )
  plist1[[i]] <- p
}

p1 <- wrap_plots(plist1, ncol = 5)
p2 <- wrap_plots(plist1, ncol = 5)
p3 <- wrap_plots(plist1, ncol = 5)
p4 <- wrap_plots(plist1, ncol = 5) 

png('/root/wangje/Project/刘老师/Myeloids/Fig/火山图/髓系大群火山图.png',height = 10000, width = 10000, res= 300)
p1/p2/p3/p4
dev.off()




# 使用ggplot2绘制火山图
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggrepel)
Plot_Volcanos <- function(data, log2FC = 0.1, q_val = 0.05, info1 = "R_Pre", info2 = "NR_Post", title = '', cols=c("#497aa2", "#ae3137")) {
  data <- data %>%
    mutate(threshold = dplyr::case_when(
      avg_log2FC > log2FC & p_val_adj < q_val ~ paste0("Up in ", info1),
      avg_log2FC < (-log2FC) & p_val_adj < q_val ~ paste0("Up in ", info2),
      TRUE ~ "NS"
    ))
  tb <- data %>% group_by(threshold) %>% count() %>% dplyr::rename_with(~c('group','n'),1:2)
  g1 <- ggplot(data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = threshold)) +
    geom_point(alpha = 0.8, size = 0.8) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype = 2, color = "grey") +
    geom_hline(yintercept = -log10(q_val), linetype = 2, color = "grey") +
    labs(title = title) +
    xlab(bquote(Log[2] * FoldChange)) +
    ylab(bquote(-Log[10] * italic(P.adj))) +
    theme(
      legend.box = "horizontal",
      legend.position = "top",
      legend.spacing.x = unit(0, "pt"),
      legend.text = element_text(margin = margin(r = 20)),
      legend.margin = margin(b = -10, unit = "pt"),
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    scale_color_manual(
      values = c(cols[1], cols[2], 'grey')
    ) + 
    ggrepel::geom_text_repel(
      data = data %>% dplyr::filter(threshold != 'NS'),
      aes(x = avg_log2FC, y = -log10(p_val_adj), label = rownames(data)),
      color = 'balck', alpha = 1
    )

  return(g1)
}

Plot_Volcanos(flist[[1]], title = 'test')
png('/root/wangje/Project/刘老师/Myeloids/Fig/火山图/test1.png', height = 2000,width = 2000,res=300)
Plot_Volcanos(flist[[1]], title = 'test')
dev.off()

test  <- flist[[1]] %>%
    as.tibble() %>%
    mutate(threshold = dplyr::case_when(
      avg_log2FC > log2FC & p_val_adj < q_val ~ paste0("Up in ", info1),
      avg_log2FC < (-log2FC) & p_val_adj < q_val ~ paste0("Up in ", info2),
      TRUE ~ "NS"
    ))

head(test)
test %>% group_by(threshold) %>% count()



#######################################################################
###################### 绘制大群热图
DefaultAssay(scRNA_seurat) <- 'RNA'
Idents(scRNA_seurat) <- scRNA_seurat$celltype
scRNA_seurat <- scRNA_seurat %>% ScaleData()
scRNA_seurat.Findmarkers <- FindAllMarkers(scRNA_seurat, only.pos = TRUE,BPPARAM = MulticoreParam(35))
Top10.genes <- scRNA_seurat.Findmarkers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
png('RNA.heatmap.png',height = 4000, width = 4000, res=300)
p4 <- DoHeatmap(subset(scRNA_seurat, downsample = 500), features = Top10.genes$gene, size = 4) + NoLegend() + theme(text = element_text(size = 6))
p4
dev.off()
