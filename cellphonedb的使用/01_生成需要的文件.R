library(Seurat)
library(qs)
library(patchwork)
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyfast)   

#<<<<<<<< 读入NK& T cells文件 >>>>>>>>>>
setwd('/root/wangje/Project/刘老师/NK_T/Data')
load("new_NKT_Seurat.RData")
# 修改celltype
scRNA_seurat$cellphonedb_celltype = ifelse(
    scRNA_seurat$new_celltype %in% c('Treg','SOX4+ CD4','TH1-like','TFH cells', 'Naive T cells') ,'CD4+ T cells',ifelse(
        scRNA_seurat$new_celltype %in% c('Proliferating cells'), 'Proliferating cells', ifelse(
            scRNA_seurat$new_celltype %in% c('NK cells'), 'NK cells', 'CD8+ T cells'
        )
    )
)
# 去除Proliferating cells
scRNA_seurat = scRNA_seurat[,scRNA_seurat$cellphonedb_celltype != 'Proliferating cells']
# 生成矩阵文件
counts.nkt=as.data.frame(scRNA_seurat@assays$RNA@counts)
# 生成meta文件
Idents(scRNA_seurat) = scRNA_seurat$cellphonedb_celltype
meta.nkt=as.data.frame(scRNA_seurat@active.ident)
meta.nkt$cell=row.names(meta.nkt)
colnames(meta.nkt) = c('cellphonedb_celltype','cell')
meta.nkt=meta.nkt[,c("cell","cellphonedb_celltype")]

