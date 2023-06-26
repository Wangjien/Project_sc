library(infercnv)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(AnnoProbe)
library()

# load data 
big.fastmnn = '/root/wangje/wangjien_new/liu-big/scRNA_FastMNN_new.RData'
load(big.fastmnn)

# 提取部分文件
# scRNA_sub = scRNA[,scRNA$celltype %in% c('Fibroblasts','Endothelials','Epithelials')]
# 使用全部的上皮数据进行分析，结果发现，内存会报错，所有将所有的上皮细胞分成小部分进行分析

step01 <- function(data = scRNA, celltype = celltype) {
    options(scipen = 100)
    if (!require("pacman")) {
        message("安装pacman")
        install.packages("pacman")
    } else {
        pacman::p_load("infercnv", "Seurat", "dplyr", "stringr", "AnnoProbe")
    } 
    # 获取表达矩阵
    if ("celltype" %in% colnames(data@meta.data)) {
        scRNA_sub <<- data[, data$celltype %in% celltype]
        scRNA_matrix <<- GetAssayData(scRNA_sub)
    } else {
        message("缺少celltype信息")
    }
    # 获得注释信息
    cellinfo <- as.data.frame(scRNA_sub$celltype)
    # 获得基因位置信息 <----- 注意这是使用的是人类信息 ----->
    geneInfor <- AnnoProbe::annoGene(rownames(scRNA_matrix), "SYMBOL", "human")
    geneInfor <- geneInfor[with(geneInfor, order(chr, start)), c(1, 4:6)]
    geneInfor <- geneInfor[!duplicated(geneInfor[, 1]), ]
    print("01----------------")
    # 去除性染色体
    new_geneInfor <- geneInfor[!(geneInfor$chr %in% c("chrM", "chrX", "chrY")), ]
    # 按照新的排列方式
    new_geneInfor$chr <- factor(new_geneInfor$chr, levels = paste0("chr", 1:22))
    # 表达矩阵按照基因文件和分组信息进行排序
    print("02-----------------")
    scRNA_matrix <- scRNA_matrix[rownames(scRNA_matrix) %in% new_geneInfor[, 1], ]
    scRNA_matrix <- scRNA_matrix[match(new_geneInfor[, 1], rownames(scRNA_matrix)), ]
    scRNA_matrix <- scRNA_matrix[, colnames(scRNA_matrix) %in% rownames(cellinfo)]
    scRNA_matrix <- scRNA_matrix[, match(rownames(cellinfo), colnames(scRNA_matrix))]
    print("03-----------------")
    # 输入部分信息
    print(sprintf("矩阵信息与基因信息相同%s", identical(rownames(scRNA_matrix), new_geneInfor$SYMBOL)))
    print(sprintf("矩阵信息与细胞信息相同%s", identical(colnames(scRNA_matrix), rownames(cellinfo))))
    save(scRNA_matrix, file = "./scRNA_matrix.RData")
    groupFiles <<- "groupFiles.txt" # nolint
    write.table(cellinfo, file = groupFiles, sep = "\t", col.names = F, row.names = T, quote = F)
    geneFile <<- "geneFile.txt"
    write.table(new_geneInfor, file = geneFile, sep = "\t", quote = F, col.names = F, row.names = F)
}
step01(scRNA,celltype = c('Fibroblasts','Endothelials','Epithelials'))

# 查看celltype数目
cell = read.table('/root/wangje/Project/刘老师/new_infer_res/groupFiles.txt')
gene = read.table('/root/wangje/Project/刘老师/new_infer_res/geneFile.txt')
load('/root/wangje/Project/刘老师/new_infer_res/scRNA_matrix.RData')

#----------------------------------------------------------------
# 在R中的infercnv，细胞数量过对多会报错，可以分批次执行
#----------------------------------------------------------------
dim(scRNA_matrix) # 34112 160958
table(cell$V2) 
# Endothelials Epithelials Fibroblasts
# 21441           111181      28336
## 从中挑选出部分细胞进行分析
## 从内皮细胞中挑选出1万个细胞，从成纤维细胞中挑选出1万个稀薄，上皮细胞分为6个批次分别运行
## 内皮细胞选出1万个细胞
cell_endo = cell %>% filter(cell$V2 == 'Endothelials')
cell_endo = cell_endo[sample(cell_endo, 10000, replace=F),]

## 成纤维细胞挑选出1万个细胞
cell_fib = cell %>% filter(cell$V2 == 'Fibroblasts')
cell_fib = cell_fib[sample(cell_fib$V2 , 10000, replace=F),]

## 上皮细胞
cell_epi = c


