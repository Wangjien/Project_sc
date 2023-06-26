library(infercnv)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(AnnoProbe)
library(qs)

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
cell_endo = cell_endo[cell_endo$V1 %in% sample(cell_endo$V1, 10000, replace=F),]

## 成纤维细胞挑选出1万个细胞
cell_fib = cell %>% filter(cell$V2 == 'Fibroblasts')
cell_fib = cell_fib[cell_fib$V1 %in% sample(cell_fib$V1 , 10000, replace=F),]

## 上皮细胞
cell_epi = cell %>% filter(cell$V2 == 'Epithelials')
cell_epi_1 = cell_epi %>% slice(1:20000)

#---------------------------------------------------------------
# dataset 1
#---------------------------------------------------------------
dataset1 = rbind(cell_endo, cell_fib, cell_epi_1)
mtx1 = scRNA_matrix[,(which(colnames(scRNA_matrix) %in% dataset1$V1))]
write.table(dataset1,file = '/root/wangje/Project/刘老师/new_infer_res/dataset1/groupFiles.txt',sep = '\t', row.names=F, col.names=F, quote=F)
write.table(gene, file = '/root/wangje/Project/刘老师/new_infer_res/dataset1/geneFiles.txt',sep = '\t', row.names=F, col.names=F, quote=F)
qsave(mtx1, file = '/root/wangje/Project/刘老师/new_infer_res/dataset1/scRNA_matrix.qs')
# infercnv
scRNA_matrix = qread('/root/wangje/Project/刘老师/new_infer_res/dataset1/scRNA_matrix.qs')
geneFile = '/root/wangje/Project/刘老师/new_infer_res/dataset1/geneFiles.txt'
groupFiles = '/root/wangje/Project/刘老师/new_infer_res/dataset1/groupFiles.txt'
ref_group_names = c("Endothelials","Fibroblasts")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=scRNA_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=ref_group_names) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="/root/wangje/Project/刘老师/new_infer_res/dataset1/", 
                             cluster_by_groups=TRUE, # 聚类
                             denoise=TRUE,
                             HMM=FALSE,
                             # 运行的线程数目
                             num_threads = 40,
                             write_expr_matrix=TRUE
                             ) 

#----------------------------------------------------------------
# dataset2
#----------------------------------------------------------------
cell_epi_2 = cell_epi %>% slice(20001:40000)
dataset2 = rbind(cell_endo, cell_fib, cell_epi_2)
mtx2 = scRNA_matrix[,(which(colnames(scRNA_matrix) %in% dataset2$V1))]
qsave(mtx2, file = './scRNA_matrix.qs')
write.table(dataset2,file = '/root/wangje/Project/刘老师/new_infer_res/dataset2/groupFiles.txt',sep = '\t', row.names=F, col.names=F, quote=F)
write.table(gene, file = '/root/wangje/Project/刘老师/new_infer_res/dataset2/geneFiles.txt',sep = '\t', row.names=F, col.names=F, quote=F)

### 
scRNA_matrix = qread('/root/wangje/Project/刘老师/new_infer_res/dataset2/scRNA_matrix.qs')
geneFile = '/root/wangje/Project/刘老师/new_infer_res/dataset2/geneFiles.txt'
groupFiles = '/root/wangje/Project/刘老师/new_infer_res/dataset2/groupFiles.txt'
ref_group_names = c("Endothelials","Fibroblasts")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=scRNA_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=ref_group_names) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="/root/wangje/Project/刘老师/new_infer_res/dataset2/", 
                             cluster_by_groups=TRUE, # 聚类
                             denoise=TRUE,
                             HMM=FALSE,
                             # 运行的线程数目
                             num_threads = 40,
                             write_expr_matrix=TRUE
                             ) 
#----------------------------------------------------------------
# dataset3
#----------------------------------------------------------------                                                                 
cell_epi_3 = cell_epi %>% slice(40001:60000)
dataset3 = rbind(cell_endo, cell_fib, cell_epi_3)
mtx3 = scRNA_matrix[,(which(colnames(scRNA_matrix) %in% dataset3$V1))]
qsave(mtx3, file = './scRNA_matrix.qs')
write.table(dataset3,file = '/root/wangje/Project/刘老师/new_infer_res/dataset3/groupFiles.txt',sep = '\t', row.names=F, col.names=F, quote=F)
write.table(gene, file = '/root/wangje/Project/刘老师/new_infer_res/dataset3/geneFiles.txt',sep = '\t', row.names=F, col.names=F, quote=F)

## 运行infercnv
scRNA_matrix = qread('/root/wangje/Project/刘老师/new_infer_res/dataset3/scRNA_matrix.qs')
geneFile = '/root/wangje/Project/刘老师/new_infer_res/dataset3/geneFiles.txt'
groupFiles = '/root/wangje/Project/刘老师/new_infer_res/dataset3/groupFiles.txt'
ref_group_names = c("Endothelials","Fibroblasts")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=scRNA_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=ref_group_names) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="/root/wangje/Project/刘老师/new_infer_res/dataset3/", 
                             cluster_by_groups=TRUE, # 聚类
                             denoise=TRUE,
                             HMM=FALSE,
                             # 运行的线程数目
                             num_threads = 40,
                             write_expr_matrix=TRUE
                             ) 
#------------------------------------------------------------------
# dataset4
#------------------------------------------------------------------                             
cell_epi_4 = cell_epi %>% slice(60001:80000)
dataset4 = rbind(cell_endo, cell_fib, cell_epi_4)
mtx4 = scRNA_matrix[,(which(colnames(scRNA_matrix) %in% dataset4$V1))]
qsave(mtx4, file = './scRNA_matrix.qs')
write.table(dataset4,file = '/root/wangje/Project/刘老师/new_infer_res/dataset4/groupFiles.txt',sep = '\t', row.names=F, col.names=F, quote=F)
write.table(gene, file = '/root/wangje/Project/刘老师/new_infer_res/dataset4/geneFiles.txt',sep = '\t', row.names=F, col.names=F, quote=F)

## 运行infercnv
scRNA_matrix = qread('/root/wangje/Project/刘老师/new_infer_res/dataset4/scRNA_matrix.qs')
geneFile = '/root/wangje/Project/刘老师/new_infer_res/dataset4/geneFiles.txt'
groupFiles = '/root/wangje/Project/刘老师/new_infer_res/dataset4/groupFiles.txt'
ref_group_names = c("Endothelials","Fibroblasts")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=scRNA_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=ref_group_names) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="/root/wangje/Project/刘老师/new_infer_res/dataset4/", 
                             cluster_by_groups=TRUE, # 聚类
                             denoise=TRUE,
                             HMM=FALSE,
                             # 运行的线程数目
                             num_threads = 40,
                             write_expr_matrix=TRUE
                             ) 

#----------------------------------------------------------------
# dataset5
#----------------------------------------------------------------
cell_epi_5 = cell_epi %>% slice(80001:100000)
dataset5 = rbind(cell_endo, cell_fib, cell_epi_5)
mtx5 = scRNA_matrix[,(which(colnames(scRNA_matrix) %in% dataset5$V1))]
qsave(mtx5, file = '/root/wangje/Project/刘老师/new_infer_res/dataset5/scRNA_matrix.qs')
write.table(dataset5,file = '/root/wangje/Project/刘老师/new_infer_res/dataset5/groupFiles.txt',sep = '\t', row.names=F, col.names=F, quote=F)
write.table(gene, file = '/root/wangje/Project/刘老师/new_infer_res/dataset5/geneFiles.txt',sep = '\t', row.names=F, col.names=F, quote=F)
## 运行infercnv
scRNA_matrix = qread('/root/wangje/Project/刘老师/new_infer_res/dataset5/scRNA_matrix.qs')
geneFile = '/root/wangje/Project/刘老师/new_infer_res/dataset5/geneFiles.txt'
groupFiles = '/root/wangje/Project/刘老师/new_infer_res/dataset5/groupFiles.txt'
ref_group_names = c("Endothelials","Fibroblasts")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=scRNA_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=ref_group_names) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="/root/wangje/Project/刘老师/new_infer_res/dataset5/", 
                             cluster_by_groups=TRUE, # 聚类
                             denoise=TRUE,
                             HMM=FALSE,
                             # 运行的线程数目
                             num_threads = 40,
                             write_expr_matrix=TRUE
                             ) 

#---------------------------------------------------------------
# dataset6
#---------------------------------------------------------------
cell_epi_6 = cell_epi %>% slice(100001:111181)
dataset6 = rbind(cell_endo, cell_fib, cell_epi_6)
mtx6 = scRNA_matrix[,(which(colnames(scRNA_matrix) %in% dataset6$V1))]
qsave(mtx6, file = '/root/wangje/Project/刘老师/new_infer_res/dataset6/scRNA_matrix.qs')
write.table(dataset6,file = '/root/wangje/Project/刘老师/new_infer_res/dataset6/groupFiles.txt',sep = '\t', row.names=F, col.names=F, quote=F)
write.table(gene, file = '/root/wangje/Project/刘老师/new_infer_res/dataset6/geneFiles.txt',sep = '\t', row.names=F, col.names=F, quote=F)
## 运行infercnv
scRNA_matrix = qread('/root/wangje/Project/刘老师/new_infer_res/dataset6/scRNA_matrix.qs')
geneFile = '/root/wangje/Project/刘老师/new_infer_res/dataset6/geneFiles.txt'
groupFiles = '/root/wangje/Project/刘老师/new_infer_res/dataset6/groupFiles.txt'
ref_group_names = c("Endothelials","Fibroblasts")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=scRNA_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=ref_group_names) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="/root/wangje/Project/刘老师/new_infer_res/dataset6/", 
                             cluster_by_groups=TRUE, # 聚类
                             denoise=TRUE,
                             HMM=FALSE,
                             # 运行的线程数目
                             num_threads = 40,
                             write_expr_matrix=TRUE
                             ) 





















