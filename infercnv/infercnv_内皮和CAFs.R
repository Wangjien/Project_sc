rm(list=ls())
gc()
setwd("/root/wangje/Project/刘老师/infercnv_使用内皮和CAFs_全部上皮")
library(Seurat)
library(infercnv)
library(tidyverse)
# 读入文件（上皮细胞）
load("/root/wangje/Project/刘老师/EPi/新结果/添加肝脏/new/new_FastMNN.RData")
EPI <- scRNA_seurat
rm(scRNA_seurat)
# 读入文件（内皮细胞）
load("/root/wangje/wangjien_new/liu-big/scRNA_FastMNN_new.RData")
scRNA@meta.data %>% 
    group_by(celltype) %>% 
    count(celltype)
#   celltype           n
#   <chr>          <int>
# 1 B cells         5273
# 2 Endothelials   21441
# 3 Epithelials   111181
# 4 Fibroblasts    28336
# 5 Hepatocytes     3036
# 6 Keratinocytes  10230
# 7 Myeloids       87482
# 8 NK & T cells  101962
# 9 Plasame cells   6397
Endo <- scRNA[,scRNA$celltype %in% "Endothelials"] # nolint
# 读入文件（成纤维细胞）
CAFs <- scRNA[,scRNA$celltype %in% "Fibroblasts"]
# 去除相同的细胞
Express <- rownames(Endo@meta.data) %in% rownames(EPI@meta.data)[(rownames(EPI@meta.data) %in% rownames(Endo@meta.data))]
Endo_sub <- Endo[,!Express]

table(rownames(EPI@meta.data) %in%  rownames(CAFs@meta.data))
Express <- rownames(CAFs@meta.data) %in% rownames(EPI@meta.data)[(rownames(EPI@meta.data) %in% rownames(CAFs@meta.data))]
CAFs_sub <- CAFs[,!Express]

# 合并数据
tmp <- merge(EPI,Endo_sub)
MergeResult <- merge(CAFs_sub,tmp)
rm(tmp);gc()

### 运行infercnv
MergeResult@meta.data[which(is.na(MergeResult$celltype)),"celltype"] <- "Epithelials"
celltype <- unique(MergeResult$celltype)
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
step01(MergeResult,celltype = celltype)

## infercnv第二步
ref_group_names = c("Endothelials","Fibroblasts")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=scRNA_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=ref_group_names) 
save(infercnv_obj,file = "./infercnv_obj.RData") 
## Infercnv第三步
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="./", 
                             cluster_by_groups=TRUE, # 聚类
                             denoise=TRUE,
                             HMM=FALSE,
                             # 运行的线程数目
                             num_threads = 25) 
## ---------------------- 安装比例进行筛选 ----------------------------
test <- MergeResult@meta.data  %>% 
            group_by(celltype) %>%
            count(celltype)
test$Prop <- apply(test[2],1,function(x) x/sum(test[[2]]))
test$num <- round(50000*test[[3]])
head(test)  
# celltype          n  Prop   num
#   <chr>         <int> <dbl> <dbl>
# 1 Endothelials  21344 0.134  6679
# 2 Epithelials  110447 0.691 34560
# 3 Fibroblasts   27998 0.175  8761
EPI_sub <- EPI[,rownames(EPI@meta.data) %in% sample(rownames(EPI@meta.data),34560)]
EPI_sub$celltype <- "Epithelials"
CAFs_sub <- CAFs_sub[,rownames(CAFs_sub@meta.data) %in% sample(rownames(CAFs_sub@meta.data),8761)]
Endo_sub <- Endo_sub[,rownames(Endo_sub@meta.data) %in% sample(rownames(Endo_sub@meta.data),6679)]
tmp <- merge(EPI_sub,CAFs_sub)
MergeResult <- merge(Endo_sub,tmp)
step01(data = MergeResult,celltype = celltype)
setwd("/root/wangje/Project/刘老师/infercnv_使用内皮和CAFs_部分上皮")
