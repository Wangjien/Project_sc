#----------------------------------------------------------------
# 将FastMNN中的Myeloids提取后再分析
#----------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(DoubletFinder)
# load data 
big.fastmnn = '/root/wangje/wangjien_new/liu-big/scRNA_FastMNN_new.RData'
load(big.fastmnn)
scRNA_sub = scRNA[,scRNA$celltype == 'Myeloids']

#----------------------------------------------------------------
# 使用FastMNN进行分析
#----------------------------------------------------------------

scRNAlist <- SplitObject(scRNA_sub,split.by="sample")
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x))

scRNA <- RunFastMNN(object.list = scRNAlist)
scRNA <- RunUMAP(scRNA, reduction = "mnn", dims = 1:50)
scRNA <- RunTSNE(scRNA, reduction = "mnn", dims = 1:50)
scRNA <- FindNeighbors(scRNA, reduction = "mnn", dims = 1:50)
scRNA <- FindClusters(scRNA,resolution=seq(0.05,1,0.05))

# 保存图片


#----------------------------------------------------------------
# 使用Monocle3进行聚类分析
#----------------------------------------------------------------
#! 获取表达矩阵
expression_matrix = scRNA@assays$RNA@data
cell_metadata = scRNA@meta.data
gene_metadata = as.data.frame(rownames(scRNA))



#----------------------------------------------------------------
# 使用Scnpy的BBKNN进行聚类
#----------------------------------------------------------------
# 01 Seurat to Anndata
library(Seurat)
library(dior)
library(optparse)

option_list = list(
   make_option(c('-f','--file'),type = 'character', default = FALSE, action = 'store', help = 'input files')
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# 访问解析的命令行参数
input_files = opt$file
if(endsWith(input_files, 'qs')){
    scRNA = qread()
}else if(endsWith(input_files, 'rds')){
    scRNA = readRDS()
}else if(endsWith(input_files, 'RDS')){
    scRNA = readRDS()
}else if(endsWith(input_files, '')){

}


## 添加信息
id1 <- c("CJME","CMDI","HDYU","HXZH","LCLI","WYJU","WZLA","ZXME","ZJLI_0116")
id2 <- c("CZYI","FHXHBS1","FYYI","HEJI","LAWE","ZEYI","ZFXI","LIPE")
id3 <- c("CJME_0707","CMDI_0624","HDYU_0720","HXZH_0220","LCLI_0623","WYJU_0122","ZXME_0223","ZJLI_0312")
scRNA[['Patient']] <- unlist(strsplit(rownames(scRNA@meta.data), split = '_[A|T|C|G]*$'))
scRNA@meta.data$Treat_assess <-ifelse(scRNA@meta.data$Patient %in% id1, "R_Pre", 
                                    ifelse(scRNA@meta.data$Patient %in% id2,"NR_Pre",
                                    ifelse(scRNA@meta.data$Patient %in% id3,"R_Post","NR_Post")))

ids4 <- c('CJME','CJME_0707','CMDI','CMDI_0624','HDYU','HDYU_0720','HXZH','HXZH_0220','LCLI','LCLI_0623','WYJU','WYJU_0122','ZXME','ZXME_0223','ZJLI_0116','ZJLI_0312',
'CZYI','CZYI_0702','FHXHBS1','FHXHBS2','HEJI','HEJX','LAWE','LAWE_0309','ZEYI','ZEYI_0204','WZLA','FYYI','ZFXI','LIPE')
ids5 <- c('RCJMEprePRLymph','RCJMEpostPRLymph','RCMDIprePRLymph','RCMDIpostPRLymph','RHDYUprePRBreast','RHDYUpostPRBrest','RHXZHprePRBreast','RHXZHpostPRBreast','RLCLIpreCRLymph','RLCLIpostChestWall','RWYJUprePRLymph','RWYJUpostPRLymph','RZXMEprePRBreast','RZXMEpostPRBreast','RZJLIprePRBreast','RZJLIpostPRBreast',
'NRCZYIprePDLiver','NRCZYIpostPDLiver','NRFHXHBSpreSDBreast','NRFHXHBSpostSDBreast','NRHEJIpreSDBreast','NRHEJXpostSDBreast','NRLAWEprePDLiver','NRLAWEpostPDLiver','NRZEYIprePDBreast','NRZEYIpostPDBreast',
'RRWZLApreChestWall','NRFYYIpreSDLiver','NRZFXIprePDBreast','NRLIPEprePDLiver')
ids6 <- c('Lymph','Lymph','Lymph','Lymph','Breast','Breast','Breast','Breast','Lymph','ChestWall','Lymph','Lymph','Breast','Breast','Breast','Breast',
'Liver','Liver','Breast','Breast','Breast','Breast','Liver','Liver','Breast','Breast','ChestWall','Liver','Breast','Liver')
ids7 <- c('R','R','R','R','R','R','R','R','R','R','R','R','R','R','R','R','NR','NR','NR','NR','NR','NR','NR','NR','NR','NR','R','NR','NR','NR')

for(i in 1:length(ids4)){
    print(ids4[i])
    print(ids5[i])
    print(ids6[i])
    print(ids7[i])
   scRNA@meta.data[which(scRNA$Patient == ids4[i]),"Rename"]  <- ids5[i]
   scRNA@meta.data[which(scRNA$Patient == ids4[i]),"Tissue"]  <- ids6[i]
   scRNA@meta.data[which(scRNA$Patient == ids4[i]),"Rseponse"]  <- ids7[i]
}

scRNA$Treat_assess = factor(scRNA$Treat_assess,levels=c('R_Pre','R_Post', 'NR_Pre', 'NR_Post'))
scRNA$Rseponse = factor(scRNA$Rseponse,levels =c('R','NR'))

