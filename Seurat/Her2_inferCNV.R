

library("dplyr")
library("Seurat")
library("SeuratObject")
library("ggplot2")
library("infercnv")
library("stringr")
library("qs")

# 生成参考文件
scRNA <- qread('/root/wangje/Project/刘老师/大群/FastMNN.qs')
table(scRNA$celltype)
# Adipocytes           B cells   Dendritic cells Endothelial cells  Epithelial cells       Fibroblasts       Hepatocytes     Keratinocytes
#      2071             12185              2526             21243            105326             26637              3301             10222
#     Macro|Mono           Neurons      NK & T cells
#         84661              5776            101390

ref1 <- scRNA[,rownames(scRNA@meta.data) %in% sample(Cells(subset(scRNA,subset = celltype == "Endothelial cells")),replace = FALSE,size = 2124)]
ref2 <- scRNA[,rownames(scRNA@meta.data) %in% sample(Cells(subset(scRNA,subset = celltype == "Fibroblasts")),replace = FALSE,size = 2664)]
ref <- merge(ref1,ref2)
qsave(ref, file = "./Ref_seurat.qs")


### 非其他细胞
###----------------------------------------------------------------------------------------------------------------------
# 读入文件
setwd('/root/wangje/Project/刘老师/大群/Epithelials/')
scRNA1 <- qread('./Epi_MNN.qs')
scRNA1
# 细胞命名
# scRNA@meta.data <- scRNA@meta.data %>% dplyr::mutate(
#     celltype = case_when(
#         RNA_snn_res.0.2 %in% c(1) ~ 'Cycling Epi',
#         RNA_snn_res.0.2 %in% c(7) ~ 'Low HER2+ custer1',
#         RNA_snn_res.0.2 %in% c(5,6) ~ 'Low HER2+ cluster2',
#         RNA_snn_res.0.2 %in% c(8,0,2,3) ~ 'Malignant epithelial cells',
#         RNA_snn_res.0.1 %in% c(4) ~ 'other cluaster1',
#         RNA_snn_res.0.1 %in% c(3) ~ 'other cluster2'
#     )
# )
scRNA1@meta.data <- scRNA1@meta.data  %>% dplyr::mutate(
    celltype = case_when(
        RNA_snn_res.0.3 %in% c(3,7) ~ "Malignant Epi 2",
        RNA_snn_res.0.3 %in% c(1) ~ 'Cycling Epi',
        RNA_snn_res.0.3 %in% c(4) ~ "other Epi 1",
        RNA_snn_res.0.3 %in% c(10,11) ~ 'other Epi 2',
        TRUE ~ "Malignant Epi 1"))


scRNA1 <- scRNA1[,!is.na(scRNA1$celltype)]
table(scRNA1$celltype)
#  Cycling Epi Malignant Epi 1 Malignant Epi 2     other Epi 1     other Epi 2
#     20410           69429           15495           10398            3117

# 提取quest细胞
allCells <- c(sample(Cells(subset(scRNA1,subset = celltype == "Malignant Epi 1")),replace = FALSE,size = 6943),
    sample(Cells(subset(scRNA1,subset = celltype == "Malignant Epi 2")),replace = FALSE,size = 1550),
    sample(Cells(subset(scRNA1,subset = celltype == "Cycling Epi")),replace = FALSE,size = 2041))
quest1 <- scRNA1[,rownames(scRNA1@meta.data) %in% allCells]

## -----> 将SeuratV5结构降为Seurat V4结构 <------ ##
quest1.meta <- quest1@meta.data
quest1 <- CreateSeuratObject(GetAssayData(quest1,layer = 'counts'))
quest1 <- AddMetaData(quest1, quest1.meta)

ref <- qread("/root/wangje/Project/刘老师/大群/Epithelials/PNG/ref_seurat.qs")
ref.meta <- ref@meta.data
ref <- CreateSeuratObject(GetAssayData(ref,layer = 'counts'))
ref <- AddMetaData(ref, ref.meta)

# 合并quest和Ref数据
sce <- merge(ref, quest1)
# sce <-JoinLayers(sce)
# 1) 表达矩阵文件
saveRDS(as.matrix(sce[["RNA"]]$counts), "sce_CNV2_counts.matrix")
matrix_counts <- readRDS("sce_CNV2_counts.matrix")
matrix_counts <- matrix_counts[,names(sce$celltype)]

# 2）celltype文件
write.table(sce$celltype, "celltype2.label.txt", sep = "\t", quote = F, col.names = F)


infercnv_obj = CreateInfercnvObject(raw_counts_matrix = matrix_counts,#count矩阵
                                    annotations_file="./celltype2.label.txt",#celltype信息
                                    delim="\t",
                                    gene_order_file="/root/wangje/software/InferCNV/hg38_gencode_v27.txt",
                                    ref_group_names=c("Endothelial cells","Fibroblasts"),#HC组正常的细胞作为reference
                                    chr_exclude=c("chrY", "chrM"))#选择不需要的染色体，查看帮助函数去除

#具体参数可参照帮助函数，一般使用默认或者看看参考文献
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="Epithe",  #分析结果out文件夹名
                             no_prelim_plot = T,
                             cluster_by_groups=T,
                             denoise=TRUE,
                             HMM=F,
                             min_cells_per_gene = 10,
                             num_threads=5,#线程数
                             write_expr_matrix = T
                             #这里要选择T，分析后将结果exp矩阵导出，出现infercnv.observations.txt结果文件
                             # 很多网上好几年前的教程，没有这个。那是因为inferCNV更新了，这里的默认参数是F
                             )

### 其他细胞
###------------------------------------------------------------------------------------------------------------
# 提取出其他细胞
quest2 <- scRNA1[,scRNA1$celltype %in% c('other Epi 1','other Epi 2')]
quest2.meta <- quest2@meta.data
quest2 <- CreateSeuratObject(GetAssayData(quest2,layer = 'counts'))
quest2 <- AddMetaData(quest2, quest2.meta)

ref <- qread("/root/wangje/Project/刘老师/大群/Epithelials/PNG/ref_seurat.qs")
ref.meta <- ref@meta.data
ref <- CreateSeuratObject(GetAssayData(ref,layer = 'counts'))
ref <- AddMetaData(ref, ref.meta)

# 合并文件
sce <- merge(ref, quest2)

# 1) 表达矩阵文件
saveRDS(as.matrix(sce[["RNA"]]$counts), "sce_CNV_counts.matrix")
matrix_counts <- readRDS("sce_CNV_counts.matrix")
matrix_counts <- matrix_counts[,names(sce$celltype)]

# 2）celltype文件
write.table(sce$celltype, "celltype.label.txt", sep = "\t", quote = F, col.names = F)


infercnv_obj = CreateInfercnvObject(raw_counts_matrix = matrix_counts,#count矩阵
                                    annotations_file="./celltype.label.txt",#celltype信息
                                    delim="\t",
                                    gene_order_file="/root/wangje/software/InferCNV/hg38_gencode_v27.txt",
                                    ref_group_names=c("Endothelial cells","Fibroblasts"),#HC组正常的细胞作为reference
                                    chr_exclude=c("chrY", "chrM"))#选择不需要的染色体，查看帮助函数去除


#具体参数可参照帮助函数，一般使用默认或者看看参考文献
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="otherEpi",  #分析结果out文件夹名
                             no_prelim_plot = T,
                             cluster_by_groups=T,
                             denoise=TRUE,
                             HMM=F,
                             min_cells_per_gene = 10,
                             num_threads=5,#线程数
                             write_expr_matrix = T
                             #这里要选择T，分析后将结果exp矩阵导出，出现infercnv.observations.txt结果文件
                             # 很多网上好几年前的教程，没有这个。那是因为inferCNV更新了，这里的默认参数是F
                             )


##=============================================================================================================================================
## 对所有的上皮细胞进行分析
##=============================================================================================================================================
for(sample in unique(Epi$Rename)){
    tryCatch({
        if(!dir.exists(sample)){dir.create(file.path("/root/wangje/Project/刘老师/大群/Epithelials/Infercnv",sample))}
        message(sprintf("------> %s <------",sample))
        setwd(file.path("/root/wangje/Project/刘老师/大群/Epithelials/Infercnv",sample))
        getwd() # 查看目录
        # 合并文件
        sce <- merge(ref,Epi[,Epi$Rename == sample])
        saveRDS(as.matrix(sce[["RNA"]]$counts), "sce_CNV_counts.matrix")
        matrix_counts <- readRDS("sce_CNV_counts.matrix")
        matrix_counts <- matrix_counts[,names(sce$celltype)]
        write.table(sce$celltype, "celltype.label.txt", sep = "\t", quote = F, col.names = F) # 输出celltype
        # 生成Infercnv第一步
        infercnv_obj = CreateInfercnvObject(
            raw_counts_matrix = matrix_counts,#count矩阵
            annotations_file="./celltype.label.txt",#celltype信息
            delim="\t",
            gene_order_file="/root/wangje/software/InferCNV/hg38_gencode_v27.txt",
            ref_group_names=c("Endothelial cells","Fibroblasts"),# reference
            chr_exclude=c("chrY", "chrM"))#选择不需要的染色体，查看帮助函数去除
        # 保存第一步文件                            
        qsave(infercnv_obj, file = file.path("/root/wangje/Project/刘老师/大群/Epithelials/Infercnv",sample,"001_infercnv_obj.qs"))
        # 生成infercnv第二步
        infercnv_obj = infercnv::run(
            infercnv_obj,
            cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
            out_dir="otherEpi",  #分析结果out文件夹名
            no_prelim_plot = T,
            cluster_by_groups=T,
            denoise=TRUE,
            HMM=F,
            min_cells_per_gene = 10,
            num_threads=20,#线程数
            write_expr_matrix = T)
            #这里要选择T，分析后将结果exp矩阵导出，出现infercnv.observations.txt结果文件
            # 很多网上好几年前的教程，没有这个。那是因为inferCNV更新了，这里的默认参数是F
        # 保存infercnv第二步的文件
        qsave(infercnv_obj, file = file.path("/root/wangje/Project/刘老师/大群/Epithelials/Infercnv",sample,"002_infercnv_obj.qs"))
    },error = function(e){
        message(e$error)
        return(NULL)
    },finally = {
        message(sprintf("%s finish....",sample))
    })
}

##-------------------------------------------------------------------------------------------------------------
## 提取数据进行Kmeans分群
setwd("/root/wangje/Project/刘老师/大群/Epithelials/Infercnv")
Dir <- "/root/wangje/Project/刘老师/大群/Epithelials/Infercnv"
out_dir <- '/root/wangje/Project/刘老师/大群/Epithelials/Infercnv/result/'
result <- list()
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(infercnv)
library(patchwork)

for(sample in list.files(path = '.', pattern = "^R|^NR",full.names = T)){
    message(sample)
    # 读入数据
    cat("数据地址: ",file.path(sample,"002_infercnv_obj.qs"))
    infercnv_obj <- qs::qread(file.path(sample,"002_infercnv_obj.qs"))
    tryCatch({
    expr <- infercnv_obj@expr.data
    #top注释-染色体
    gene_pos <- read.delim("/root/wangje/software/InferCNV/hg38_gencode_v27.txt", header = F)
    gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
    new_cluster <- unique(gene_pos$V2)
    top_color <- HeatmapAnnotation(
        cluster = anno_block(
            labels = gsub("chr", "", new_cluster),
            gp = gpar(col = "white"),
            labels_gp = gpar(cex = 1, col = "black"),
            height = unit(5,"mm"))) 
    #提取reference cell位置
    ref_cell <- c(infercnv_obj@reference_grouped_cell_indices$`Endothelial cells`,
                infercnv_obj@reference_grouped_cell_indices$`Fibroblasts`)
    #提取测试细胞的位置
    obs_cell <- c(infercnv_obj@observation_grouped_cell_indices$`Epithelail cells`)
    # 添加注释信息
    cell_anno =data.frame(cell_id =c(colnames(expr)[ref_cell],colnames(expr)[obs_cell]),
                      group =c(rep("Ref",length(ref_cell)),rep("Epithelial",length(obs_cell))))
    print("01")
    # 进行kmeans聚类分群
    set.seed(123)
    kmeans.result <- kmeans(t(expr), 8)
    kmeans_df <- data.frame(kmeans.result$cluster)   
    colnames(kmeans_df) <- "k_cluster"
    kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
    kmeans_df=kmeans_df%>%inner_join(cell_anno,by="cell_id")%>%arrange(k_cluster)
    kmeans_df$k_cluster=as.factor(kmeans_df$k_cluster) 
    #注释
    # print(kmeans_df)
    annotation_row = data.frame(
    k_cluster =kmeans_df$k_cluster,
    group = kmeans_df$group)
    row.names(annotation_row) <- kmeans_df$cell_id
    color_cluster=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#4dac8e","#d75030","#834127")
    names(color_cluster)=as.character(1:8)
    left_anno <- rowAnnotation(df = annotation_row,
                            col=list(group=c("Ref"="#00A0877F","Epithelial" = "#E64B357F"),
                            k_cluster=color_cluster),
                            show_annotation_name = F)
    #绘制热图
    filename = paste0(paste0(out_dir,stringr::str_split_fixed(sample,pattern = '\\./',n=2)[,2],'_CNV_haetmap.png'))
    png(file =filename,width = 15,height = 10,units = 'in', res = 300)
    ht = Heatmap(t(log2(expr))[rownames(annotation_row),],
                col = colorRamp2(c(-0.5,0,0.5), c("#2166ac","white","#b2182b")),
                cluster_rows = F,
                cluster_columns = F,
                show_column_names = F,
                show_row_names = F,
                column_split = factor(gene_pos$V2, new_cluster),
                heatmap_legend_param = list(title = "inferCNV",
                                            direction = "vertical",
                                            title_position = "leftcenter-rot",
                                            legend_height = unit(3, "cm")),
                
                left_annotation = left_anno, 
                row_title = NULL,
                column_title = NULL,
                top_annotation = top_color,
                border = T)
    # draw(ht)
    draw(ht, heatmap_legend_side = "right")
    dev.off()
    qs::qsave(kmeans.result,paste0(paste0(out_dir,stringr::str_split_fixed(sample,pattern = '\\./',n=2)[,2],'_kmeans.result.qs')))
    # 保存kmeans结果
    result[[stringr::str_split_fixed(sample,'\\./',n=2)[,2]]] <- kmeans.result
    },error = function(e){
        message(e$error)
        return(NULL)
    },finally = {
        message("finish... ...")
    })
}


##-------------------------------------------------------------------------------------------------------------
## 确定出正常细胞和恶性肿瘤细胞
# 每个样本kmeans聚类结果，提取出正常上皮细胞和非正常上皮细胞
cluster_list <- list(
    NRCZYIpostPD_Liver=c(1,2,5,7,8),
    NRCZYIprePD_Live=c(1,4,5,7,8),
    NRFHXHBSpostSDBreast=c(1,5,7,8),
    NRFHXHBSpreSDBreast=c(1,2,5,8),
    NRHEJIpreSDBreast =c(1,2,4,5,7),
    NRHEJXpostSDBreast=c(3,8))

# head(as.data.frame(result[[1]]$cluster))                                             
#                                    result[[1]]$cluster                                   
# CJME_0707_AACGCTTAATAGCGACAAACATCG                   2                                   
# CJME_0707_AACGTGATACATTGGCAGCACCTC                   7                                   
# CJME_0707_AACGTGATATAGCGACATGCCTAA                   7                                   
# CJME_0707_AACGTGATCGGATTGCCACCTTAC                   7                                   
# CJME_0707_AAGACGGAAAGAGATCACACGACC                   1                                   
# CJME_0707_AAGACGGACGCATACATATCAGCA                   2                                  
                                                       
##-------------------------------------------------------------------------------------------------------------
plist <- list()
for(sample in list.files(path = '.', pattern = "^R|^NR",full.names = T)){
    sp <- stringr::str_split_fixed(sample,pattern = "\\./",n=2)[,2]
    message(sp)
    # 读入数据
    cat("数据地址: ",file.path(sample,"002_infercnv_obj.qs"))
    infercnv_obj <- qs::qread(file.path(sample,"002_infercnv_obj.qs"))
    tryCatch({
        expr <- infercnv_obj@expr.data
        CNV_score=as.data.frame(colMeans((expr-1)^2))
        colnames(CNV_score)="CNV_score"
        CNV_score$cell_id=rownames(CNV_score)
        ## 获得kmeans的结果。上一步已经把kmeans的结果放在list中
        kmeans_df <- as.data.frame(result[[sp]]$cluster)
        colnames(kmeans_df) <- "k_cluster"
        kmeans_df$cell_id <- rownames(kmeans_df)
        CNV_score=CNV_score%>%inner_join(kmeans_df,by="cell_id")#结合分组和kemans聚类结果
        CNV_score$k_cluster <- as.character(CNV_score$k_cluster)
        print(head(CNV_score))
        ## 对每一个样本的kmeans结果绘制小提琴图
        p <- ggplot(CNV_score, aes(k_cluster,CNV_score))+
        geom_violin(aes(fill=k_cluster),color="black")+
        scale_fill_manual(values = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#E78AC3", "#333399", "#A6D854", "#E5C494"))+
        theme_bw()+
        ggtitle(sp)+
        theme(axis.text.x = element_text(colour = "black", size = 12),
                axis.text.y = element_text(colour = "black", size = 10),
                axis.title = element_text(color = 'black', size = 12),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
        ggpubr::stat_compare_means()+
        stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), 
                    geom = "pointrange", color = "black", size=0.5)+
        stat_summary(fun.y = "mean", geom = "line", aes(group = 1),color="orange",linewidth = 1.1) # 添加趋势线
        print(p)
        plist[[sp]] <- p
    },error = function(e){
        message(e$error)
        return(NULL)
    })
}

## 写出绘图的文件
png('./不同样本kmens小提琴图.png',height = 20,width = 35, units = 'in',res=300)
ggpubr::ggarrange(plotlist = plist, common.legend = T, ncol=7,nrow = 4)
dev.off()

## 挑选出每个样本中的肿瘤细胞
k_cluster_select <- list(
    "NRCZYIpostPDLiver" =c(3,4,6),"NRCZYIprePDLiver" =c(2,3,6),"NRFHXHBSpostSDBreast"=c(3,4,6),
    "NRFHXHBSpreSDBreast"=c(3,4,6,7),"NRHEJIpreSDBreast"=c(3,7),"NRHEJXpostSDBreast"=c(1,3,6,7,8),
    "NRLAWEpostPDLiver"=c(2,3,4,7), "NRLAWEprePDLiver" =c(1,2,3,8),"NRLIPEprePDLiver"=c(2,3,6),
    "NRZEYIpostPDBreast"=c(2,3,4,5,8), "NRZEYIprePDBreast"=c(2,3,5,6),"NRZFXIprePDBreast"=c(1,2,3,4,6,7),
    "RCJMEpostPRLymph"=c(2,3,5,6,8), "RCJMEprePRLymph"=c(2,3,8), "RCMDIprePRLymph"=c(2,3,7,8),
    "RHDYUpostPRBrest"=c(2,3,4,6), "RHDYUprePRBreast"=c(2,3,4,5,6,8), "RHXZHpostPRBreast"=c(2,3,5,6),
    "RHXZHprePRBreast"=c(1,3,4),   "RLCLIpostChestWall"=c(2,5,8), "RRWZLApreChestWall"=c(3,8),
    "RWYJUpostPRLymph"=c(1,2,3,8),  "RWYJUprePRLymph"=c(3,4,5,8), "RZJLIprePRBreast"=c(2,3),
    "RZXMEpostPRBreast"=c(4,6,7),    "RZXMEprePRBreast"=c(1,2,3,5,8)
)
## 提取出恶性肿瘤细胞和正常肿瘤细胞
flist_tumor <- list()
flist_normal <- list()

for(sample in list.files(path = '.', pattern = "^R|^NR",full.names = T)){
    sp <- stringr::str_split_fixed(sample,pattern = "\\./",n=2)[,2]
    message(sp)
    kmeans_df <- as.data.frame(result[[sp]]$cluster)
    colnames(kmeans_df) <- "k_cluster"
    kmeans_df$cell_id <- rownames(kmeans_df)
    flist_tumor[[sp]] <- kmeans_df %>% dplyr::filter(k_cluster %in% k_cluster_select[[sp]]) %>% dplyr::pull(cell_id)
    flist_normal[[sp]] <-  kmeans_df %>% dplyr::filter(!k_cluster %in% k_cluster_select[[sp]]) %>% dplyr::pull(cell_id)
}

## 分样本绘制肿瘤细胞和非肿瘤细胞以及参考细胞的热图
# dir.create("/root/wangje/Project/刘老师/大群/Epithelials/Infercnv/result/Exper")
setwd('/root/wangje/Project/刘老师/大群/Epithelials/Infercnv/')

out_dir <- "/root/wangje/Project/刘老师/大群/Epithelials/Infercnv/result/Exper"
flist2 <- list()
for(sample in list.files(path = '.', pattern = "^R|^NR",full.names = T)){
    sp <- stringr::str_split_fixed(sample,pattern = "\\./",n=2)[,2]
    message(sp)
    # 读入数据
    cat("数据地址: ",file.path(sample,"002_infercnv_obj.qs"))
    infercnv_obj <- qs::qread(file.path(sample,"002_infercnv_obj.qs"))
    print(dim(infercnv_obj@expr.data))
    # 提取其中的肿瘤细胞和非肿瘤上皮细胞以及参考细胞
    cells1 <-colnames(infercnv_obj@expr.data)[infercnv_obj@observation_grouped_cell_indices$`Epithelail cells`][colnames(infercnv_obj@expr.data)[infercnv_obj@observation_grouped_cell_indices$`Epithelail cells`] %in% flist_tumor[[sp]]]
    tmuor <- infercnv_obj@expr.data[, cells1]
    ref_cell <- c(infercnv_obj@reference_grouped_cell_indices$`Endothelial cells`,
            infercnv_obj@reference_grouped_cell_indices$`Fibroblasts`)
    Ref <- infercnv_obj@expr.data[,ref_cell]
    # 预测的正常上皮细胞
    cells2 <-colnames(infercnv_obj@expr.data)[infercnv_obj@observation_grouped_cell_indices$`Epithelail cells`][!colnames(infercnv_obj@expr.data)[infercnv_obj@observation_grouped_cell_indices$`Epithelail cells`] %in% flist_tumor[[sp]]]
    epi <- infercnv_obj@expr.data[,cells2]
    # 生成统计文件
    tmp <- data.frame(
        sample =sp,
        `Malignant Epithelial Cells Num` = dim(tmuor)[2],
        `Normal Epithelial Cells Num` = dim(epi)[2],
        ReferenceCellsNum = length(ref_cell))
    flist2[[sp]] <- tmp
    # 保存上面的矩阵文件
    qs::qsave(tmuor, file = file.path(out_dir,paste0(sp,'_tmour_exper.qs')))
    qs::qsave(epi, file = file.path(out_dir,paste0(sp,'_epi_exper.qs')))
    qs::qsave(Ref, file = file.path(out_dir,paste0(sp,'_Ref_exper.qs')))
}

##　写出结果文件
df <- do.call(what = "rbind",flist2)
rownames(df) <- NULL
df[nrow(df) + 1, ] <- list("Sum", sum(df[[2]]),sum(df[[3]]), NA) # 增加一行总和
# 绘制表格
tbody.style = tbody_style(color = "black",
   fill = c("#e8f3de", "#d3e8bb"), hjust=1, x=0.9)
ggtexttable(df, rows = NULL,
           theme = ttheme(
             colnames.style = colnames_style(color = "white", fill = "#8cc257"),
             tbody.style = tbody.style)) -> t_table
# 保存统计表格
ggsave(filename = './kmeans区分恶性细胞统计.png', height = 14, width = 12, 
    plot = t_table,limitsize = FALSE,bg = "white")

## ----------------------------------------------------------------------------------------------------
## 绘制热图
library(ComplexHeatmap)
library(dplyr)
library(patchwork)
library(ggplot2)
library(infercnv)

setwd('/root/wangje/Project/刘老师/大群/Epithelials/Infercnv/result/Exper')
flist <- list()
for(i in list.files(path = './', pattern = '_epi_exper.qs$')){
    flist[[i]] <- qs::qread(i)
}
# 修改名称
names(flist) <- stringr::str_split_fixed(names(flist),'_',n= 5)[,1]

flist2 <- list()
for(i in list.files(path = './', pattern = '_tmour_exper.qs$')){
    flist2[[i]] <- qs::qread(i)
}
# 修改名称
names(flist2) <- stringr::str_split_fixed(names(flist2),'_',n= 5)[,1]

flist3 <- list()
for(i in list.files(path = './', pattern = '_Ref_exper.qs$')){
    flist3[[i]] <- qs::qread(i)
}
# 修改名称
names(flist3) <- stringr::str_split_fixed(names(flist3),'_',n= 5)[,1]
ids5 <- c('RCJMEprePRLymph','RCJMEpostPRLymph','RCMDIprePRLymph','RCMDIpostPRLymph','RHDYUprePRBreast',
            'RHDYUpostPRBrest','RHXZHprePRBreast','RHXZHpostPRBreast','RLCLIpreCRLymph','RLCLIpostChestWall',
            'RWYJUprePRLymph','RWYJUpostPRLymph','RZXMEprePRBreast','RZXMEpostPRBreast','RZJLIprePRBreast',
            'RZJLIpostPRBreast','NRCZYIprePDLiver','NRCZYIpostPDLiver','NRFHXHBSpreSDBreast','NRFHXHBSpostSDBreast',
            'NRHEJIpreSDBreast','NRHEJXpostSDBreast','NRLAWEprePDLiver','NRLAWEpostPDLiver','NRZEYIprePDBreast',
            'NRZEYIpostPDBreast','RRWZLApreChestWall','NRFYYIpreSDLiver','NRZFXIprePDBreast','NRLIPEprePDLiver')

ids6 <- c(
    "P1-Pre", "P1-Post", "P2-Pre", "P2-Post", "P3-Pre", "P3-Post", "P4-Pre", "P4-Post", "P5-Pre", "P5-Post", "P6-Pre", "P6-Post",
    "P7-Pre", "P7-Post", "P8-Pre", "P8-Post", "P9-Pre","P10-Pre", "P10-Post", "P11-Pre", "P11-Post", "P12-Pre", "P12-Post", "P13-Pre", "P13-Post",
    "P14-Pre", "P14-Post", "P15-Pre", "P16-Pre", "P17-Pre")

lst <- list()
for (i in 1:length(ids5)) {
    lst[[ids5[i]]] <- ids6[[i]]
}

Names <- vector(mode = "character",length = length(flist))
for(i in 1:length(flist)){
    Names[i] <- lst[[names(flist[i])]]
}

names(flist) <- Names
names(flist2) <- Names
names(flist3) <- Names


## 使用RColorBrawser颜色绘制热图
library(RColorBrewer)
colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100)
# colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(20)
values <- seq(-0.5, 0.5, length.out = 101)[-101]
col_fun <- colorRamp2(values, colors)

htlist <- list()
for(i in names(flist)){
    print(i)
    test <- cbind(flist3[[i]],flist[[i]],flist2[[i]]) # 合并文件
    print(dim(test))
    gene_pos <- read.delim("/root/wangje/software/InferCNV/hg38_gencode_v27.txt", header = F) # 基因坐标信息
    cell_anno =data.frame(cell_id =colnames(test),
                        group =c(rep("Reference",dim(flist3[[i]])[2]),rep("Normal epithelial cells",dim(flist[[i]])[2]),rep('Malignant epithelial cells',dim(flist2[[i]])[2]))) # 细胞注释信息
    gene_pos <- gene_pos[gene_pos$V1 %in% rownames(test),]
    new_cluster <- unique(gene_pos$V2)
    top_color <- HeatmapAnnotation(
        cluster = anno_block(
            labels = gsub("chr", "", new_cluster),
            gp = gpar(col ="white"),
            labels_gp = gpar(cex = 1, col = "black"),
            labels_rot = 45,
            height = unit(5,"mm"))) 

    annotation_row = data.frame(celltype = c(rep("Reference",dim(flist3[[i]])[2]),rep("Normal epithelial cells",dim(flist[[i]])[2]),rep('Malignant epithelial cells',dim(flist2[[i]])[2])))
    rownames(annotation_row) <- colnames(test)
    left_anno <- rowAnnotation(df = annotation_row,
                            col=list(ncol=3,celltype=c("Reference"="#2166ac","Normal epithelial cells" = "#b2182b","Malignant epithelial cells"="#ce9d3f"), title_position = "topcenter"),
                            show_annotation_name = F)
    png(paste0('003_',i,"_CNV_heatmap.png"),height = 12,width = 15,units = 'in',res=300)
    ht = Heatmap(t(log2(test))[rownames(annotation_row),],
                # col = colorRamp2(c(-0.6,0,0.6), c("#253b70","white","#7c1520")),
                col = col_fun,
                cluster_rows = F,
                cluster_columns = F,
                show_column_names = F,
                show_row_names = F,
                column_title = i,
                column_title_gp = gpar(fontsize = 14,color="black"),
                height = unit(8, "cm"),
                width = unit(18, "cm"),
                column_split = factor(gene_pos$V2, new_cluster),
                heatmap_legend_param = list(title = "inferCNV",
                                            direction = "horizontal",
                                            title_position = "lefttop",
                                            legend_height = unit(3, "cm")),
                
                left_annotation = left_anno, 
                use_raster = FALSE,
                # row_split = factor(annotation_row$celltype, levels = c('Reference','Epithelial','Tumour')),
                row_title = NULL,
                top_annotation = top_color,
                border = T)
    draw(ht, heatmap_legend_side = "top")
    dev.off()
    htlist[[i]] <-  draw(ht, heatmap_legend_side = "top")
}

### 合并在一起进行绘制 ############
# 先对数据列表中的每一个样本进行分析维度，查看出最低维度的样本，以这个维度为基础进行分析
# 对列表中的没一个样本进行筛选
test <- list()
for (i in names(flist)) {
    print(i)
    test[[i]] <- cbind(flist3[[i]], flist[[i]], flist2[[i]])
}
## 需要筛选出所有样本共有的基因型
genes <- lapply(unique(names(test)), function(x) rownames(test[[x]]))
common.genes <- Reduce("intersect", genes)

test.sub <- list()
for (i in names(test)) {
    test.sub[[i]] <- test[[i]][common.genes, ]
}

## 对上面的flist,flist2,flist3文件分别进行筛选
### flist normal tumor
for (i in names(flist)) {
    flist[[i]] <- flist[[i]][common.genes, ]
}
### flist2 tumor
for (i in names(flist2)) {
    flist2[[i]] <- flist2[[i]][common.genes, ]
}
### flist3 Reference
for (i in names(flist3)) {
    flist3[[i]] <- flist3[[i]][common.genes, ]
}

flist2.2 <- list()
for (i in names(flist2)) {
    print(i)
    if (dim(flist2[[i]])[2] > 500) {
        flist2.2[[i]] <- flist2[[i]][, names(sort(colSums(flist2[[i]]), decreasing = T)[1:500])]
    } else {
        flist2.2[[i]] <- flist2[[i]]
    }
}


flist3.2 <- list()
for (i in names(flist3)) {
    print(i)
    if (dim(flist3[[i]])[2] > 500) {
        flist3.2[[i]] <- flist3[[i]][, names(sort(colSums(flist3[[i]]), decreasing = F)[1:500])]
    } else {
        flist3.2[[i]] <- flist3[[i]]
    }
}

#### 合并数据 ####################
test <- list()
for (i in names(flist)) {
    print(i)
    test[[i]] <- cbind(flist3.2[[i]], flist[[i]], flist2.2[[i]])
}
tmp <- Reduce('cbind',test)
dim(tmp) # [1]  4359 25764
tmp <- t(log2(tmp))

## 生成排序数据
df1 <- list()
df2 <- list()
for (i in names(flist3.2)) {
    print(i)
    cell_anno <- data.frame(cell_id = colnames(cbind(flist3.2[[i]], flist[[i]], flist2.2[[i]])), 
        group = c(rep("Reference", dim(flist3.2[[i]])[2]), rep("Normal epithelial cells", dim(flist[[i]])[2]), rep("Malignant epithelial cells", dim(flist2.2[[i]])[2])))
    cell_anno$sample <- i
    df1[[i]] <- cell_anno
    annotation_row <- data.frame(celltype = c(rep("Reference", dim(flist3.2[[i]])[2]), rep("Normal epithelial cells", dim(flist[[i]])[2]), rep("Malignant epithelial cells", dim(flist2.2[[i]])[2])))
    annotation_row$sample <- i
    rownames(annotation_row) <- colnames(cbind(flist3.2[[i]], flist[[i]], flist2.2[[i]]))
    df2[[i]] <- annotation_row
}
# 合并出数据框
top_anno <- Reduce('rbind',df1)
annotation_row <- Reduce('rbind',df2)


gene_pos <- read.delim("/root/wangje/software/InferCNV/hg38_gencode_v27.txt", header = F) # 基因坐标信息
gene_pos <- gene_pos[gene_pos$V1 %in% colnames(tmp),]
new_cluster <- unique(gene_pos$V2)
column_split = factor(gene_pos$V2, new_cluster)


sample_col = cols[1:length(names(flist))]
names(sample_col) <- names(flist)
left_anno <- rowAnnotation(df = annotation_row,
                        col=list(sample = sample_col,ncol=3,celltype=c("Reference"="#2166ac","Normal epithelial cells" = "#b2182b","Malignant epithelial cells"="#ce9d3f"), title_position = "topcenter"),
                        show_annotation_name = F)

top_color <- HeatmapAnnotation(
    cluster = anno_block(
        labels = gsub("chr", "", new_cluster),
        gp = gpar(col ="white"),
        labels_gp = gpar(cex = 1, col = "black"),
        labels_rot = 45,
        height = unit(5,"mm"))) 


ht <- Heatmap(tmp,
    # col = colorRamp2(c(-0.6,0,0.6), c("#253b70","white","#7c1520")),
    col = col_fun,
    cluster_rows = F,
    cluster_columns = F,
    show_column_names = F,
    show_row_names = F,
    column_title = NULL,
    column_title_gp = gpar(fontsize = 14, color = "black"),
    height = unit(14, "cm"),
    width = unit(18, "cm"),
    column_split = factor(gene_pos$V2, new_cluster),
    heatmap_legend_param = list(
        title = "inferCNV",
        direction = "horizontal",
        title_position = "lefttop",
        legend_height = unit(3, "cm")
    ),
    left_annotation = left_anno,
    use_raster = FALSE,
    # row_split = factor(annotation_row$celltype, levels = c('Reference','Epithelial','Tumour')),
    row_title = NULL,
    top_annotation = top_color,
    border = T
)

##################################################################################################################################################################################
tt <- factor(unique(names(flist)),levels = c(
        "P1-Pre", "P1-Post", "P2-Pre", "P2-Post", "P3-Pre", "P3-Post", "P4-Pre", "P4-Post", "P5-Pre", "P5-Post", "P6-Pre", "P6-Post",
        "P7-Pre", "P7-Post", "P8-Pre", "P8-Post", "P9-Pre","P10-Pre", "P10-Post", "P11-Pre", "P11-Post", "P12-Pre", "P12-Post", "P13-Pre", "P13-Post",
        "P14-Pre", "P14-Post", "P15-Pre", "P16-Pre", "P17-Pre"))
tt <- levels(tt)

ilist1 <- list()
for(i in tt){
    ilist1[[i]] <- flist[[i]]
}

ilist2 <- list()
for(i in tt){
    ilist2[[i]] <- flist2.2[[i]]
}

ilist3 <- list()
for(i in tt){
    ilist3[[i]] <- flist3.2[[i]]
}


test <- list()
for (i in tt) {
    print(i)
    test[[i]] <- cbind(ilist3[[i]], ilist2[[i]], ilist2[[i]])
}
###############################################
tmp <- Reduce('cbind',test)
dim(tmp) # [1]  4359 25764
tmp <- t(log2(tmp))


df1 <- list()
df2 <- list()
for (i in names(ilist1)) {
    print(i)
    cell_anno <- data.frame(cell_id = colnames(cbind(ilist3[[i]], ilist1[[i]], ilist2[[i]])), 
        group = c(rep("Reference", dim(ilist3[[i]])[2]), rep("Normal epithelial cells", dim(ilist1[[i]])[2]), rep("Malignant epithelial cells", dim(ilist2[[i]])[2])))
    cell_anno$sample <- i
    df1[[i]] <- cell_anno
    annotation_row <- data.frame(celltype = c(rep("Reference", dim(ilist3[[i]])[2]), rep("Normal epithelial cells", dim(ilist1[[i]])[2]), rep("Malignant epithelial cells", dim(ilist2[[i]])[2])))
    annotation_row$sample <- i
    rownames(annotation_row) <- colnames(cbind(ilist3[[i]], ilist1[[i]], ilist2[[i]]))
    df2[[i]] <- annotation_row
}
# 合并出数据框
top_anno <- Reduce('rbind',df1)
annotation_row <- Reduce('rbind',df2)

gene_pos <- read.delim("/root/wangje/software/InferCNV/hg38_gencode_v27.txt", header = F) # 基因坐标信息
gene_pos <- gene_pos[gene_pos$V1 %in% colnames(tmp),]
new_cluster <- unique(gene_pos$V2)
column_split = factor(gene_pos$V2, new_cluster)


sample_col = cols[1:length(tt)]
names(sample_col) <- tt
left_anno <- rowAnnotation(df = annotation_row,
                        col=list(sample = sample_col,ncol=3,celltype=c("Reference"="#2166ac","Normal epithelial cells" = "#b2182b","Malignant epithelial cells"="#ce9d3f"), title_position = "topcenter"),
                        show_annotation_name = F)

top_color <- HeatmapAnnotation(
    cluster = anno_block(
        labels = gsub("chr", "", new_cluster),
        gp = gpar(col ="white"),
        labels_gp = gpar(cex = 1, col = "black"),
        labels_rot = 45,
        height = unit(5,"mm"))) 









cols <- c("#cb7c77", "#68d359", "#6a7dc9", "#c9d73d", "#c555cb",
"#d7652d", "#7cd5c8", "#c49a3f", "#507d41", "#5d8d9c",
"#90353b", "#674c2a", "#1B9E77", "#c5383c", "#0081d1",
"#ffd900", "#502e71", "#c8b693", "#aed688", "#f6a97a",
"#c6a5cc", "#798234", "#6b42c8", "#cf4c8b", "#666666",
"#feb308", "#ff1a1a", "#1aff1a", "#1a1aff", "#ffff1a",
"#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
"#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
"#437BFE", "#FEC643", "#43FE69", "#FE6943", "#C643FE",
"#43D9FE", "#B87A3D", "#679966", "#993333", "#7F6699",
"#E78AC3", "#333399", "#A6D854", "#E5C494",
"#87b3d4", "#d5492f", "#6bd155", "#683ec2", "#c9d754",
"#d04dc7", "#81d8ae", "#d34a76", "#607d3a", "#6d76cb",
"#ce8d3f", "#81357a", "#d3c3a4", "#3c2f5a", "#b96f49",
"#4e857e", "#6e282c", "#d293c8", "#393a2a", "#997579",
"#87b3d4", "#d5492f", "#6bd155", "#683ec2", "#c9d754",
"#d04dc7", "#81d8ae", "#d34a76", "#607d3a", "#6d76cb",
"#ce9d3f", "#81357a", "#d3c3a4", "#3c2f5a", "#b96f49",
"#4e857e", "#6e282c", "#d293c8", "#393a2a", "#997579",
"#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
"#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
"#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
"#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
"#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
"#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
"#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
"#e598d7", "#343b20", "#bbb5d9", "#975223", "#576e8b",
"#d97f5e", "#253e44", "#de959b", "#417265", "#712b5b",
"#8c6d30", "#a56c95", "#5f3121", "#8f846e", "#8f5b5c")



htlist <- list()
for(i in names(flist)[1:5]){
    print(i)
    test <- cbind(flist3[[i]],flist[[i]],flist2[[i]]) # 合并文件
    print(dim(test))
    gene_pos <- read.delim("/root/wangje/software/InferCNV/hg38_gencode_v27.txt", header = F) # 基因坐标信息
    cell_anno =data.frame(cell_id =colnames(test),
                        group =c(rep("Reference",dim(flist3[[i]])[2]),rep("Normal epithelial cells",dim(flist[[i]])[2]),rep('Malignant epithelial cells',dim(flist2[[i]])[2]))) # 细胞注释信息
    gene_pos <- gene_pos[gene_pos$V1 %in% rownames(test),]
    new_cluster <- unique(gene_pos$V2)
    top_color <- HeatmapAnnotation(
        cluster = anno_block(
            labels = gsub("chr", "", new_cluster),
            gp = gpar(col ="white"),
            labels_gp = gpar(cex = 1, col = "black"),
            labels_rot = 45,
            height = unit(5,"mm"))) 

    annotation_row = data.frame(sample = rep(i,c(dim(flist[[i]])[2]+dim(flist2[[i]])[2]+dim(flist3[[i]])[2])),
                                celltype = c(rep("Reference",dim(flist3[[i]])[2]),
                                        rep("Normal epithelial cells",dim(flist[[i]])[2]),
                                        rep('Malignant epithelial cells',dim(flist2[[i]])[2])))
    ## 设置样本颜色
    sample_col = cols[1:length(names(flist))]
    names(sample_col) <- names(flist)
    rownames(annotation_row) <- colnames(test)
    left_anno <- rowAnnotation(df = annotation_row,
                            col=list(sample= sample_col[i],ncol=3,celltype=c("Reference"="#2166ac","Normal epithelial cells" = "#b2182b","Malignant epithelial cells"="#ce9d3f"), title_position = "topcenter"),
                            show_annotation_name = F)
    png(paste0('003_',i,"_CNV_heatmap.png"),height = 12,width = 15,units = 'in',res=300)
        ht = Heatmap(t(log2(test))[rownames(annotation_row),],
                    # col = colorRamp2(c(-0.6,0,0.6), c("#253b70","white","#7c1520")),
                    col = col_fun,
                    cluster_rows = F,
                    cluster_columns = F,
                    show_column_names = F,
                    show_row_names = F,
                    # column_title = i,
                    column_title_gp = gpar(fontsize = 14,color="black"),
                    height = unit(3, "cm"),
                    width = unit(18, "cm"),
                    column_split = factor(gene_pos$V2, new_cluster),
                    heatmap_legend_param = list(title = "inferCNV",
                                                direction = "horizontal",
                                                title_position = "lefttop",
                                                legend_height = unit(3, "cm")),
                    
                    left_annotation = left_anno, 
                    use_raster = FALSE,
                    # row_split = factor(annotation_row$celltype, levels = c('Reference','Epithelial','Tumour')),
                    row_title = NULL,
                    top_annotation = top_color,
                    border = T)
        draw(ht, heatmap_legend_side = "top")
        dev.off()
        htlist[[i]] <-  ht
    }

##### 将全部样本的infercnv数据合并成一个文件进行绘图
# 合并列表中的文件为一个大文件
tmp <- Reduce('cbind',test.sub) #








## 分开绘制
annotation_row <- data.frame(celltype = c(rep("Ref",dim(flist3[[i]])[2])))
left_anno <- rowAnnotation(df = annotation_row,col = list(celltype = c("Ref"="#2166ac")))
ht1 <- Heatmap(t(log2(flist3[[1]])),
                col = colorRamp2(c(-0.5,0,0.5), c("#2166ac","white","#b2182b")),
                cluster_rows = F,
                cluster_columns = F,
                height = unit(2,'cm'),
                width = unit(26,'cm'),
                show_column_names = F,
                border = TRUE,
                left_annotation = left_anno,
                border_gp = gpar(color = "black",lwd = unit(1,'cm')),
                show_row_names = F)

ht2 <- Heatmap(t(log2(flist2[[1]])),
                col = colorRamp2(c(-0.5,0,0.5), c("#2166ac","white","#b2182b")),
                cluster_rows = F,
                cluster_columns = F,
                height = unit(5,'cm'),
                width = unit(26,'cm'),
                show_column_names = F,
                show_row_names = F)

## 合在一起注释

htList <- list()
for(i in names(flist)) {
    print(i)
    
    test1 <- cbind(flist3[[i]], flist[[i]])  # 合并文件
    test2 <- flist2[[i]]                     # 合并文件
    
    print(dim(test1))
    print(dim(test2))
    
    gene_pos <- read.delim("/root/wangje/software/InferCNV/hg38_gencode_v27.txt", header = FALSE) # 基因坐标信息
    
    cell_anno1 <- data.frame(
        cell_id = colnames(test1),
        group = c(rep("Ref", dim(flist3[[i]])[2]), rep("Epithelial", dim(flist[[i]])[2]))) # 细胞注释信息
    cell_anno2 <- data.frame(
        cell_id = colnames(test2),
        group = c(rep("Tumour", dim(flist2[[i]])[2]))) # 细胞注释信息
    
    # Ensure `test` is correctly referenced; assuming it should be `test1` or `test2`
    gene_pos <- gene_pos[gene_pos$V1 %in% rownames(test1), ]
    new_cluster <- unique(gene_pos$V2)
    
    top_color <- HeatmapAnnotation(
        cluster = anno_block(
            labels = gsub("chr", "", new_cluster),
            gp = gpar(col = "white"),
            labels_gp = gpar(cex = 1, col = "black"),
            height = unit(5, "mm")))

    annotation_row1 <- data.frame(
        celltype = c(rep("Reference", dim(flist3[[i]])[2]), rep("Epithelial", dim(flist[[i]])[2]))
    )
    rownames(annotation_row1) <- colnames(test1)
    
    annotation_row2 <- data.frame(
        celltype = c(rep("Tumour", dim(flist2[[i]])[2]))
    )
    rownames(annotation_row2) <- colnames(test2)
    
    # Remove unnecessary filtering since all rows already have the desired cell types
    # annotation_row1 <- annotation_row1[annotation_row1$celltype %in% c("Reference", "Epithelial"), ]

    left_anno1 <- rowAnnotation(
        df = annotation_row1,
        col = list(celltype = c("Reference" = "#347fb9", "Epithelial" = "#e51316")),
        show_annotation_name = FALSE
    )


    ht1 <- Heatmap(
        t(log2(test1))[rownames(annotation_row1), ],
        col = colorRamp2(c(-0.5, 0, 0.5), c("#2166ac", "white", "#b2182b")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        column_title = i,
        name = "CNV score",
        show_row_names = FALSE,
        column_split = factor(gene_pos$V2, levels = new_cluster),
        column_gap = unit(0, "mm"),
        heatmap_legend_param = list(
            title = "inferCNV",
            direction = "vertical",
            title_position = "topcenter",
            legend_height = unit(3, "cm")
        ),
        height = unit(3, "cm"),
        width = unit(18, "cm"),
        left_annotation = left_anno1,
        row_title = NULL,
        top_annotation = top_color,
        border = TRUE
    )
    
    # Remove unnecessary filtering since all rows already have the desired cell type
    # annotation_row2 <- annotation_row2[annotation_row2$celltype %in% c("Tumour"), ]

    left_anno2 <- rowAnnotation(
        df = annotation_row2,
        col = list(celltype = c("Tumour" = "#4cb049")),
        show_annotation_name = FALSE
    )
    
    ht2 <- Heatmap(
        t(log2(test2))[rownames(annotation_row2), ],
        col = colorRamp2(c(-0.5, 0, 0.5), c("#2166ac", "white", "#b2182b")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        column_title = i,
        name = "CNV score",
        show_row_names = FALSE,
        column_split = factor(gene_pos$V2, levels = new_cluster),
        column_names_gp = gpar(fontsize = 0), 
        column_gap = unit(0, "mm"),
        heatmap_legend_param = list(
            title = "inferCNV",
            direction = "vertical",
            title_position = "topcenter",
            legend_height = unit(3, "cm")
        ),
        height = unit(5.5, "cm"),
        width = unit(18, "cm"),
        left_annotation = left_anno2,
        row_title = NULL,
        top_annotation = top_color,
        border = TRUE
    )
    
    png(
        paste0('002', i, "_CNV_heatmap.png"),
        height = 8,
        width = 18,
        units = 'in',
        res = 300
    )
    
    draw(ht1 %v% ht2, ht_gap = unit(0.2, "cm"))
    dev.off()
    
    htList[[i]] <- draw(ht1 %v% ht2, ht_gap = unit(0.2, "cm"))
}



## --------------------------------------------------------------------------------------------------------------------------------------------------
## 角质细胞和肝细胞
library(Seurat)
library(dplyr)
library(patchwork)
library(infercnv)
library(ComplexHeatmap)
library(data.table)
setwd('/root/wangje/Project/刘老师/大群/Epithelials/otherEpi')

# 读入infercnv的结果
infercnv_obj <- readRDS('run.final.infercnv_obj')
str(infercnv_obj)
class(infercnv_obj)
# [1] "infercnv"
# attr(,"package")
# [1] "infercnv"
slotNames(infercnv_obj)
# [1] "expr.data"                        "count.data"                       "gene_order"                       "reference_grouped_cell_indices"
# [5] "observation_grouped_cell_indices" "tumor_subclusters"                "options"                          ".hspike"

# 提取参考细胞和测试细胞
ref.cells.endo <- colnames(infercnv_obj@expr.data)[infercnv_obj@reference_grouped_cell_indices$`Endothelial cells`]
ref.cells.fib <- colnames(infercnv_obj@expr.data)[infercnv_obj@reference_grouped_cell_indices$Fibroblasts]
test.cells.epi1 <- colnames(infercnv_obj@expr.data)[infercnv_obj@observation_grouped_cell_indices$`other Epi 1`]
test.cells.epi2 <- colnames(infercnv_obj@expr.data)[infercnv_obj@observation_grouped_cell_indices$`other Epi 2`]

# cell name
cells <- data.frame(
    cells = c(ref.cells.endo,ref.cells.fib,test.cells.epi1,test.cells.epi2),
    group = c(rep("Endothelial",length(ref.cells.endo)),rep("Fibroblasts",length(ref.cells.fib)),
            rep("otherEpi1",length(test.cells.epi1)),rep("otherEpi2",length(test.cells.epi2)))
)
head(cells)
#                                cells       group
# 1 CJME_0707_AACGCTTAATAGCGACAAACATCG Endothelial
# 2 CJME_0707_AACGTGATACATTGGCAGCACCTC Endothelial
# 3 CJME_0707_AACGTGATATAGCGACATGCCTAA Endothelial
# 4 CJME_0707_AACGTGATCGGATTGCCACCTTAC Endothelial
# 5 CJME_0707_AAGACGGAAAGAGATCACACGACC Endothelial
# 6 CJME_0707_AAGACGGACGCATACATATCAGCA Endothelial

# 添加CNV score
expr <- infercnv_obj@expr.data
CNV_score=as.data.frame(colMeans((expr-1)^2)) 
colnames(CNV_score) <- "CNV_score"
CNV_score <- CNV_score %>% tibble::rownames_to_column(var = "cells")

CNV_score <- dplyr::left_join(CNV_score,cells, by = "cells")
head(CNV_score)
#                                cells    CNV_score       group
# 1 CJME_0707_AACGCTTAATAGCGACAAACATCG 0.0013940044 Endothelial
# 2 CJME_0707_AACGTGATACATTGGCAGCACCTC 0.0005838440 Endothelial
# 3 CJME_0707_AACGTGATATAGCGACATGCCTAA 0.0022108258 Endothelial
# 4 CJME_0707_AACGTGATCGGATTGCCACCTTAC 0.0011338771 Endothelial
# 5 CJME_0707_AAGACGGAAAGAGATCACACGACC 0.0009057562 Endothelial
# 6 CJME_0707_AAGACGGACGCATACATATCAGCA 0.0015402194 Endothelial

# 绘制小提琴图
ggplot(CNV_score, aes(group,CNV_score))+
    geom_violin(aes(fill=group),color="black")+
    scale_fill_manual(values = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#E78AC3", "#333399", "#A6D854", "#E5C494"))+
    theme_bw()+
    ggtitle("Other Epithelial cells")+
    theme(axis.text.x = element_text(colour = "black", size = 12),
            axis.text.y = element_text(colour = "black", size = 10),
            axis.title = element_text(color = 'black', size = 12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
    ggpubr::stat_compare_means()+
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), 
                geom = "pointrange", color = "black", size=0.5)+
    stat_summary(fun.y = "mean", geom = "line", aes(group = 1),color="orange",linewidth = 1.1) # 添加趋势线


gg1 <- ggplot(data=CEC, aes(x = variable, y = value,fill=variable)) +
  geom_boxplot(alpha =0.5,size=0.5,outlier.shape = NA)+
  scale_fill_manual(values=c("#E29827","#922927"))+
  stat_compare_means(method = "t.test",paired = F, 
                     comparisons=list(c("HC_CEC", "AEH_CEC")))+
  geom_jitter(alpha = 0.3,size=3, shape=21)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  theme_bw() + 
  theme(panel.grid =element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')+
  labs(y='Pathway activity score')



## 添加其他上皮细胞数据
setwd('/root/wangje/Project/刘老师/大群/Epithelials/Infercnv/')
flist <- list()
for(sample in list.files(path = '.', pattern = "^R|^NR",full.names = T)){
    sp <- stringr::str_split_fixed(sample,pattern = "\\./",n=2)[,2]
    message(sp)
    # 读入数据
    cat("数据地址: ",file.path('/root/wangje/Project/刘老师/大群/Epithelials/Infercnv',sample,"002_infercnv_obj.qs"))
    infercnv_obj <- qs::qread(file.path('/root/wangje/Project/刘老师/大群/Epithelials/Infercnv',sample,"002_infercnv_obj.qs"))
    expr <- infercnv_obj@expr.data
    CNV_score=as.data.frame(colMeans((expr-1)^2))
    colnames(CNV_score) <- "CNV_score"
    CNV_score <- CNV_score %>% tibble::rownames_to_column(var = 'cells')
    flist[[sp]] <- CNV_score
}
AllCells <- do.call("rbind",flist)

## 读入预测的上皮细胞
tmuor <- read.table('/root/wangje/Project/刘老师/大群/Epithelials/Infercnv/result/01_预测肿瘤上皮细胞.txt')
epi <- read.table('/root/wangje/Project/刘老师/大群/Epithelials/Infercnv/result/01_预测正常上皮细胞.txt')
# 区分出预测的正常上皮和肿瘤细胞
AllCells <- AllCells[AllCells$cells %in% c(tmuor$tumour,epi$epi),]
AllCells <- AllCells %>% dplyr::mutate(
    celltype = case_when(
        cells %in% tmuor$tumour ~ 'Malignant epithelial cells',
        cells %in% epi$epi ~ "Epithelial cells"))
 head(AllCells)
#                                   cells   CNV_score                   celltype
# 4789 CZYI_0702_AAACATCGAACGCTTAAGCACCTC 0.007591962 Malignant epithelial cells
# 4790 CZYI_0702_AAACATCGAACTCACCAACGCTTA 0.006146286 Malignant epithelial cells
# 4791 CZYI_0702_AAACATCGAAGACGGAGGTGCGAA 0.006143767 Malignant epithelial cells
# 4792 CZYI_0702_AAACATCGAAGAGATCACAGCAGA 0.007797905 Malignant epithelial cells
# 4793 CZYI_0702_AAACATCGAAGAGATCCCGACAAC 0.006631365 Malignant epithelial cells
# 4794 CZYI_0702_AAACATCGAAGAGATCCGACTGGA 0.005786138 Malignant epithelial cells

colnames(CNV_score)[3] <- 'celltype'
df <- rbind(CNV_score,AllCells)

ggplot(df,aes(celltype,CNV_score))+
  geom_jitter(alpha = 0.1,size = 0.01)+
    geom_violin(aes(fill=celltype),color="black",alpha = 0.9)+
    scale_fill_manual(values = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#E78AC3", "#333399", "#A6D854", "#E5C494"))+
    theme_void()+
    # ggtitle("Other Epithelial cells")+
    theme(axis.text.x = element_text(colour = "black", size = 15,angle = 45, hjust = 1,vjust = 1),
            axis.text.y = element_text(colour = "black", size = 15),
            axis.title.y = element_text(color = 'black', size = 17,angle = 90),
            legend.text = element_text(size = 11),
            legend.title = element_blank(),
            panel.border = element_rect(linewidth = 1.3,color = "black"),
            axis.title.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
    # ggpubr::stat_compare_means()+
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), 
                geom = "pointrange", color = "black", size=0.1)+
    stat_summary(fun.y = "mean", geom = "line", aes(group = 1),color="orange",linewidth = 0.5) # 添加趋势线

## 读取每个样本的交集进行分析
# 使用数据为上面读取的flist,flist2,flist3




## ---------------------------------------------------------------------------------------------------------------------------------------------------
## 使用相关性和阈值去找恶性细胞


#读入inferCNV运行完成后的两个数据，其实本质是reference cell和观测细胞CNV表达矩阵
obs <- read.table("./infercnv.observations.txt", header=T, check.names=F)
ref <- read.table("./infercnv.references.txt", header=T, check.names=F)
#这里引用一个函数，但是本质和上面计算CNV_score是一样的，只不过将reference和观测细胞分开做了而已
#此外还计算了相关性

estimateCNV <- function(obs,#读入的obs矩阵
                        ref, #读入的ref矩阵
                        score_threshold, #CNV_score阈值，大于多少定义为Malignant cell
                        #个人人为这个参数阈值有很大的主观性
                        cor_threshold#相关性阈值，也是一个主观性参数
                        ){
  #obs和ref的CNV_score计算
  cnv_obs <- colMeans((obs-1)^2)
  cnv_ref <- colMeans((ref-1)^2) 
  #挑选具有top cnv的cell，我这里挑选观测细胞前50%的细胞CNV,并计算每个细胞的CNV平均值
  cell_top <- names(sort(cnv_obs, decreasing=T))[1:round(length(cnv_obs)*0.5)]
  cnv_top <- rowMeans(obs[, cell_top])
  #分别计算top cnv与obs/ref CNV的相关性
  cor_obs <- apply(obs, 2, function(x)cor(x, cnv_top))
  cor_ref <- apply(ref, 2, function(x)cor(x, cnv_top))
  
  cnv <- data.frame(score=c(cnv_obs, cnv_ref), 
                    cor=c(cor_obs, cor_ref), 
                    barcode=c(colnames(obs), colnames(ref)))
  
  #设定阈值，鉴定Malignant cell
  
  cnv$type <- 'Other'
  cnv$type[cnv$score>score_threshold & cnv$cor>cor_threshold] <- 'Malignant'
  cnv$type[cnv$score<score_threshold & cnv$cor<cor_threshold] <- 'Not Malignant'
  return(cnv)
}


CNV_score_cor <- estimateCNV(obs, ref, 
                             score_threshold= 0.003, 
                             cor_threshold= 0.3)
colnames(CNV_score_cor)[3] <- "cell_id"
#与之前的分组数据结合
CNV_score_cor=CNV_score_cor%>%inner_join(kmeans_df,by="cell_id")
#因为我们需要判断EEC肿瘤中的Malignant cell，所以我们可以只提取EEC组的数据作图

CNV_score_cor_EEC <- subset(CNV_score_cor, group=='EEC')

#作图
library(ggplot2)
ggplot(CNV_score_cor_EEC, aes(x=score,y=cor,color=type))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title = element_text(color = 'black', size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="CNV Score",
       y="CNV Correlation")+
  scale_color_manual(values =c("#009E73","grey","black") )


library(ggplot2)
library(qs)
library(infercnv)
library(dplyr)
library(patchwork)
library(ggpubr)

plotList <- list()
for(sample in list.files(path = '.', pattern = "^R|^NR",full.names = T)){
    sp <- stringr::str_split_fixed(sample,pattern = "\\./",n=2)[,2]
    message(sp)
    # 读入数据
    cat("数据地址: ",file.path('/root/wangje/Project/刘老师/大群/Epithelials/Infercnv',sample,"002_infercnv_obj.qs"))
    infercnv_obj <- qs::qread(file.path('/root/wangje/Project/刘老师/大群/Epithelials/Infercnv',sample,"002_infercnv_obj.qs"))
    expr <- infercnv_obj@expr.data
    #提取reference cell位置
    ref_cell <- c(infercnv_obj@reference_grouped_cell_indices$`Endothelial cells`,
                infercnv_obj@reference_grouped_cell_indices$`Fibroblasts`)
    #提取测试细胞的位置
    obs_cell <- c(infercnv_obj@observation_grouped_cell_indices$`Epithelail cells`)
    # 细胞注释信息
    cell_anno =data.frame(cell_id =c(colnames(expr)[ref_cell],colnames(expr)[obs_cell]),
                      group =c(rep("Endothelial cells",length(infercnv_obj@reference_grouped_cell_indices$`Endothelial cells`)),
                      rep("Fibroblasts",length(infercnv_obj@reference_grouped_cell_indices$`Fibroblasts`)),rep("Epithelial",length(obs_cell))))
    # 读入需要的数据
    obs <- read.table(file.path('/root/wangje/Project/刘老师/大群/Epithelials/Infercnv',sample,"otherEpi","infercnv.observations.txt"), header=T, check.names=F)
    ref <- read.table(file.path('/root/wangje/Project/刘老师/大群/Epithelials/Infercnv',sample,"otherEpi","infercnv.references.txt"), header=T, check.names=F)
    CNV_score_cor <- estimateCNV(obs, ref, 
                             score_threshold= 0.003, 
                             cor_threshold= 0.3)
    CNV_score_cor <- CNV_score_cor %>% tibble::rownames_to_column(var = "cell_id")
    # 合并细胞注释结果和绘图数据
    CNV_score_cor <- CNV_score_cor %>% dplyr::left_join(cell_anno, by = "cell_id")

    #作图
    p <- ggplot(CNV_score_cor, aes(x=score,y=cor,color=group))+
        geom_point()+
        theme_bw()+
        ggtitle(sp)+
        theme(axis.text.x = element_text(colour = "black", size = 12),
                axis.text.y = element_text(colour = "black", size = 10),
                axis.title = element_text(color = 'black', size = 12),
                plot.title = element_text(size = 14,color = 'black'),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
        labs(x="CNV Score",y="CNV Correlation")+
        scale_color_manual(values =c("#d7652d", "#7cd5c8", "#c49a3f", "#507d41"))
    plotList[[sp]] <- p
    }

## 保存绘图的结果
png("./CNV相关性结果01.png", height = 20,width = 32, units = 'in', res=300)
ggpubr::ggarrange(plotlist = plotList, ncol = 6, nrow = 5)
dev.off()


#--------------------------------------------------------------------------------------------------------------------------
prefix <- "scRNA"
reduction="SeuratUMAP2D"
scRNA$sample <- stringr::str_split_fixed(scRNA$new_Rename,'-',n=2)[,1]
scRNA$sample <- factor(scRNA$sample, levels = paste0('P',1:17))
scRNA$Treatment <- scRNA$Treat_assess
p1 <- Seurat::DimPlot(scRNA, group.by = "Seuratpca_SNN_res.0.4", label = T, repel = T,reduction = reduction) + SCP::theme_blank(xlab='UMAP1',ylab='UMAP2')
p2 <- Seurat::DimPlot(scRNA, group.by = c('sample','Treatment',"Tissue"),ncol = 3,label = F, reduction = reduction) & SCP::theme_blank(xlab='UMAP1',ylab='UMAP2')

# 保存图片
design <- "ABBB"
png(paste0("Seurat","_DimPlot.png"), height = 4, width = 18,units = 'in', res = 300)
p1 + p2 + plot_layout(design = design)+ plot_annotation(tag_levels = "A")
dev.off()


plotFeaturePlot <- function(srt, marks, reduction){
  plist <- list()
    for(i in names(marks)){
      for(j in marks[[i]]){
        res <- tryCatch({
          FeaturePlot(srt, features = j, raster = TRUE, reduction = reduction) +
          SCP::theme_blank(xlab='UMAP1',ylab='UMAP2')+
          labs(title = paste0(i, '_', j))
        }, error = function(e){
          message(paste0(i, '_', j, ' Error...'))
          message(e$message)
          return(NULL)
        }, finally = {
          print(paste0(i, '_', j, ' Finish...'))
        })
        plist[[paste0(i, '_', j)]] <- res
        }
      }
      return(plist)
}



# FeaturePlot
markers <- list(
  Epithelial = c("CD24", "KRT19", "EPCAM",'KRT18'),
  Cycling = c('MKI67','TOP2A'),
  HER2 = c('ERBB2')
  )
png(paste0("Seurat","_FeaturePlot.png"), height = 8, width = 25,units = 'in', res = 300)
patchwork::wrap_plots(plotFeaturePlot(scRNA,marks = markers, reduction = "SeuratUMAP2D"),ncol = 4)
dev.off()



png('./test.png', height = 35, width = 60, units = 'in', res=300)
Reduce('%v%', htlist)
dev.off()
