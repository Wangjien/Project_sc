library(Seurat)
library(qs)
library(patchwork)
library(dplyr)
library(ggplotify)
library(ggplot2)
library(data.table)
library(tidyfast)  
library(stringr) 

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
qsavem(meta.nkt,counts.nkt,file = '/root/wangje/Project/刘老师/NK_T/Data/nkt_cellphonedb.qs')

#<<<<<<<<< 读入Bcells >>>>>>>>>>>>>>>>>>>>>>>
setwd("/root/wangje/Project/刘老师/Bcells/Data")
load('new_Bcells_CCA.RData')
# 新的celltype
scRNA_seurat$cellphonedb_celltype = ifelse(
    scRNA_seurat$celltype %in% c('Plasma cells_1','Plasma cells_2') ,'Plasma cells', 'B cells'
)
# 提取表达矩阵
counts.bcells = as.data.frame(scRNA_seurat@assays$RNA@counts)
# 提取celltype
Idents(scRNA_seurat) = scRNA_seurat$cellphonedb_celltype
meta.bcells=as.data.frame(scRNA_seurat@active.ident)
meta.bcells$cell=row.names(meta.bcells)
colnames(meta.bcells) = c('cellphonedb_celltype','cell')
meta.bcells=meta.bcells[,c("cell","cellphonedb_celltype")]
qsavem(meta.bcells,counts.bcells,file = 'bcells_cellphonedb.qs')

#<<<<<<<<<<<<<< 读入Endothelials cells >>>>>>>>>>>>>>>>>
setwd('"/root/wangje/Project/刘老师/Endothelials"')
load('Liu_Endothelials_seurat.RData')
scRNA_seurat$cellphonedb_celltype = 'Endothelial cells'

# 提取表达矩阵
counts.endo = as.data.frame(scRNA_seurat@assays$RNA@counts)
# 提取celltype
Idents(scRNA_seurat) = scRNA_seurat$cellphonedb_celltype
meta.endo=as.data.frame(scRNA_seurat@active.ident)
meta.endo$cell=row.names(meta.endo)
colnames(meta.endo) = c('cellphonedb_celltype','cell')
meta.endo=meta.endo[,c("cell","cellphonedb_celltype")]
qsavem(meta.endo,counts.endo,file = 'endo_cellphonedb.qs')

#<<<<<<<<<<<<<<<<< 读入Fibroblasts >>>>>>>>>>>>>>>>>>>>>
setwd('/root/wangje/Project/刘老师/Fibroblasts')
load('fibroblasts_cca.RData')
scRNA_seurat$cellphonedb_celltype = 'Fibroblasts'
# 提取表达矩阵
counts.fib = as.data.frame(scRNA_seurat@assays$RNA@counts)
# 提取celltype
Idents(scRNA_seurat) = scRNA_seurat$cellphonedb_celltype
meta.fib=as.data.frame(scRNA_seurat@active.ident)
meta.fib$cell=row.names(meta.fib)
colnames(meta.fib) = c('cellphonedb_celltype','cell')
meta.fib=meta.fib[,c("cell","cellphonedb_celltype")]
qsavem(meta.fib,counts.fib,file = 'fib_cellphonedb.qs')

#<<<<<<<<<<<<<<<<<< 读入Myeloids cells >>>>>>>>>>>>>>>>>>>>
setwd('/root/wangje/Project/刘老师/Myeloids/Data')
load('Myeloids_CCA.RData')
# 添加新的celltype
scRNA_seurat$cellphonedb_celltype = ifelse(
    scRNA_seurat$celltype %in% unique(grep('^Macro',scRNA_seurat$celltype, value =T)),'Macrophage',ifelse(
        scRNA_seurat$celltype %in% "Mono_CD14", 'Monoc-CD14',ifelse(
            scRNA_seurat$celltype %in% 'Proliferating', 'Proliferating', ifelse(
                scRNA_seurat$celltype %in% 'Mast cells', 'Mast cells','DC'
            )
        )
    )
)
scRNA_seurat = scRNA_seurat[, scRNA_seurat$cellphonedb_celltype != 'Proliferating']
# 提取表达矩阵
counts.mye = as.data.frame(scRNA_seurat@assays$RNA@counts)
# 提取celltype
Idents(scRNA_seurat) = scRNA_seurat$cellphonedb_celltype
meta.mye=as.data.frame(scRNA_seurat@active.ident)
meta.mye$cell=row.names(meta.mye)
colnames(meta.mye) = c('cellphonedb_celltype','cell')
meta.mye=meta.mye[,c("cell","cellphonedb_celltype")]
qsavem(meta.mye,counts.mye,file = 'mye_cellphonedb.qs')

#<<<<<<<<<<<<<<<<<< 读入cancer cells >>>>>>>>>>>>>>>>>>>>>
setwd('/root/wangje/Project/刘老师/EPi/新结果/添加肝脏')
load('CCA_EPI.RData')
scRNA_seurat$cellphonedb_celltype = 'Cancer cells'
# 提取表达矩阵
counts.cancer = as.data.frame(scRNA_seurat@assays$RNA@counts)
# 提取celltype
Idents(scRNA_seurat) = scRNA_seurat$cellphonedb_celltype
meta.cancer=as.data.frame(scRNA_seurat@active.ident)
meta.cancer$cell=row.names(meta.cancer)
colnames(meta.cancer) = c('cellphonedb_celltype','cell')
meta.cancer=meta.cancer[,c("cell","cellphonedb_celltype")]
qsavem(meta.cancer,counts.cancer,file = 'cancer_cellphonedb.qs')

#<<<<<<<<<<<<<<<<<< 合并文件和分样本分开 >>>>>>>>>>>>>>>>>>>
mye = qread('/root/wangje/Project/刘老师/Myeloids/Data/mye_cellphonedb.qs')
nkt = qread('/root/wangje/Project/刘老师/NK_T/Data/nkt_cellphonedb.qs')
bcell = qread('/root/wangje/Project/刘老师/Bcells/Data/bcells_cellphonedb.qs')
fib = qread('/root/wangje/Project/刘老师/Fibroblasts/fib_cellphonedb.qs')
endo = qread('/root/wangje/Project/刘老师/Endothelials/endo_cellphonedb.qs')
cancer = qread('/root/wangje/Project/刘老师/EPi/新结果/添加肝脏/cancer_cellphonedb.qs')

all_meta = rbind(meta.nkt,mye$meta.mye, bcell$meta.bcells, fib$meta.fib, endo$meta.endo, cancer$meta.cancer)
meta = all_meta[all_meta$cell[!duplicated(all_meta$cell)],1:2]
counts = cbind(counts.nkt,mye$counts.mye, bcell$counts.bcells, fib$counts.fib, endo$counts.endo, cancer$counts.cancer)
sub_counts = counts %>% select(meta$cell)

setwd('/root/wangje/Project/刘老师/new_cpdb')
meta$sample = str_split_fixed(meta$cell, '_[A|T|G|C].*', n= 2)[,1]
for(sp in unique(meta$sample)){
    print(sp)
    dir.create(paste0('./',sp))
    sp_meta = meta[which(meta$sample == sp),c('cell','cellphonedb_celltype')]
    fwrite(sp_meta, file = paste0('./',sp,'/',sp,'_meta.txt'), row.names=F, quote = F, sep = '\t')
}

sub_counts$gene = rownames(sub_counts)
for(sp in unique(meta$sample)){
    print(sp)
    sp_counts = sub_counts %>% select(gene,matches(paste0('^',sp,'_[A|T|G|C].*')))
    fwrite(sp_counts, file = paste0('./',sp,'/',sp,'_counts.txt'), row.names=F, quote = F, sep = '\t')
}
 
#<<<<<<<<<<<<<<<<<< 运行cellphonedb >>>>>>>>>>>>>>>>>>>>>>>>
mtx_path=/root/wangje/Project/刘老师/new_cpdb/CJME/CJME_counts.txt
obs_path=/root/wangje/Project/刘老师/new_cpdb/CJME/CJME_meta.txt
output_path=/root/wangje/Project/刘老师/new_cpdb/CJME/
n_jobs=12

cellphonedb method statistical_analysis \
  --counts-data gene_name \
  --output-path ${output_path} \
  --threshold 0.1 \
  --threads ${n_jobs} \
  --iterations 1000 \
  --output-format csv \
  ${obs_path} \
  ${mtx_path}

####### for 循环 #########
for sample in $(dir /root/wangje/Project/刘老师/new_cpdb/)
do 
    echo "*************** ${sample} ********************"
    mtx_path=/root/wangje/Project/刘老师/new_cpdb/${sample}/${sample}_counts.txt
    obs_path=/root/wangje/Project/刘老师/new_cpdb/${sample}/${sample}_meta.txt
    output_path=/root/wangje/Project/刘老师/new_cpdb/${sample}/
    n_jobs=40

    cellphonedb method statistical_analysis \
    --counts-data gene_name \
    --output-path ${output_path} \
    --threshold 0.1 \
    --subsampling-log true \
    --threads ${n_jobs} \
    --iterations 1000 \
    --output-format csv \
    ${obs_path} \
    ${mtx_path}  
done

#--subsampling-log true #对于没有log转化的数据，还要加这个参数

#<<<<<<<<<<<<<<<<< 使用cellphonedb自带的绘图函数 >>>>>>>>>>>>>>>>>>>
path=/root/wangje/Project/刘老师/new_cpdb/Result
for sp in $(dir ${path})
do 
    echo "****************** ${sp} ***********************"
    cellphonedb plot dot_plot --means-path ${path}/${sp}/means.csv --pvalues-path ${path}/${sp}/pvalues.csv --output-path  ${path}/${sp}
    cellphonedb plot heatmap_plot  ${path}/${sp}/${sp}_meta.txt --pvalues-path ${path}/${sp}/pvalues.csv --output-path ${path}/${sp}

done

#<<<<<<<<<<<<<<<<< 使用R绘制cellphonedb热图 >>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(dplyr)
library(patchwork)
library(stringr)
library(tibble)
library(tidyfast)
library(tidyverse)
library(pheatmap)

## 多个样本一起分析
id1 <- c("CJME","CMDI","HDYU","HXZH","LCLI","WYJU","WZLA","ZXME","ZJLI_0116") # R_Pre
id2 <- c("CZYI","FHXHBS1","FYYI","HEJI","LAWE","ZEYI","ZFXI","LIPE")  # NR_Pre
id3 <- c("CJME_0707","CMDI_0624","HDYU_0720","HXZH_0220","LCLI_0623","WYJU_0122","ZXME_0223","ZJLI_0312") # R_Post
id4 <- c("CZYI_0702","FHXHBS2","HEJX","LAWE_0309","ZEYI_0204") # NR_Post
###############################################################################################
# 这里的程序是对每一个样本的cellphonedb结果中的pvale进行计数，计算pvalue的值小于0.05的，
# 计算两种celltype之间的配受体数目， 但是在合并的时候，由于不同样本间的dim不同，合并比较有难度
# 于是是有下面的程序，保留长格式后再合并
###############################################################################################
count_int = function(data) {
    tmp = data %>% select(12:last_col())
    tmp_long = tmp %>% 
        tidyr::pivot_longer(everything(),
               names_to = "group",
               values_to = "value")
    tmp_long$sig = tmp_long$value < 0.05
    tmp_wide = tmp_long %>% select(group,sig) %>% group_by(group) %>% dplyr::summarise(value = sum(sig))    
    tmp_wide$source = str_split_fixed(tmp_wide$group,'\\|', n = 2)[,1]
    tmp_wide$target = str_split_fixed(tmp_wide$group,'\\|', n = 2)[,2]
    count_final = tmp_wide %>% select(source, target, value)
    colnames(count_final) = c("SOURCE", "TARGET", "COUNT")
    count_final = count_final %>% tidyr::pivot_wider(names_from = TARGET, values_from =COUNT)
    count_mat = count_final %>% as.data.frame() 
    count_mat = tibble::column_to_rownames(count_mat, 'SOURCE')
    dcm <- diag(as.matrix(count_mat))
    count_mat <- count_mat + t(count_mat)
    diag(count_mat) <- dcm
    return(count_mat)
}
###############################################################################################
count_long = count_int = function(data) {
    tmp = data %>% select(12:last_col())
    tmp_long = tmp %>% 
        tidyr::pivot_longer(everything(),
               names_to = "group",
               values_to = "value")
    tmp_long$sig = tmp_long$value < 0.05
    tmp_wide = tmp_long %>% select(group,sig) %>% group_by(group) %>% dplyr::summarise(value = sum(sig))    
    tmp_wide$source = str_split_fixed(tmp_wide$group,'\\|', n = 2)[,1]
    tmp_wide$target = str_split_fixed(tmp_wide$group,'\\|', n = 2)[,2]
    count_final = tmp_wide %>% select(source, target, value)
    colnames(count_final) = c("SOURCE", "TARGET", "COUNT")
    return(count_final)
}

tidy_long = function(data){
    if(is.list(data)){
        res = Reduce(rbind, data)
        res = res %>% group_by(SOURCE, TARGET) %>% dplyr::summarise(vaue = sum(COUNT))
        colnames(res) = c('SOURCE', 'TARGET', 'COUNT')
        count_final = res %>% tidyr::pivot_wider(names_from = TARGET, values_from =COUNT)
        count_mat = count_final %>% as.data.frame() 
        count_mat = tibble::column_to_rownames(count_mat, 'SOURCE')
        dcm <- diag(as.matrix(count_mat))
        count_mat <- count_mat + t(count_mat)
        diag(count_mat) <- dcm
    }
    return(count_mat)
}

#%%%%%%%%%%%%% 查看不同组合的数目 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
generate_pvales = function(sample){
    flist_long = list()
    for(pwd in sample){
        print(pwd)
        org_path = "/root/wangje/Project/刘老师/new_cpdb/Reault"
        pfile = paste0()
        pvals = read.csv(,stringsAsFactors = F, check.names = F, header = T)
        # 保存长数据
        flist_long[[pwd]] = count_long(data = pvals)
        # 生成单个样本的宽数据
        mat = count_int(pvals)
        if (dir.exists(paste0(org_path,'/',pwd,'/out'))) {
        dir.create(paste0(org_path,'/',pwd,'/out'))
        }else{
            next
        }
        write.table(mat, file = paste0(org_path,'/',pwd,'/out/heatmap_count.txt'), sep = '\t', row.names =T, quote = F)
    }
    return(tidy_long(data = flist_long))
}

#<<<<<<<<<<<< 生成不同组合的数据 并且绘图>>>>>>>>>>>>>>>
R_Pre = generate_pvales(sample = id1)
R_Post = generate_pvales(sample = id3)
NR_Pre = generate_pvales(sample = id2)
NR_Post = generate_pvales(sample = id4)
#<<<<<<<<<<<<< 绘图 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
p_r_pre = pheatmap(
    R_Pre,
    main = 'R_Pre',
    cluster_rows = F,
    cluster_cols = F,
    display_numbers = TRUE,
    fondsize = 17
)

p_r_post = pheatmap(
    R_Post,
    main = 'R_Post',
    cluster_rows = F,
    cluster_cols = F,
    display_numbers = TRUE,
    fondsize = 17
)

p_nr_pre = pheatmap(
    NR_Pre,
    main = 'NR_Pre',
    cluster_rows = F,
    cluster_cols = F,
    display_numbers = TRUE,
    fondsize = 17
)

p_nr_post = pheatmap(
    NR_Post,
    main = 'NR_Post',
    cluster_rows = F,
    cluster_cols = F,
    display_numbers = TRUE,
    fondsize = 17
)

g1 = as.ggplot(p_r_pre)
g2 = as.ggplot(p_r_post)
g3 = as.ggplot(p_nr_pre)
g4 = as.ggplot(p_nr_post)

ggsave(filename = './cpdb_pheatmap.png', height =  8, width =  32, plot = g1|g2|g3|g4, bg = 'white')

# <<<<<<<<<<<<<<<<<< 绘制单个样本的热图 >>>>>>>>>>>>>>>>>>>>
ids4 <- c('CJME','CJME_0707','CMDI','CMDI_0624','HDYU','HDYU_0720','HXZH','HXZH_0220','LCLI','LCLI_0623','WYJU','WYJU_0122',
            'ZXME','ZXME_0223','ZJLI_0116','ZJLI_0312','CZYI','CZYI_0702','FHXHBS1','FHXHBS2','HEJI','HEJX','LAWE','LAWE_0309','ZEYI',
            'ZEYI_0204','WZLA','FYYI','ZFXI','LIPE')
ids5 <- c('RCJMEprePRLymph','RCJMEpostPRLymph','RCMDIprePRLymph','RCMDIpostPRLymph','RHDYUprePRBreast','RHDYUpostPRBrest','RHXZHprePRBreast',
            'RHXZHpostPRBreast','RLCLIpreCRLymph','RLCLIpostChestWall','RWYJUprePRLymph','RWYJUpostPRLymph','RZXMEprePRBreast','RZXMEpostPRBreast','RZJLIprePRBreast','RZJLIpostPRBreast',
            'NRCZYIprePDLiver','NRCZYIpostPDLiver','NRFHXHBSpreSDBreast','NRFHXHBSpostSDBreast','NRHEJIpreSDBreast','NRHEJXpostSDBreast','NRLAWEprePDLiver','NRLAWEpostPDLiver','NRZEYIprePDBreast','NRZEYIpostPDBreast',
            'RRWZLApreChestWall','NRFYYIpreSDLiver','NRZFXIprePDBreast','NRLIPEprePDLiver')
id = data.frame(ids4, ids5)

flist = list()
for(pwd in dirname){
    print(pwd)
    file = read.csv(paste0('/root/wangje/Project/刘老师/new_cpdb/',pwd,'/out/heatmap_count.txt'), sep = '\t', check.names = F)
    colnames(file) = c("SOURCE", "TARGET", "COUNT")
    count_final = file %>% tidyr::pivot_wider(names_from = TARGET, values_from =COUNT)
    count_mat = count_final %>% as.data.frame() 
    count_mat = tibble::column_to_rownames(count_mat, 'SOURCE')
    dcm <- diag(as.matrix(count_mat))
    count_mat <- count_mat + t(count_mat)
    diag(count_mat) <- dcm
    p = pheatmap(
        count_mat,
        main = id[which(id$ids4 == pwd),'ids5'],
        cluster_rows = F,
        cluster_cols = F,
        display_numbers = TRUE,
        fondsize = 17)
    flist[[pwd]] = as.ggplot(p)
}
p = wrap_plots(flist, ncol= 6)
ggsave(filename = './cpdb_sample.png', height =  30, width =  36, plot = p, bg = 'white')

# <<<<<<<<<<<<<<<<<<< 绘制网络图 >>>>>>>>>>>>>>>>>>>>>>>>>
library(psych)   # 这边一个心理学用的比较多的R包，可以绘制网络图
library(qgraph)
library(igraph)
library(dplyr)
library(tidyverse)

# 读入上一步生成的count-network文件
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

mynet<-read.delim('count_network.txt',check.names=FALSE)#读取数据
table(mynet$count)
mynet%>%filter(count>0)->mynet#过滤count为o的数据（有零会报错）
head(mynet)
net<-graph_from_data_frame(mynet)#构建net对象
plot(net)

karate_groups<-cluster_optimal(net)#统计每个端点的和
coords <- layout_in_circle(net,order=order(membership(karate_groups)))#设置网络布局
E(net)$width<- E(net)$count/10#根据count值设置边的宽
plot(net,edge.arrow.size=.1,
        edge.curved=0,
        vertex.color=allcolour,
        vertex.frame.color="#555555",
        vertex.label.color="black",
        layout = coords,
        vertex.label.cex=.7)


net2<-net#复制一份备用
for(i in 1:length(unique(mynet$SOURCE))){#配置发出端的颜色
        print(i)
        E(net)[map(unique(mynet$SoURCE),function(x){
        get.edge.ids(net,vp= c(unique(mynet$SOURCE)[i],x))
    })%>% unlist()]$color <- allcolour[i]
}

#这波操作谁有更好的解决方案？
plot(net,edge.arrow.size=.1,
        edge.curved=0.2,
        vertex.color=allcolour,
        vertex.frame.color="#555555",
        vertex.label.color="black",
        layout= coords,
        vertex.label.cex=.7)

################################################################    
### 循环绘制网络图
plot_net = function(path){
    for(sample in dir(path)){
        if(file.exists(paste0(sample,"/")))
    }

}




#####################################################################################################
##### 使用R包ktplots绘图
#<<<<<<<<<< 注意输入文件需要有Seurat对象挥着SingleCellExperienment对象 >>>>>>>>>>
### 生成Seurat对象
library(Seurat)
# library(ktplots) 
library(dplyr)
library(tidyverse)
library(stringr)
library(data.table)

## 单个样本不去除批次效应
counts = fread("CJME_counts.txt",data.table = F, header = T, sep = '\t')
rownames(counts) = counts$gene
counts = counts[,-1]
counts <- as(as.matrix(counts), "dgCMatrix")
scRNA = CreateSeuratObject(counts = counts)
scRNA[['percent.mt']] = PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA = scRNA %>% NormalizeData() %>%  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData()
scRNA = RunPCA(scRNA) 
scRNA <- FindNeighbors(scRNA, dims = 1:50)
scRNA <- FindClusters(scRNA, resolution = 0.5)
scRNA <- RunUMAP(scRNA, dims = 1:50)
# 添加celltype
cells = read.table('CJME_meta.txt', header = T, sep = '\t')
scRNA$celltype = cells[match(rownames(scRNA@meta.data),cells$cell),"cellphonedb_celltype"]
## 绘制和弦图

p <- plot_cpdb3(cell_type1 = 'CD8+ T cells|CD4+ T cells|Macrophage|Mast cells|DC', cell_type2 = 'Cancer cells',
    scdata = scRNA,
    idents = 'celltype', # column name where the cell ids are located in the metadata
    means = means,
    pvals = pvalues,
    deconvoluted = deconvoluted, # new options from here on specific to plot_cpdb3
    keep_significant_only = TRUE,
    standard_scale = TRUE,
    remove_self = TRUE
    )







scTCR_list.singleR <- GetAssayData(scRNA, solt = 'data')
clusters <- scRNA$RNA_snn_res.0.5
scTCR_list.singleRAll <- SingleR(
    test = scTCR_list.singleR,
    ref = list(encode, hema, hpca, immune, monaImmune),
    labels = list(
    encode$label.main,
    hema$label.main,
    hpca$label.main,
    immune$label.main,
    monaImmune$label.main),
    clusters = clusters,
    BPPARAM = BiocParallel::MulticoreParam(6))
celltype = data.frame(ClusterID = rownames(scTCR_list.singleRAll), celltype=scTCR_list.singleRAll$labels, stringsAsFactors = F)

## <<<<<<<<<<< Seurat v5 >>>>>>>>>>>>>>>>>>>>













############################################################
### 尝试直接合并pvalues文件
# combine_cpdb <- function(...) {
#     output <- list(...)
#     anames <- c("id_cp_interaction", "interacting_pair", "partner_a", "partner_b",
#         "gene_a", "gene_b", "secreted", "receptor_a", "receptor_b", "annotation_strategy",
#         "is_integrin")
#     bnames <- c("gene_name", "uniprot", "is_complex", "protein_name", "complex_name",
#         "id_cp_interaction")
#     if (all(colnames(output[[1]])[1:11] == anames)) {
#         out <- output %>% reduce(full_join, by = anames)
#     } else if (all(colnames(output[[1]])[1:6] == bnames)) {
#         out <- output %>% reduce(full_join, by = bnames)
#     }
#     return(out)
# }




#############################################################################################################################################################################################################
#### cellphoneDB v4
#############################################################################################################################################################################################################
# 创建虚拟环境，并且在虚拟环境中暗转需要的模块
















#=====================================================================================================
# 计算不同组合共有的，独有的interaction类型
get_prop_num = function(data1, data2){
    df = data.frame(matrix(nrow = 144, ncol = 4))
    colnames(df) = c('interaction','intersect_pair','data1_uniq','data2_uniq')
    df$interaction = unique(data1$group)
    # 计算不同的比例数目
        for (celltype in unique(data1$group)) {
            tmp1 = data1 %>% filter(group == celltype)
            tmp2 = data2 %>% filter(group == celltype)
            # cat('Data1维度: ', dim(tmp1),'\t','Data2维度: ', dim(tmp2), '\n')
            # 共有的interaction数目
            interaction_pair = intersect(tmp1$interacting_pair, tmp2$interacting_pair) %>% length()
            # 查看data1中独特的interaction pair数目
            diff_data1 = setdiff(tmp1$interacting_pair, tmp2$interacting_pair) %>% length()
            # 查看data2中独
            diff_data2 = setdiff(tmp2$interacting_pair, tmp1$interacting_pair) %>% length()
            # cat(celltype,interaction_pair, diff_data2, diff_data1, '\n')
            df[which(df$interaction==celltype),'intersect_pair'] = interaction_pair
            df[which(df$interaction == celltype),'data1_uniq'] = diff_data1
            df[which(df$interaction == celltype), 'data2_uniq'] = diff_data2
    }
    df = df %>% mutate(
    interaction_1 = str_split_fixed(interaction, '\\|', n= 2)[,1],
    interaction_2 = str_split_fixed(interaction, '\\|', n =2)[,2])
    return(df)
}

#+========================================================================================================
# 分别合并指定的配对，例如CD4+|CD8+,CD8+|CD4+这两个应该算一种，合并他们的数目
tidy_wider = function(data, colNum = 2){
    tmp = data %>% select(interaction_1, interaction_2, colNum)
    colnames(tmp) = c('SOURCE','TARGET','COUNT')
    df = tmp %>% pivot_wider(names_from = SOURCE, values_from = COUNT) %>% 
        as.data.frame() %>% 
        tibble::column_to_rownames('TARGET')
    df = df %>% select(rownames(df))
    dcm <- diag(as.matrix(df))
    count_mat <- df + t(df)
    diag(count_mat) <- dcm
    return(count_mat)
}

# 分别选择并且计算
merged_type = function(data){
    # 合并shared tcr
    shareed_interaction = tidy_wider(data, column = 2)
    shareed_interaction = shared_interaction %>% tibble::rownames_to_column('SOURCE') %>% pivot_longer(-1,names_to = 'TARGET',values_to = 'Shared')
    # data1中独有的interaction 
    data1_interaction = tidy_wider(df, colNum = 3)
    data1_interaction = data1_interaction %>% tibble::rownames_to_column('SOURCE') %>% pivot_longer(-1,names_to = 'TARGET',values_to = 'data1')
    # data2独有个体
    data2_interaction = tidy_wider(df, colNum = 4)
    data2_interaction = data2_interaction %>% tibble::rownames_to_column('SOURCE') %>% pivot_longer(-1,names_to = 'TARGET',values_to = 'data2')
    # 合并样本
    tmp =  merge(shared_interaction, data1_interaction, by = c('SOURCE','TARGET'))
    tmp = merge(tmp, data2_interaction, by = c('SOURCE','TARGET'))
    return(tmp)
}

#===============================================================================================================
# 去除重复的数据
Remove_Dup = function(df){
    df$sorted_combination <- apply(df[, c("SOURCE", "TARGET")], 1, function(x) paste(sort(x), collapse = " - "))
    df = df[!duplicated(df$sorted_combination), ]
    df$SUM = rowSum(df[c(3,4,5)])
    df = dplyr::arrange(df, -SUM) 
    colnames(df)[c(4,5)] = c('data1','data2')

    # 转换为长数据 
    tmp_long = tidy_R_Pre_NR_Pre %>% select(1:5) %>% pivot_longer(
    cols = c('Shared','data1','data2'),
    names_to = 'variable',
    values_to = 'value'
    )
    return(tmp.long)

}

#==================================================================================================================
# 绘图
plot_col = function(data){
    col = ifelse(grepl('Cancer cells', levels(data$x)), 'Red','black')
    p = data %>% ggplot(aes(x = x, y = value, fill = variable)) + 
    geom_col()+
    coord_flip() + 
    theme_classic() + 
    theme(axis.text.y = element_text(size = 15))+
    theme(axis.text.y = element_text(color = col),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 15, colour = 'black'),
          axis.title = element_text(size = 18, colour = 'black'),
          plot.title = element_text(hjust = 0.5, size=20))+
    scale_y_continuous(expand = c(0,0)) + 
    labs(title = 'R_Pre vs NR_Pre') 
    return(p)
}

#======================================================================================================================
# 绘制气泡图
## 转换pvalues格式
pvales = read.table('./pvalues.csv', sep = ',', header = T, check.names = F)
pvales_sub = pvales %>% select(interacting_pair,12:last_col())
pvales_sub_long = pvales_sub %>% pivot_longer(2:last_col())

## 转换means.txt
means = read.table('./means.csv', sep = ',', header = T, check.names = F)
means_sub = means %>% select(interacting_pair,12:last_col())
means_sub_long = means_sub %>% pivot_longer(2:last_col())


## 查找需要的组合
id1 <- c("CJME","CMDI","HDYU","HXZH","LCLI","WYJU","WZLA","ZXME","ZJLI_0116")
id2 <- c("CZYI","FHXHBS1","FYYI","HEJI","LAWE","ZEYI","ZFXI","LIPE")
id3 <- c("CJME_0707","CMDI_0624","HDYU_0720","HXZH_0220","LCLI_0623","WYJU_0122","ZXME_0223","ZJLI_0312")
id4 = c("CZYI_0702","FHXHBS2","HEJX","LAWE_0309","ZEYI_0204") #NR_Post
#// means 
flist = list()
for(index1 in id3){
    # 读入文件
    means = read.table('./means.csv', sep = ',', header = T, check.names = F)
    means_sub = means %>% select(interacting_pair,12:last_col())
    means_sub_long = means_sub %>% pivot_longer(2:last_col())
    flist[[index1]] = means_sub_long
}

flist1 = list()
for(i in 1:length(flist)){
    name = names(flist[i])
    # print(name)
    flist1[[name]] = flist[[i]][flist[[i]]$interacting_pair %in% list1$V1,]
    # 查找指定的组合
    flist1[[name]] = flist1[[name]] %>% filter(name %in% c('CD8+ T cells|DC','CD8+ T cells|Macrophage','CD8+ T cells|Monoc-CD14'))
}

#// pvalues
flist2 = list()
for(index1 in id1){
    # 读入文件
    means = read.table('./pvalues.csv', sep = ',', header = T, check.names = F)
    means_sub = means %>% select(interacting_pair,12:last_col())
    means_sub_long = means_sub %>% pivot_longer(2:last_col())
    flist2[[index1]] = means_sub_long
}

flist3 = list()
for(i in 1:length(flist2)){
    name = names(flist2[i])
    # print(name)
    flist3[[name]] = flist2[[i]][flist2[[i]]$interacting_pair %in% list1$V1,]
    # 查找指定的组合
    flist3[[name]] = flist3[[name]] %>% filter(name %in% c('CD8+ T cells|DC','CD8+ T cells|Macrophage','CD8+ T cells|Monoc-CD14'))
}

#// 合并两个文件
merge_flist = function(file_list){
    for(i in 1:length(file_list)){
        print(names(file_list[i]))
        file_list[[i]]$sample = names(file_list[i])
    }
    tmp = Reduce(rbind, file_list)
    return(tmp) 
}

tmp1 = merge_flist(flist1) # means
tmp2 = merge_flist(flist3) # pvalues
colnames(tmp1)[3] = 'means'
colnames(tmp2)[3] = 'pvalues'

tmp = full_join(tmp1, tmp2, by = c('interacting_pair','name','sample')) 

p1 = ggplot(tmp, aes(x = name, y = interacting_pair, size = -log10(pvalues), color = means))


## // 写成函数 //
plot_dotplot = function(ids){
    flist = list()
    for(index1 in ids){
        # 读入文件
        filename = paste0('/root/wangje/Project/刘老师/new_cpdb/Result/',index1,'/means.csv')
        means = read.table(filename, sep = ',', header = T, check.names = F)
        means_sub = means %>% select(interacting_pair,12:last_col())
        means_sub_long = means_sub %>% pivot_longer(2:last_col())
        flist[[index1]] = means_sub_long
    }

    flist1 = list()
    for(i in 1:length(flist)){
        name = names(flist[i])
        # print(name)
        flist1[[name]] = flist[[i]][flist[[i]]$interacting_pair %in% list1$V1,]
        # 查找指定的组合
        flist1[[name]] = flist1[[name]] %>% filter(name %in% c('CD8+ T cells|DC','CD8+ T cells|Macrophage','CD8+ T cells|Monoc-CD14'))
    }

    flist2 = list()
    for(index1 in ids){
        # 读入文件
        filename = paste0('/root/wangje/Project/刘老师/new_cpdb/Result/',index1,'/pvalues.csv')
        means = read.table(filename, sep = ',', header = T, check.names = F)
        means_sub = means %>% select(interacting_pair,12:last_col())
        means_sub_long = means_sub %>% pivot_longer(2:last_col())
        flist2[[index1]] = means_sub_long
    }

    flist3 = list()
    for(i in 1:length(flist2)){
        name = names(flist2[i])
        # print(name)
        flist3[[name]] = flist2[[i]][flist2[[i]]$interacting_pair %in% list1$V1,]
        # 查找指定的组合
        flist3[[name]] = flist3[[name]] %>% filter(name %in% c('CD8+ T cells|DC','CD8+ T cells|Macrophage','CD8+ T cells|Monoc-CD14'))
    }

    # 合并文件
    tmp1 = merge_flist(flist1) # means
    tmp2 = merge_flist(flist3) # pvalues
    colnames(tmp1)[3] = 'means'
    colnames(tmp2)[3] = 'pvalues'
    tmp = full_join(tmp1, tmp2, by = c('interacting_pair','name','sample')) 
    return(tmp)
}

data1 = plot_dotplot(id1)
data2 = plot_dotplot(id2)
data3 = plot_dotplot(id3)
data4 = plot_dotplot(id4)

p1 = ggplot(data1, aes(x = name, y = interacting_pair, size = -log10(pvalues), color = means)) + geom_point()+ theme(axis.text.x=element_text(angle=90, hjust=1,vjust=1)) + facet_wrap(.~sample, nrow = 1)
p2 = ggplot(data3, aes(x = name, y = interacting_pair, size = -log10(pvalues), color = means)) + geom_point() + theme(axis.text.x=element_text(angle=90, hjust=1,vjust=1)) +  facet_wrap(.~sample, nrow = 1)
p3 = ggplot(data2, aes(x = name, y = interacting_pair, size = -log10(pvalues), color = means)) + geom_point() + theme(axis.text.x=element_text(angle=90, hjust=1,vjust=1)) +  facet_wrap(.~sample, nrow = 1)
p4 = ggplot(data4, aes(x = name, y = interacting_pair, size = -log10(pvalues), color = means)) + geom_point() + theme(axis.text.x=element_text(angle=90, hjust=1,vjust=1)) +  facet_wrap(.~sample, nrow = 1)

# 保存图片
ggsave(filename = '/root/wangje/Project/刘老师/new_cpdb/R_pre.png',  width = 18, height = 10, plot = p1 + labs(title = 'R_Pre'))
ggsave(filename = '/root/wangje/Project/刘老师/new_cpdb/R_post.png',  width = 18, height = 10, plot = p2 + labs(title = 'R_Post'))
ggsave(filename = '/root/wangje/Project/刘老师/new_cpdb/NR_pre.png',  width = 18, height = 10, plot = p3 + labs(title = 'NR_Pre'))
ggsave(filename = '/root/wangje/Project/刘老师/new_cpdb/NR_post.png',  width = 18, height = 10, plot = p4 + labs(title = 'NR_Post'))





# pvals <- read.delim("pvalues.txt", check.names = FALSE)
# means <- read.delim("means.txt", check.names = FALSE)

# I've provided an example dataset
data(cpdb_output)

# pvals <- read.delim("pvalues.csv", check.names = FALSE, sep =',')
# means <- read.delim("means.csv", check.names = FALSE, sep = ',')
# test1 = plot_cpdb(cell_type1 = 'B cell', cell_type2 = 'CD4+ T cell', scdata = as.SingleCellExperiment(test),
# 	idents = 'celltype', # column name where the cell ids are located in the metadata
# 	split.by = 'Experiment', # column name where the grouping column is. Optional.
# 	means = means, pvals = pvals,
# 	genes = c("XCR1", "CXCL10", "CCL5")) +
# small_axis(fontsize = 12) + small_grid() + small_guide() + small_legend(fontsize = 12) # some helper functions included in ktplots to help with the plotting

# plot_cpdb(cell_type1 = 'B cells', cell_type2 = 'CD4+ T cells', scdata = as.SingleCellExperiment(test),
#     idents = 'celltype', means = means, pvals = pvals,
#     gene.family = 'chemokines') + small_guide() + small_axis() + small_legend(keysize=.5)

