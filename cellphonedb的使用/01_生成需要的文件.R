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
for sp in $(dir /root/wangje/Project/刘老师/new_cpdb/)
do 
    echo "****************** ${sp} ***********************"
    cellphonedb plot dot_plot 
        --means-path /root/wangje/Project/刘老师/new_cpdb/${sp}/means.csv\ 
        --pvalues-path /root/wangje/Project/刘老师/new_cpdb/${sp}/pvalues.csv\ 
        --output-path  /root/wangje/Project/刘老师/new_cpdb/${sp}
    cellphonedb plot heatmap_plot 
        /root/wangje/Project/刘老师/new_cpdb/${sp}/${sp}_meta.txt\
        --pvalues-path /root/wangje/Project/刘老师/new_cpdb/${sp}/pvalues.csv\
        --output-path /root/wangje/Project/刘老师/new_cpdb/${sp}

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
        org_path = "/root/wangje/Project/刘老师/new_cpdb"
        pvals = read.csv(paste0(org_path,'/',pwd,'/pvalues.csv'),stringsAsFactors = F, check.names = F, header = T)
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





############################################################
### 尝试直接合并pvalues文件
flist1 = list()
for(pwd in list.files()){
    print(pwd)
    org_path = "/root/wangje/Project/刘老师/new_cpdb"
    pvals = read.csv(paste0(org_path,'/',pwd,'/pvalues.csv'),stringsAsFactors = F, check.names = F, header = T)
    flist1[[pwd]] = pvals
}