#########################################
# Myeloids 差异基因GO富集分析
#########################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(Seurat)
library(tidyverse)
library(ReactomePA)
library(reactome.db)

#%%%%% 寻找差异基因
Find_DEGs <- function(data,
                      ident.1 = 'R_Post',
                      ident.2 = 'R_Pre',
                      group = 'Treat_assess',
                      idents = 'celltype') {
    flist = list()
    if (DefaultAssay(data) != 'RNA') {
        DefaultAssay(data) = 'RNA'
    }
    Idents(data) = idents
    for (celltype in unique(Idents(data))) {
        print(celltype)
        res = FindMarkers(
            object = data,
            ident.1 = ident.1,
            ident.2 = ident.2,
            group.by = group,
            subset.ident = celltype,
            logfc.threshold = 0.05,
            # 进行GSEA分析的时候可以将这个参数设置为0
            # min.pct = 0,
            # thresh.use = 0.99,
            # test.use = MAST
            only.pos = F
        )
        print(dim(res))
        flist[[celltype]] <- res
    }
    return(flist)
}

# ------> Myeloids
R_Post_R_Pre <- Find_DEGs(data = scRNA_seurat, group = 'Treat_assess')
NR_Post_NR_Pre <- Find_DEGs(data = scRNA_seurat, ident.1 = 'NR_Post', ident.2 = 'NR_Pre',group = 'Treat_assess')
NR_Post_R_Post <- Find_DEGs(data = scRNA_seurat, ident.1 = 'NR_Post', ident.2 = 'R_Post',group = 'Treat_assess')
NR_Pre_R_Pre <- Find_DEGs(data = scRNA_seurat, ident.1 = 'NR_Pre', ident.2 = 'R_Pre',group = 'Treat_assess')

#%%%%%%%%%%%%% GO 富集
go_enrichment_analysis <- function(data,
                                   title1 = '',
                                   title2 = '') {
    plist_up <- list()
    plist_down <- list()
    
    for (i in names(data)) {
        print(i)
        # 区分出上调和下调的基因 avg_log2FC
        data_up = data[[i]] %>% dplyr::filter(avg_log2FC > 0 &
                                                  p_val_adj < 0.05)
        data_down = data[[i]] %>% dplyr::filter(avg_log2FC < 0 &
                                                    p_val_adj < 0.05)
        
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


## 运行并且保存图片
### R_Post vs R_Pre
R_Post_R_Pre_go <- go_enrichment_analysis(R_Post_R_Pre, title1 = 'R_Post vs R_Pre up', title2 = 'R_Post vs R_Pre down')
png('/root/wangje/Project/刘老师/Myeloids/Fig/GO富集/20230424_R_Post_R_Pre_大群_GO.png',height = 8000,width = 10000, res=300)
wrap_plots(c(R_Post_R_Pre_go[[1]],R_Post_R_Pre_go[[1]]),ncol=5)
dev.off()


### NR_Post vs NR_Pre
NR_Post_NR_Pre_go <- go_enrichment_analysis(NR_Post_NR_Pre, title1 = 'NR_Post vs NR_Pre up', title2 = 'NR_Post vs NR_Pre down')
png('/root/wangje/Project/刘老师/Myeloids/Fig/GO富集/20230424_NR_Post_NR_Pre_大群_GO_up.png',height = 8000,width = 10000, res=300)
wrap_plots(c(NR_Post_NR_Pre_go[[1]],NR_Post_NR_Pre_go[[2]]),ncol=5)
dev.off()


### NR_Post vs R_Post
NR_Post_R_Post_go <- go_enrichment_analysis(NR_Post_R_Post, title1 = 'NR_Post vs R_Post up', title2 = 'NR_Post vs R_Post down')
png('/root/wangje/Project/刘老师/Myeloids/Fig/GO富集/20230424_NR_Post_R_Post_大群_GO.png',height = 8000,width = 10000, res=300)
wrap_plots(c(NR_Post_R_Post_go[[1]],NR_Post_R_Post_go[[2]]),ncol=5)
dev.off()

### NR_Pre vs R_Pre
NR_Pre_R_Pre_go <- go_enrichment_analysis(NR_Post_R_Post, title1 = 'NR_Post vs R_Post up', title2 = 'NR_Post vs R_Post down')
png('/root/wangje/Project/刘老师/Myeloids/Fig/GO富集/20230424_NR_Pre_R_Pre_大群_GO.png',height = 8000,width = 10000, res=300)
wrap_plots(c(NR_Pre_R_Pre_go[[1]],NR_Pre_R_Pre_go[[2]]),ncol=5)
dev.off()

####################################################
# Myeloids差异基因 KEGGE富集
####################################################
kegg_enrichmernt_analysis <- function(data, title1 = '', title2 = ''){
    if(!typeof(data) == 'list'){
        cat('输入数据不是列表结构')
    }
    plist_up <- list()
    plist_down <- list()
    flist_up <- list()
    flist_down <- list()
    for (i in names(data)) {
        print(i)
        # 提取出上下调的基因
        tmp_up <- data[[i]] %>% dplyr::filter(avg_log2FC > 0 & p_val_adj < 0.05)
        tmp_down <- data[[i]] %>% dplyr::filter(avg_log2FC < 0 & p_val_adj < 0.05)
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
        print(nrow(kegg_up))
        print(nrow(kegg_down))
        if(is.null(kegg_up)){
            next
        }else if(nrow(as.data.frame(kegg_up)) >=15){
            plist_up[[i]] <- barplot(kegg_up,showCategory = 15) + labs(title = paste0(i, '--->',title1))
        }else if(nrow(kegg_up)>0)
        {
            plist_up[[i]] <- barplot(kegg_up,showCategory = nrow(as.data.frame(kegg_up))) + labs(title = paste0(i, '--->',title1))
        }else if(nrow(keeg)==0){
            next
        }
        
        if(is.null(kegg_down)){
            next
        }else if(nrow(as.data.frame(kegg_down)) >= 15){
            plist_down[[i]] <- barplot(kegg_down,showCategory = 15) + labs(title = paste0(i, '--->',title2))
        }else if(nrow(kegg_down) > 0)
        {
            plist_down[[i]] <- barplot(kegg_down,showCategory = nrow(as.data.frame(kegg_down))) + labs(title = paste0(i, '--->',title2))
        }else if(nrow(kegg_down) == 0) {
            next
        }
        flist_up[[i]] <- as.data.frame(kegg_up)
        flist_down[[i]] <- as.data.frame(kegg_down)
    }
    result_list <- list(plist_up, plist_down,flist_up,flist_down)
    return(result_list)
}

#%%%%%%%% 运行并且保存图片
#%%%%%%%%%%%%%%%%%% R_Post vs R_Pre
R_Post_R_Pre_kegg <- kegg_enrichmernt_analysis(R_Post_R_Pre,title1 = 'R_Post vs R_Pre Up', title2 = 'R_Post vs R_Pre down')
p1 <- wrap_plots(R_Post_R_Pre_kegg[[1]],ncol = 6)
p2 <- wrap_plots(R_Post_R_Pre_kegg[[2]],ncol = 6)
png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/20230424_R_Post_R_Pre_KEGG_大群上调.png',height = 3000,width = 10000,res=300)
p1
dev.off()

png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/20230424_R_Post_R_Pre_KEGG_大群下调.png',height = 3000,width = 10000,res=300)
p2
dev.off()

#%%%%%%%%%%%%%%%%%% NR_Post vs NR_Pre
NR_Post_NR_Pre_kegg <- kegg_enrichmernt_analysis(NR_Post_NR_Pre,title1 = 'NR_Post vs NR_Pre Up', title2 = 'NR_Post vs NR_Pre down')
p1 <- wrap_plots(NR_Post_NR_Pre_kegg[[1]],ncol = 6)
p2 <- wrap_plots(NR_Post_NR_Pre_kegg[[2]],ncol = 6)
png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/20230424_NR_Post_NR_Pre_KEGG_大群上调.png',height = 3000,width = 10000,res=300)
p1
dev.off()

png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/20230424_NR_Post_NR_Pre_KEGG_大群下调.png',height = 3000,width = 10000,res=300)
p2
dev.off()

#%%%%%%%%%%%%%%%%%%% NR_Post vs R_Post
NR_Post_R_Post_kegg <- kegg_enrichmernt_analysis(NR_Post_R_Post,title1 = 'NR_Post vs R_Post Up', title2 = 'NR_Post vs R_Post down')
p1 <- wrap_plots(NR_Post_R_Post_kegg[[1]],ncol = 6)
p2 <- wrap_plots(NR_Post_R_Post_kegg[[2]],ncol = 6)
png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/20230424_NR_Post_R_Post_KEGG_大群上调.png',height = 3000,width = 10000,res=300)
p1
dev.off()

png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/20230424_NR_Post_R_Post_KEGG_大群下调.png',height = 3000,width = 10000,res=300)
p2
dev.off()

#%%%%%%%%%%%%%%%%%%%% NR_Pre vs R_Pre 
NR_Pre_R_Pre_kegg <- kegg_enrichmernt_analysis(NR_Pre_R_Pre,title1 = 'NR_Pre vs R_Pre Up', title2 = 'NR_Pre vs R_Pre down')
p1 <- wrap_plots(NR_Pre_R_Pre_kegg[[1]],ncol = 6)
p2 <- wrap_plots(NR_Pre_R_Pre_kegg[[2]],ncol = 6)
png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/20230424_NR_Pre_R_Pre_KEGG_大群上调.png',height = 3000,width = 10000,res=300)
p1
dev.off()

png('/root/wangje/Project/刘老师/Myeloids/Fig/KEGG富集/20230424_NR_Pre_R_Pre_KEGG_大群下调.png',height = 3000,width = 10000,res=300)
p2
dev.off()

###############################################
# Myeloids Reactome富集分析
###############################################
reactome_enrichmernt_analysis <- function(data,
                                          title1,
                                          title2){
    if(!typeof(data) == 'list'){
        cat('输入数据不是列表结构')
    }
    plist_up <- list()
    plist_down <- list()
    flist_up <- list()
    flist_down <- list()
    for (i in names(data)) {
        print(i)
        # 提取出上下调的基因
        tmp_up <- data[[i]] %>% dplyr::filter(avg_log2FC > 0 ,p_val < 0.05)
        tmp_down <- data[[i]] %>% dplyr::filter(avg_log2FC < 0 ,p_val < 0.05)
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
        if(is.null(reactome_up)){
            next
        }else if(nrow(as.data.frame(reactome_up)) >= 15){
            plist_up[[i]] <- barplot(reactome_up, showCategory = 15) + labs(title = paste0(i, '--->',title1))
        }else if(nrow(as.data.frame(reactome_up))>0 & nrow(as.data.frame(reactome_up)) < 15){
            plist_up[[i]] <- barplot(reactome_up, showCategory = nrow(as.data.frame(reactome_up))) + labs(title = paste0(i, '--->',title1))
        } else{
            next
        }
        
        if(is.null(reactome_down)){
            next
        }else if(nrow(as.data.frame(reactome_down))>=15){
            plist_down[[i]] <- barplot(reactome_down, showCategory = 15) + labs(title = paste0(i, '--->',title2))
        }else if(nrow(as.data.frame(reactome_down)) >0 & nrow(as.data.frame(reactome_down)) < 15){
            plist_down[[i]] <- barplot(reactome_down, showCategory = nrow(as.data.frame(reactome_down))) + labs(title = paste0(i, '--->',title2))
        }else{
            next
        } 
        flist_up[[i]] <- as.data.frame(reactome_up)
        flist_down[[i]] <- as.data.frame(reactome_down)
    }
    result_list <- list(plist_up,plist_down,flist_up, flist_down)
    return(result_list)
}
R_Post_R_Pre_Reactome <- reactome_enrichmernt_analysis(data = R_Post_R_Pre,title1 = 'R_Post va R_Pre up',title2 = 'R_Post vs R_Pre down')
# R_Post vs R_Pre
png('/root/wangje/Project/刘老师/Myeloids/Fig/Reactome富集/20230424_R_Post_R_pre_Reactome_Padj0.05_大群.png',height = 10000,width = 12000,res=300)
wrap_plots(c(R_Post_R_Pre_Reactome[[1]],R_Post_R_Pre_Reactome[[1]]),ncol=5)
dev.off()

# NR_Post vs NR_Pre
NR_Post_NR_Pre_Reactome <- reactome_enrichmernt_analysis(data = NR_Post_NR_Pre,title1 = 'NR_Post vs NR_Pre up',title2 = 'NR_Post vs NR_Pre down')

png('/root/wangje/Project/刘老师/Myeloids/Fig/Reactome富集/20230424_NR_Post_NR_pre_Reactome_Padj0.05_大群.png',height = 10000,width = 12000,res=300)
wrap_plots(c(NR_Post_R_Post_Reactome[[1]],NR_Post_R_Post_Reactome[[2]]),ncol=5)
dev.off()

# NR_Post vs R_Post
NR_Post_R_Post_Reactome <- reactome_enrichmernt_analysis(data = NR_Post_R_Post,title1 = 'NR_Post vs R_Post up',title2 = 'NR_Post vs R_Post down')
png('/root/wangje/Project/刘老师/Myeloids/Fig/Reactome富集/20230424_NR_Post_R_post_Reactome_Padj0.05_大群.png',height = 10000,width = 12000,res=300)
wrap_plots(c(NR_Post_R_Post_Reactome[[1]],NR_Post_R_Post_Reactome[[2]]),ncol =5)
dev.off()

# NR_Pre vs R_Pre
NR_Pre_R_pre_Reactome <- reactome_enrichmernt_analysis(data = NR_Pre_R_Pre,title1 = 'NR_Pre vs R_Pre up',title2 = 'NR_Pre vs R_Pre down')
png('/root/wangje/Project/刘老师/Myeloids/Fig/Reactome富集/20230424_NR_Pre_R_pre_Reactome_Padj0.05_大群.png',height = 10000,width = 12000,res=300)
wrap_plots(c(NR_Post_R_Post_Reactome[[1]],NR_Post_R_Post_Reactome[[2]]),ncol=5)
dev.off()

#%%%%%%%%%%%%%%%%%%% 输出差异基因文件
add_info <- function(data){
for(i in seq_len(12))
    data[[i]]$group = names(data[i])
    return(data)
}
    
Reduce(rbind,add_info(R_Post_R_Pre)) %>% write.csv(file = '/root/wangje/Project/刘老师/Myeloids/Data/R_Post_R_pre_差异基因.csv',row.names = T,quote = F)   
Reduce(rbind,add_info(NR_Post_NR_Pre)) %>% write.csv(file = '/root/wangje/Project/刘老师/Myeloids/Data/NR_Post_NR_pre_差异基因.csv',row.names = T,quote = F) 
Reduce(rbind,add_info(NR_Pre_R_Pre)) %>% write.csv(file = '/root/wangje/Project/刘老师/Myeloids/Data/NR_Pre_R_pre_差异基因.csv',row.names = T,quote = F) 
Reduce(rbind,add_info(NR_Post_R_Post)) %>% write.csv(file = '/root/wangje/Project/刘老师/Myeloids/Data/NR_Post_R_post_差异基因.csv',row.names = T,quote = F) 
