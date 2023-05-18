library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(clusterProfiler)


stewd("/root/wangje/Project/刘老师/Bcells/Data")
load('new_Bcells_CCA.RData')
########################################################################
# 获得差异基因
#########################################################################
Find_DEGs <- function(
  data,
  ident.1 = 'R_Post',
  ident.2 = 'R_Pre',
  group = 'Treat_assess',
  idents = 'celltype'
) {
  flist <- list()
  if (DefaultAssay(data) != 'RNA') {
    DefaultAssay(data) <- 'RNA'
    Idents(data) <- idents
  }
  tryCatch({
    for (celltype in unique(Idents(data))) {
      print(celltype)
      res <- FindMarkers(
        object = data,
        ident.1 = ident.1,
        ident.2 = ident.2,
        group.by = group,
        subset.ident = celltype,
        logfc.threshold = 0.05, # 进行GSEA分析的时候可以将这个参数设置为0
        # min.pct = 0,
        # thresh.use = 0.99,
        # test.use = MAST
        only.pos = FALSE
      )
      print(dim(res))
      flist[[celltype]] <- res
    }
  }, error = function(err) {
    cat("\033[1;31m出现错误\033[0m\n")
    return(NULL)
  }, finally = {
    print("运行结束")
  })
  return(flist)
}

# 比较组合
R_Post_R_Pre <- Find_DEGs(data = scRNA_seurat, group = 'Treat_assess')
NR_Post_NR_Pre <- Find_DEGs(data = scRNA_seurat, ident.1 = 'NR_Post', ident.2 = 'NR_Pre',group = 'Treat_assess')
NR_Post_R_Post <- Find_DEGs(data = scRNA_seurat, ident.1 = 'NR_Post', ident.2 = 'R_Post',group = 'Treat_assess')
NR_Pre_R_Pre <- Find_DEGs(data = scRNA_seurat, ident.1 = 'NR_Pre', ident.2 = 'R_Pre',group = 'Treat_assess')

###################################################################
# 进行GO富集分析
###################################################################

go_enrichment_analysis <- function(
    data,
    title1 = '',
    title2 = '') 
    {
        plist_up <- list()
        plist_down <- list()
        tryCatch({
                
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
        
        # plist_up[[i]] = barplot(tmp_up, showCategory = 15) + labs(title = paste0(i, '-->', title1))
        # plist_down[[i]] = barplot(tmp_down, showCategory = 15) + labs(title = paste0(i, '-->', title2))
        if (dim(as.data.frame(tmp_up))[1] >=15) {
            tmp_up = as.data.frame(tmp_up)
            tmp_up$LOG10padj = -log10(tmp_up$p.adjust)
            tmp_up_sub = tmp_up %>% arrange(desc(LOG10padj)) %>% dplyr::slice(1:15)   
            p <- ggplot(tmp_up_sub,aes(reorder(Description,LOG10padj),LOG10padj,fill=Description))+
                geom_col(fill = '#a54947') + 
                labs(title = stringr::str_c(title1,' ',i))+
                theme(axis.text  = element_text(size = 10),
                plot.title = element_text(face = 'bold'))+
                coord_flip()
            plist_up[[i]] = p
        }else{
            tmp_up = as.data.frame(tmp_up)
            tmp_up$LOG10padj = -log10(tmp_up$p.adjust)
            tmp_up_sub = tmp_up %>% arrange(desc(LOG10padj)) %>% dplyr::slice(1:nrow(tmp_up))   
            p <- ggplot(tmp_up_sub,aes(reorder(Description,LOG10padj),LOG10padj,fill=Description))+
                geom_col(fill = '#a54947') + 
                labs(title = stringr::str_c(title1,' ',i))+
                theme(axis.text  = element_text(size = 10),
                plot.title = element_text(face = 'bold'))+
                coord_flip()
            plist_up[[i]] = p
        }

        if(dim(as.data.frame(tmp_down))[1] >=15){
            tmp_down = as.data.frame(tmp_down)
            tmp_down$LOG10padj = -log10(tmp_down$p.adjust)
            tmp_down_sub = tmp_down %>% arrange(desc(LOG10padj)) %>% dplyr::slice(1:15)   
            p1 <- ggplot(tmp_down_sub,aes(reorder(Description,LOG10padj),LOG10padj,fill=Description))+
                geom_col(fill = '#a54947') + 
                labs(title = stringr::str_c(title2,' ',i))+
                theme(axis.text  = element_text(size = 10),
                plot.title = element_text(face = 'bold'))+
                coord_flip()
            plist_down[[i]] = p1
        }else{
            tmp_down = as.data.frame(tmp_down)
            tmp_down_sub = tmp_down %>% arrange(desc(LOG10padj)) %>% dplyr::slice(1:nrow(tmp_down))
            p1 <- ggplot(tmp_down_sub,aes(reorder(Description,LOG10padj),LOG10padj,fill=Description))+
                geom_col(fill = '#a54947') + 
                labs(title = stringr::str_c(title2,' ',i))+
                theme(axis.text  = element_text(size = 10),
                plot.title = element_text(face = 'bold'))+
                coord_flip()
            plist_down[[i]] = p1
        }
    }
        # },error = function(err){
        #     cat("\033[1;31m出现错误\033[0m\n")
        #     return(NULL)
        # },finally={
        #     print("运行结束")
        })

    result_list <- list(plist_up, plist_down)
    return(result_list)
}

R_Post_R_Pre_go = go_enrichment_analysis(data = R_Post_R_Pre, title1 = 'R_Post vs R_Pre up', title2 = 'R_Post vs R_Pre down')
NR_Post_NR_Pre_go = go_enrichment_analysis(data = NR_Post_NR_Pre, title1 = 'NR_Post vs NR_Pre up', title2 = 'NR_Post vs NR_Pre down')
NR_Post_R_Post_go = go_enrichment_analysis(data = NR_Post_R_Post, title1 = 'NR_Post vs R_Post up', title2 = 'NR_Post vs R_Post down')
NR_Pre_R_Pre_go = go_enrichment_analysis(data = NR_Pre_R_Pre, title1 = 'NR_Pre vs R_Pre up', title2 = 'NR_Pre vs R_Pre down')


## 保存图片
pwd = "/root/wangje/Project/刘老师/Bcells/Fig/富集结果/"
ggsave(filename=paste0(pwd,'R_Post_R_Pre_Go富集分析.png'), height = 16, width = 42, plot=wrap_plots(c(R_Post_R_Pre_go[[1]], R_Post_R_Pre_go[[2]]),ncol=4), bg = 'white')
ggsave(filename=paste0(pwd,'NR_Post_NR_Pre_Go富集分析.png'), height = 16, width = 42, plot=wrap_plots(c(NR_Post_NR_Pre_go[[1]], NR_Post_NR_Pre_go[[2]]),ncol=4), bg = 'white')
ggsave(filename=paste0(pwd,'NR_Post_R_Post_Go富集分析.png'), height = 16, width = 42, plot=wrap_plots(c(NR_Post_R_Post_go[[1]], NR_Post_R_Post_go[[2]]),ncol=4), bg = 'white')
ggsave(filename=paste0(pwd,'NR_Pre_R_Pre_Go富集分析.png'), height = 16, width = 42, plot=wrap_plots(c(NR_Pre_R_Pre_go[[1]], NR_Pre_R_Pre_go[[2]]),ncol=4), bg = 'white')

#################################################################
# KEGG 富集
#################################################################
kegg_enrichmernt_analysis <- function(data, title1 = '', title2 = ''){
    if(!typeof(data) == 'list'){
        stop('输入数据不是列表结构')
    }
        plist_up <- list()
        plist_down <- list()
        tryCatch({
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
            kegg_up = as.data.frame(kegg_up)
            kegg_down = as.data.frame(kegg_down)
            if (dim(kegg_up)[1] >=15) {
                kegg_up$LOG10padj = -log10(kegg_up$p.adjust)
                kegg_up_sub = kegg_up %>% arrange(desc(LOG10padj)) %>% dplyr::slice(1:15)   
                # print(kegg_up_sub)
                p <- ggplot(kegg_up_sub,aes(reorder(Description,LOG10padj),LOG10padj,fill=Description))+
                    geom_col(fill = '#a54947') + 
                    labs(title = stringr::str_c(title1,' ',i))+
                    theme(axis.text  = element_text(size = 10),
                    plot.title = element_text(face = 'bold'))+
                    coord_flip()
                plist_up[[i]] = p
            }else{
                # print(kegg_up)
                kegg_up$LOG10padj = -log10(kegg_up$p.adjust)
                kegg_up_sub = kegg_up %>% arrange(desc(LOG10padj)) %>% dplyr::slice(1:nrow(kegg_up))  
                print(kegg_down_sub) 
                p <- ggplot(kegg_up_sub,aes(reorder(Description,LOG10padj),LOG10padj,fill=Description))+
                    geom_col(fill = '#a54947') + 
                    labs(title = stringr::str_c(title1,' ',i))+
                    theme(axis.text  = element_text(size = 10),
                    plot.title = element_text(face = 'bold'))+
                    coord_flip()
                plist_up[[i]] = p
            }

            if(dim(kegg_down)[1] >=15){
                kegg_down$LOG10padj = -log10(kegg_down$p.adjust)
                kegg_down_sub = kegg_down %>% arrange(desc(LOG10padj)) %>% dplyr::slice(1:15)   
                p1 <- ggplot(kegg_down_sub,aes(reorder(Description,LOG10padj),LOG10padj,fill=Description))+
                    geom_col(fill = '#a54947') + 
                    labs(title = stringr::str_c(title2,' ',i))+
                    theme(axis.text  = element_text(size = 10),
                    plot.title = element_text(face = 'bold'))+
                    coord_flip()
                plist_down[[i]] = p1
            }else{
                kegg_down$LOG10padj = -log10(kegg_down$p.adjust)
                kegg_down_sub = kegg_down %>% arrange(desc(LOG10padj)) %>% dplyr::slice(1:nrow(kegg_down))
                p1 <- ggplot(kegg_down_sub,aes(reorder(Description,LOG10padj),LOG10padj,fill=Description))+
                    geom_col(fill = '#a54947') + 
                    labs(title = stringr::str_c(title2,' ',i))+
                    theme(axis.text  = element_text(size = 10),
                    plot.title = element_text(face = 'bold'))+
                    coord_flip()
                plist_down[[i]] = p1
            }   
        }
    })
    result_list <- list(plist_up, plist_down)
    return(result_list)
}

R_Post_R_Pre_kegg = kegg_enrichmernt_analysis(data = R_Post_R_Pre, title1 = 'R_Post vs R_Pre up', title2 = 'R_Post vs R_Pre down')
NR_Post_NR_Pre_kegg = kegg_enrichmernt_analysis(data = NR_Post_NR_Pre, title1 = 'NR_Post vs NR_Pre up', title2 = 'NR_Post vs NR_Pre down')
NR_Post_R_Post_kegg = kegg_enrichmernt_analysis(data = NR_Post_R_Post, title1 = 'NR_Post vs R_Post up', title2 = 'NR_Post vs R_Post down')
NR_Pre_R_Pre_kegg = kegg_enrichmernt_analysis(data = NR_Pre_R_Pre, title1 = 'NR_Pre vs R_Pre up', title2 = 'NR_Pre vs R_Pre down')
# 保存图片
pwd = "/root/wangje/Project/刘老师/Bcells/Fig/富集结果/"
ggsave(filename=paste0(pwd,'R_Post_R_Pre_kegg富集分析.png'), height = 16, width = 42, plot=wrap_plots(c(R_Post_R_Pre_kegg[[1]], R_Post_R_Pre_kegg[[2]]),ncol=4), bg = 'white')
ggsave(filename=paste0(pwd,'NR_Post_NR_Pre_kegg富集分析.png'), height = 16, width = 42, plot=wrap_plots(c(NR_Post_NR_Pre_kegg[[1]], NR_Post_NR_Pre_kegg[[2]]),ncol=4), bg = 'white')
ggsave(filename=paste0(pwd,'NR_Post_R_Post_kegg富集分析.png'), height = 16, width = 42, plot=wrap_plots(c(NR_Post_R_Post_kegg[[1]], NR_Post_R_Post_kegg[[2]]),ncol=4), bg = 'white')
ggsave(filename=paste0(pwd,'NR_Pre_R_Pre_kegg富集分析.png'), height = 16, width = 42, plot=wrap_plots(c(NR_Pre_R_Pre_kegg[[1]], NR_Pre_R_Pre_kegg[[2]]),ncol=4), bg = 'white')

### Reactome 富集
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
    tryCatch({
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
        )}
        # 绘图
        reactome_up = as.data.frame(reactome_up)
        reactome_down = as.data.frame(reactome_down)
        print(head(reactome_up))
        print(head(reactome_down))
        if(dim(reactome_up)[1] >=15){
            reactome_up$LOG10padj = -log10(reactome_up$p.adjust)
            reactome_up_sub = reactome_up %>% arrange(desc(LOG10padj)) %>% dplyr::slice(1:15)   
            # print(kegg_up_sub)
            p <- ggplot(reactome_up_sub,aes(reorder(Description,LOG10padj),LOG10padj,fill=Description))+
                geom_col(fill = '#a54947') + 
                labs(title = stringr::str_c(title1,' ',i))+
                theme(axis.text  = element_text(size = 10),
                plot.title = element_text(face = 'bold'))+
                coord_flip()
            plist_up[[i]] = p
        }else{
            reactome_up$LOG10padj = -log10(reactome_up$p.adjust)
            reactome_up_sub = reactome_up %>% arrange(desc(LOG10padj)) %>% dplyr::slice(1:nrow(reactome_up))  
            # print(reactome_up_sub) 
            p <- ggplot(reactome_up,aes(reorder(Description,LOG10padj),LOG10padj,fill=Description))+
                geom_col(fill = '#a54947') + 
                labs(title = stringr::str_c(title1,' ',i))+
                theme(axis.text  = element_text(size = 10),
                plot.title = element_text(face = 'bold'))+
                coord_flip()
            plist_up[[i]] = p

        }

        if(dim(reactome_down)[1] >=15){
        reactome_down$LOG10padj = -log10(reactome_down$p.adjust)
        reactome_down_sub = kegg_up %>% arrange(desc(LOG10padj)) %>% dplyr::slice(1:15)   
            # print(kegg_up_sub)
            p1 <- ggplot(reactome_down_sub,aes(reorder(Description,LOG10padj),LOG10padj,fill=Description))+
                geom_col(fill = '#a54947') + 
                labs(title = stringr::str_c(title2,' ',i))+
                theme(axis.text  = element_text(size = 10),
                plot.title = element_text(face = 'bold'))+
                coord_flip()
            plist_up[[i]] = p1
        }else{
            reactome_down$LOG10padj = -log10(reactome_down$p.adjust)
            reactome_down_sub = reactome_down %>% arrange(desc(LOG10padj)) %>% dplyr::slice(1:nrow(reactome_down))  
            # print(reactome_up_sub) 
            p1 <- ggplot(reactome_down_sub,aes(reorder(Description,LOG10padj),LOG10padj,fill=Description))+
                geom_col(fill = '#a54947') + 
                labs(title = stringr::str_c(title2,' ',i))+
                theme(axis.text  = element_text(size = 10),
                plot.title = element_text(face = 'bold'))+
                coord_flip()
            plist_up[[i]] = p1
        }

        # print(reactome_down)
        # print(reactome_up)
    # }, error = function(error){
    #     cat(i,"出现错误！！！","\n")
    #     return(NULL)
    # })

    })

    result_list <- list(plist_up, plist_down)
    return(result_list)
}


R_Post_R_Pre_Reactome <- reactome_enrichmernt_analysis(R_Post_R_Pre,title1 = 'R_Post vs R_Pre Up', title2 = 'R_Post vs R_Pre down')
NR_Post_NR_Pre_Reactome <- reactome_enrichmernt_analysis(NR_Post_NR_Pre,title1 = 'NR_Post vs NR_Pre Up', title2 = 'NR_Post vs NR_Pre down')
NR_Post_R_Post_Reactome <- reactome_enrichmernt_analysis(NR_Post_R_Post,title1 = 'NR_Post vs R_Post Up', title2 = 'NR_Post vs R_Post down')
NR_Pre_R_Pre_Reactome <- reactome_enrichmernt_analysis(NR_Pre_R_Pre,title1 = 'NR_Pre vs R_Pre Up', title2 = 'NR_Pre vs R_Pre down')
# 保存图片
ggsave(filename=paste0(pwd,'R_Post_R_Pre_Reactome富集分析.png'), height = 16, width = 42, plot=wrap_plots(c(R_Post_R_Pre_Reactome[[1]], R_Post_R_Pre_Reactome[[2]]),ncol=4), bg = 'white')
ggsave(filename=paste0(pwd,'NR_Post_NR_Pre_Reactome富集分析.png'), height = 16, width = 42, plot=wrap_plots(c(NR_Post_NR_Pre_Reactome[[1]], NR_Post_NR_Pre_Reactome[[2]]),ncol=4), bg = 'white')
ggsave(filename=paste0(pwd,'NR_Post_R_Post_Reactome富集分析.png'), height = 16, width = 42, plot=wrap_plots(c(NR_Post_R_Post_Reactome[[1]], NR_Post_R_Post_Reactome[[2]]),ncol=4), bg = 'white')
ggsave(filename=paste0(pwd,'NR_Pre_R_Pre_Reactome富集分析.png'), height = 16, width = 42, plot=wrap_plots(c(NR_Pre_R_Pre_Reactome[[1]], NR_Pre_R_Pre_Reactome[[2]]),ncol=4), bg = 'white')


