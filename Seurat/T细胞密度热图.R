######################### T细胞绘制密度热图 ##############################
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(SCpubr)
## 读入文件
# load("/data/Fat01home/wangje/Clincal_data/NKT_seurat_CCA.RData")
# load("/data/Fat01home/wangje/Clincal_data/scRNA_seurat.RData")
# load("/data/Fat01home/wangje/Clincal_data/NKT_scRNA_harmony.RData")
## 绘制marker密度热图
Tcells <-c('CD3G','CD4','CD8A','CD8B')
NKcells <-c('KLRC1','XCL1','AREG','TYROBP','FCGR3A','FGFBP2')
Treg <- c('FOXP3','SAT1','IL2RA','IKZF2',"CTLA4")
TFH  <- c('CXCL13','BCL6','EBI3','CD200','ICOS','TOX','CXCR5','STAT3','TNFSF4','SLAMF1')
TH1 <-c('CXCL13','CXCR3','IFNG',"BHLHE40","GZMB","PDCD1","HAVCR2","ICOS","IGFLR1",'CCR1','TBX21','RUNX2','RUNX3','STAT1','TNF')
TH17<- c('IL23R','RORC','IL17A','IL17B','FURIN','CTSH','KLRB1','CCR6','CTSH','CAPG','FURIN','ITGAE','CD5L')
TEX <- c('GZMB','LAG3','HAVCR2','PDCD1','ENTPD1','ITGAE','PRDM1')
Naive <- c('CCR7','SELL','TCF7','CD27','CD28','S1PR1','IL7R')
Proliferating <-c('MKI67','TOP2A','STMN1','TUBA1B','HIST1H4C')
T_activate <-c('TNFRSF9','TNFRSF18')
Tumor_reactive <-c('ENTPD1','ITGAE')
PreDysfunctional <-c('GZMK','FYN','KLRG1','FCRL6','TXNIP','LYAR','CXCR3','EOMES','ZNF683')
DysfunctionalTF <-c('EOMES','TBX21','MAF','TOX','TCF7')

test <- list(Tcells=Tcells,
    NKcells=NKcells,
    Treg=Treg,
    TFH=TFH,
    TH17=TH17,
    TH1=TH1,
    TEX=TEX,
    Naive=Naive,
    Proliferating=Proliferating,
    Tumor_reactive=Tumor_reactive,
    PreDysfunctional=PreDysfunctional,
    DysfunctionalTF=DysfunctionalTF)


##-----------------------------------------------------------------
# 绘制FeaturePlot
##-----------------------------------------------------------------
plotFeature <- function(scRNA_data = scRNA_data,
                        choose = "Feature",
                        col_num = 6, marker.list = marker.list,...) {
    pacman::p_load("Seurat", "ggplot2", "tidyverse")
    DefaultAssay(scRNA_data) <- "RNA"
    plist <- list()
    if (is.null(choose)) {
        message("请选择绘图类型")
    } else if (choose == "Feature") {
        for (i in names(marker.list)) {
            for (j in marker.list[[i]]) {
                #    print(paste0(i,"_",j))
                tmp <- tryCatch(
                    {
                        FeaturePlot(scRNA_data, features = j,raster=F) +
                            theme(legend.position = "right") +
                            labs(title = paste0(i, "_", j))
                    },
                    error = function(e) {
                        message("Error @ ", j)
                        return(NA)
                    },
                    finally = {
                        message(paste0(i, "_", j, "_next..."))
                    }
                )
                plist[[paste0(i, "_", j)]] <- tmp
            }
        }
        p_new <- Filter(Negate(anyNA), plist)
        p <- wrap_plots(p_new, bycol = T, ncol = col_num)
        return(p)
    } else if (choose == "SCpubr") {
        for (i in names(marker.list)) {
            for (j in marker.list[[i]]) {
                #    print(paste0(i,"_",j))
                tmp <- tryCatch(
                    {
                        SCpubr::do_NebulosaPlot(scRNA_data, features = j, viridis_color_map = "H", pt.size = 0.02) +
                            theme(legend.position = "righty") +
                            labs(title = paste0(i, "_", j))
                    },
                    error = function(e) {
                        message("Error @ ", j)
                        return(NA)
                    },
                    finally = {
                        message(paste0(i, "_", j, "_next..."))
                    }
                )
                plist[[paste0(i, "_", j)]] <- tmp
            }
        }
        p_new <- Filter(Negate(anyNA), plist)
        p <- wrap_plots(p_new, bycol = T, ncol = col_num)
        return(p)
    }
    ...
}

png('./NK-T_FeaturePlot.png', height = 9000, width = 10000, res =300)
plotFeature(scRNA, col_num=10)
dev.off()






























### Seurat CCA
plist <- list()
for(i in names(test)){
    for (j in test[[i]]) {
    #    print(paste0(i,"_",j))
         tmp <- tryCatch(
            {SCpubr::do_NebulosaPlot(scRNA_seurat,features=j,viridis_color_map="H", pt.size=0.02)+
                theme(legend.position="none")+
                labs(title=paste0(i,"_",j))
                },
            error = function(e) { message('Error @ ',j) ; return(NA) },
            finally = { message(paste0(i,"_",j,'_next...'))}
        )
        plist[[paste0(i,"_",j)]] <- tmp
    }
}
p_new <- Filter(Negate(anyNA),plist)
p <- wrap_plots(p_new, bycol = T, ncol = 10)

png("GSM4288842_CCA_SCpubr_NEW.png",height = 9000,width=10000,res=300)
p
dev.off()

#### FeaturePlot
plist <- list()
for(i in names(test)){
    for (j in test[[i]]) {
    #    print(paste0(i,"_",j))
         tmp <- tryCatch(
            {FeaturePlot(scRNA_seurat,features=j)+
                theme(legend.position="none")+
                labs(title=paste0(i,"_",j))
                },
            error = function(e) { message('Error @ ',j) ; return(NA) },
            finally = { message(paste0(i,"_",j,'_next...'))}
        )
        plist[[paste0(i,"_",j)]] <- tmp
    }
}
p_new <- Filter(Negate(anyNA),plist)
p <- wrap_plots(p_new, bycol = T, ncol = 10)

png("GSM4288842_FastMNN_FeaturePlot.png",height = 9000,width=10000,res=300)
p
dev.off()

## 不同分辨率的DimPlot
plist <- lapply(paste0("integrated_snn_res.",seq(0.05,1,0.05)), function(res){
    print(res)
    DimPlot(scRNA_seurat,group.by=res,label=T)+
        ggplot2::theme(legend.position="none")
})
p <- wrap_plots(plist, bycol = T, ncol = 5)
png("Merge_Her2_DimPlot.png",height = 4000,width=5000,res=300)
p
dev.off()
### FastMNN
plist <- lapply(paste0("RNA_snn_res.",seq(0.05,1,0.05)), function(res){
    print(res)
    DimPlot(scRNA,group.by=res,label=T)+
        ggplot2::theme(legend.position="none")
})
p <- wrap_plots(plist, bycol = T, ncol = 5)
png("merge_FastMNN_DimPlot.png",height = 4000,width=5000,res=300)
p
dev.off()

save(scRNA_seurat,file="merge_3HER2_6HER2_5HER2_CD4.RData")
### conos
flist <- list()
for(i in unique(scRNA_seurat$Patient)){
    # tmp <- scRNA_seurat[,scRNA_seurat$Patient %in% i]@assays$RNA@counts
    tmp <- scRNA_seurat[,scRNA_seurat$Patient %in% i]
    flist[[i]] <- tmp
}
# panel.preprocessed.seurat <- lapply(flist, basicSeuratProc)
con <- Conos$new(flist, n.cores=10) 
space='PCA' # 可以选择 PCA, CPCA,CCA
con$buildGraph(k=30, k.self=5, 
               space=space,  # PCA, CPCA,CCA
               ncomps=30, 
               n.odgenes=2000, 
               matching.method='mNN', 
               metric='angular', 
               score.component.variance=TRUE, 
               verbose=TRUE)
con$findCommunities(method=leiden.community, resolution=seq(0.05,1,0.05)) # 相当于Seurat包中的FindClusters函数



adata_seurat = read_h5ad(file = './GSE156625_HCCscanpyobj.h5ad', 
        target.object = 'seurat', 
        assay_name = 'RNA')



#### FeaturePlot
Features <- c('METTL3','YTHDF1','YTHDF2','YTHDF3','METTL14','WTAP')
plist <- list()
for(i in Features){
#    print(paste0(i,"_",j))
        tmp <- tryCatch(
        {FeaturePlot(scRNA_seurat,features=i)+
            theme(legend.position="none")+
            labs(title=paste0(i))
            },
        error = function(e) { message('Error @ ',i) ; return(NA) },
        finally = { message(paste0(i,"_",'next...'))}
    )
    plist[[i]] <- tmp
}
p_new <- Filter(Negate(anyNA),plist)
p <- wrap_plots(p_new, bycol = T, ncol = 5)
png(file="./NKT_3HER2_FastMNN_FeaturePlot.png",height=1000,width=6000,res=300)
p
dev.off()
