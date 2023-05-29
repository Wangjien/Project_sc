Tcells <- c('CD3D','CD3G','CD2')
NK_cells <- c('KLRC1','KLRC3','TRDV9')
Fibroblasts <- c('COL1A1','DCN','LUM')
Myeloids <- c('LYZ','CD68','TYROBP')
Epithelial <- c('CD24','KRT19','EPCAM')
Bcells <- c('CD79A','CD19','MS4A1')
Endothelial <- c('CLDN5','FLT1','RAMP2')
Plasma <- c('IGHG1','JCHAIN','MZB1')
Hepatocytes <- c('ALB','APOB','HP')
Keratinocytes <- c("KRT5","KRT14","FABP5")
DC <- c("LILRA4","CXCR3","IRF7")
Mast <- c("CPA3","TPSABT","TPSB2")

marker.list <- list("NK cell"=NK_cells,"T cell"=Tcells,'B cell'=Bcells,"Plasmas" =Plasma,Myeloids=Myeloids,Fibroblasts=Fibroblasts,
                    Epithelials=Epithelial,Endothelials=Endothelial,Hepatocytes=Hepatocytes,Keratinocytes=Keratinocytes,
                    DC = DC, Mast = Mast)
                    
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
png("/root/wangje/Project/刘老师/3HER2数据/5HER2文章中的3HER2数据/Fig/Myeloids_多个大群markerFeaturePlot.png",height =6000,width = 6000,res=300)
plotFeature(scRNA_data=scRNA_CCA,choose="Feature",col_num=6,marker.list=marker.list)
dev.off()

# SCpubr
png("/root/wangje/Project/刘老师/3HER2数据/5HER2文章中的3HER2数据/Fig/Myeloids_多个大群marker_scpubr.png",height =6000,width = 6000,res=300)
plotFeature(scRNA_data=scRNA_CCA,choose="SCpubr",col_num=6,marker.list=marker.list)
dev.off()
                    
