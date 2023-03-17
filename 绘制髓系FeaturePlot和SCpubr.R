library(Seurat)
library(ggplot2)
library(SCpubr)
library(dplyr)
library(patchwork)

MKIT <- c('KIT','TPSAB1','CPA3')
pDC_LILRA4 <- c('CXCR3','LILRA4','LILRB1','GZMB','IL3RA')
migDC <- c('CCR7','CCL19','CCL17')
LanghDC <- c('CD1A','CD207','RXRA')
CD1_CLEC9A <- c('CLEC9A','FLT3','IDO1','XCR1','BATF3')
CD2_CD1C <- c('CD1C','FCER1A','HLA-DQA1','CLEC10A','ITGAM')
CDC3_LAMP3 <- c('LAMP3','CCR7','FSCN1')
ASDC <- c('AXL','SIGLEC6','CD22')
Mono_CD14 <- c('CD14','FCN1','S100A9','S100A8')
Mono_CD16 <- c('FCGR3A','LST1','LILRB2')
Macro_INHBA <- c('IMHBA','IL1RA','CCL4')
Macro_NLRP3 <- c('NLRP3','EREG','IL1B')
Macro_LYVE1 <- c('LYVE1','PLTP','SEPP1')
Macro_C1QC <- c('C1QC','C1QA','APOE')
Keratinocyte <- c('KRT1','KRT10','KRTDAP')
Proloiferating <- c('MKI67','TOP2A','STMN1')
Check_point <- c('VSIR','VSIG4','LGALS9','CD274','CD273','PDCD1','SIGLEC10','PDCD1LG2')

marker.list <- list(MKIT = MKIT,
                    pDC_LILRA4 = pDC_LILRA4,
                    migDC = migDC,
                    LanghDC = LanghDC,
                    CD1_CLEC9A = CD1_CLEC9A,
                    CD2_CD1C = CD2_CD1C,
                    CDCD3_LAMP3 = CDC3_LAMP3,
                    ASDC = ASDC,
                    Mono_CD14 = Mono_CD14,
                    Mono_CD16 = Mono_CD16,
                    Macro_INHBA = Macro_INHBA,
                    Macro_NLRP3 = Macro_NLRP3,
                    Macro_LYVE1 = Macro_LYVE1,
                    Macro_C1QC = Macro_C1QC,
                    Keratinocyte = Keratinocyte,
                    Proliferating = Proliferating,
                    Check_point = Check_point)
    
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
                        FeaturePlot(scRNA_data, features = j) +
                            theme(legend.position = "none") +
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
                            theme(legend.position = "none") +
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
# FeaturePlot
png("/root/wangje/Project/刘老师/Myeloids/Fig/Myeloids_多个markerFeaturePlot.png",height =6000,width = 10000,res=300)
plotFeature(scRNA_data=scRNA_seurat,choose="Feature",col_num=10,marker.list=marker.list)
dev.off()

# SCpubr
png("/root/wangje/Project/刘老师/Myeloids/Fig/Myeloids_多个marker_scpubr.png",height =6000,width = 10000,res=300)
plotFeature(scRNA_data=scRNA_seurat,choose="SCpubr",col_num=10,marker.list=marker.list)
dev.off()
