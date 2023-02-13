# 查看时候混杂了其他的细胞
NKT <- c("CD3D", "CD3G", "CD2")
Fibroblasts <- c("COL1A1", "DCN", "LUM")
Myeloids <- c("LYZ", "CD68", "TYROBP")
Epithelial <- c("CD24", "KRT19", "EPCAM")
Bcells <- c("CD79A", "CD19", "MS4A1")
Endothelial <- c("CLDN5", "FLT1", "RAMP2")
Plasma <- c("IGHG1", "JCHAIN", "MZB1")
Hepatocytes <- c("ALB", "APOB", "HP")
Mast <- c('CPA3','TPSAB1','TPSB2')
DC <- c('LILRA4','CXCR3','IRF7')
levels <- c("NK&T cells", "B cells", "Plasmas", "Myeloids", "Epithelials", "Endothelials", "Fibroblasts", "Hepatocytes")
marker.list <- list(
    "NK&T cell" = NKT, "B cells" = Bcells, Plasmas = Plasma, Myeloids = Myeloids,
    Epithelials = Epithelial, Endothelials = Endothelial, Fibroblasts = Fibroblasts, Hepatocytes = Hepatocytes
)
# ----------------------- 写成函数 ---------------------------
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
png("/root/wangje/Project/刘老师/Fibroblasts/fastMNN/查看成纤维细胞是否有其他的细胞_FeaturePlot.png",height =4000,width = 6000,res=300)
plotFeature(scRNA_data=scRNA,choose="Feature",col_num=6,marker.list=marker.list)
dev.off()
