
#@ analysis.R
#@ 分析脚本合集
suppressMessages({
    library(Seurat)
    library(ggplot2)
    library(SCP)
    library(patchwork)
    library(scCustomize)
    library(viridis)
    library(SCpubr)
    library(ggpubr)
    library(rstatix)
    library(SingleR)
    library(dplyr)
})

#@ 对Seurat进行singleR注释
#@ data为Seurat对象
#@ cluster为需要注释的分辨率,例如data$RNA_snn_res.0.1或data$integrated_snn_res.0.1
#@ resolution为注释的分辨率，例如:RNA_snn_res.0.1或integrated_snn_res.0.1

Run_singleR <- function(data, cluster, resolution){
    stopifnot(file.exists("~/singleR.RData"))
    load('~/singleR.RData')
    ref_list <- list(encode, hema, hpca, immune, monaImmune)
	labels_list <- list(
    encode$label.main,
    hema$label.main,
    hpca$label.main,
    immune$label.main,
    monaImmune$label.main)
	sce.data <- GetAssayData(data, solt = 'data')
	sce.singleR <- SingleR(
    test = sce.data,
    ref = ref_list,
    labels = labels_list,
    clusters = cluster,
    BPPARAM = BiocParallel::MulticoreParam(6)
  )
	# 提取 celltype 数据
  celltype <- data.frame(
    ClusterID = rownames(sce.singleR),
    singleR = sce.singleR$labels,
    stringsAsFactors = FALSE
  )
  # 将注释信息添加到seurat对象中
    if(is.null(resolution)) stop("missing resolution files.")
    tmp = data@meta.data %>% select(resolution)
	colnames(tmp) = "ClusterID"
	tmp1 = left_join(tmp, celltype, by = "ClusterID")
	rownames(tmp1) = rownames(tmp)
	sce_new = AddMetaData(data, tmp1)
	return(sce_new)
}

#@ 进行CCA聚类
Run_Seurat_CCA <- function(data, 
  split.by, 
  k.weight,
  vars.to.regress = c('nCount_RNA','percenmt.mt')){
  scRNAlist <- Seurat::SplitObject(data, split.by = split.by)
  scRNAlist <- lapply(scRNAlist, FUN = function(x){
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
  })
  features <- SelectIntegrationFeatures(object.list = scRNAlist)
  anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = features)
  anchors <<- anchors
  scRNA_CCA <- IntegrateData(anchorset = anchors, k.weight = k.weight)
  DefaultAssay(scRNA_CCA) <- "integrated"
  # Run the standard workflow for visualization and clustering
  scRNA_CCA <- ScaleData(scRNA_CCA,vars.to.regress = vars.to.regress, verbose = FALSE)
  scRNA_CCA <- RunPCA(scRNA_CCA, npcs = 40, verbose = FALSE)
  scRNA_CCA <- RunUMAP(scRNA_CCA, reduction = "pca", dims = 1:40)
  scRNA_CCA <- FindNeighbors(scRNA_CCA, reduction = "pca", dims = 1:40)
  scRNA_CCA <- FindClusters(scRNA_CCA, resolution = 0.5)
  return(scRNA_CCA)
}


plist <- list()
for (i in names(marker.list)) {
    for (j in marker.list[[i]]) {
        #    print(paste0(i,"_",j))
        tmp <- tryCatch({
            scCustomize::FeaturePlot_scCustom(
                scRNA,
                features = j,
                order = T,
                colors_use = colorRampPalette(c(
                            "#3288BD", "white", "red"
                        ))(50)
                ) +  
                theme(legend.position = "right") +
                labs(title = paste0(i, "_", j)) 
        },
        error = function(e) {
            message("Error @ ", j)
            return(NA)
        },
        finally = {
            message(paste0(i, "_", j, "_next..."))
        })
        plist[[paste0(i, "_", j)]] <- tmp
    }
}
plist <- lapply(plist, FUN = function(x) x + Seurat::NoAxes())
pp <- wrap_plots(plist, ncol = 6)


plist <- list()
for (i in names(marker.list)) {
    for (j in marker.list[[i]]) {
        #    print(paste0(i,"_",j))
        tmp <- tryCatch({
            scCustomize::FeaturePlot_scCustom(
                scRNA,
                features = j,
                order = T,
                colors_use = colorRampPalette(c(
                            "#3288BD", "white", "red"
                        ))(50)
                ) +  
                theme(legend.position = "right") +
                labs(title = paste0(i, "_", j)) 
        },
        error = function(e) {
            message("Error @ ", j)
            return(NA)
        },
        finally = {
            message(paste0(i, "_", j, "_next..."))
        })
        plist[[paste0(i, "_", j)]] <- tmp
    }
}
plist <- lapply(plist, FUN = function(x) x + Seurat::NoAxes())
pp <- wrap_plots(plist, ncol = 6)


plist <- list()
for (i in names(marker.list)) {
    for (j in marker.list[[i]]) {
        #    print(paste0(i,"_",j))
        tmp <- tryCatch({
            scCustomize::Plot_Density_Custom(
                scRNA,
                features = j,
                viridis_palette = "viridis"
                ) +  
                theme(legend.position = "right") +
                labs(title = paste0(i, "_", j)) 
        },
        error = function(e) {
            message("Error @ ", j)
            return(NA)
        },
        finally = {
            message(paste0(i, "_", j, "_next..."))
        })
        plist[[paste0(i, "_", j)]] <- tmp
    }
}
plist <- lapply(plist, FUN = function(x) x + Seurat::NoAxes())
pp <- wrap_plots(plist, ncol = 6)
