library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(readr)
library(cowplot)
library(reticulate)

# 导入python模块
plot.list = list()
scanpy = import("scanpy")
celltypist = import("celltypist")
pandas <- import("pandas")
numpy = import("numpy")

# 批量注释函数
celltypist_vis = function(adata,
                          seurat_Data,
                          model = Immune_All_Low.pkl,
                          title = "Immune_All_Low"){
  ### 开始预测
  predictions = celltypist$annotate(adata, model = model, majority_voting = T)
  # predictions$predicted_labels %>% head()
  
  #### 把这些信息加入到seurat对象中去
  seurat_Data  = AddMetaData(seurat_Data , predictions$predicted_labels) 
  head(seurat_Data )
  
  seurat_Data  = SetIdent(seurat_Data , value = "majority_voting")
  p.umap = DimPlot(seurat_Data, reduction = "umap", label = TRUE, pt.size = 0.5,label.box = F,repel = T) + 
        ggtitle(title)
  return(p.umap)
}

celltypist_anno = function(sce){
    model_type = list.files("/root/wangje/Reference/cellTypist_Ref/Human/", pattern="pkl$")
    names(model_type) = str_split(string = model_type,pattern = "\\.", simplify = T)[,1]
    pwd = "/root/wangje/Reference/cellTypist_Ref/Human/"
    model_list = lapply(model_type, function(x){
        celltypist$models$Model$load(model = paste0(pwd,x))})
    if(Reduce('|', is(sce) == "Seurat")){
        # 将Seurat对象转换为assdata对象
        adata = scanpy$AnnData(X = numpy$array(as.matrix(t(as.matrix(sce[['RNA']]@counts)))),
                       obs = pandas$DataFrame(sce@meta.data),
                       var = pandas$DataFrame(data.frame(gene = rownames(sce[['RNA']]@counts),
                                                         row.names = rownames(sce[['RNA']]@counts)))
        )
        scanpy$pp$normalize_total(adata, target_sum=1e4)
        scanpy$pp$log1p(adata)

        for (i in 1:length(model_list)) {
            print(i)
            plot.list[[i]] = celltypist_vis(adata = adata,
                                    seurat_Data = sce,
                                    model = model_list[[i]],
                                    title = names(model_list[i]))  
            }
        plot.list = plot.list[is.na(plot.list)]
        plot.list[[18]] = DimPlot(sce, raster=F)
        return(plot.list)        
    }
}

# ----------------------------------------------------------------
# 读入Seurat对象
setwd('')