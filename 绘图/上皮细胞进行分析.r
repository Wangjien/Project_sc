library(seurat)
library(patchwork)
library(dplyr)
library(ggplot2)

# ----------------------------------------------------------------
# 使用fastmnn 的大群结果进行分分析
# ----------------------------------------------------------------

plist = list()
for(i in names(flist)){
    plist = DimPlot(flist[i], raster = F) + labs(title = i)
}