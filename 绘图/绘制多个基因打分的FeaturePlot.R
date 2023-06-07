library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)

##<<<<<<<<<<<<<<<<<<< 大群FeaturePlot >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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

markers.big = list(
    NK_Tcells = c(Tcells, NK_cells),
    B_Plasma = c(Bcells, Plasma),
    Myeloids = c(Myeloids,Mast, DC),
    Fibroblasts = Fibroblasts,
    Epithelial = Epithelial,
    Endothelial = Endothelial
)

plot_Feature = function(data, marker.list){
    if(!(DefaultAssay(data) == 'RNA')){
        DefaultAssay(data) = 'RNA'
    }
    if(length(marker.list) == 0){
        cat("\033[31mmarker list is empty!\033[0m\n")
    }else{
        for(i in seq_len(length(marker.list))){
           data =  
        }
    }


}

