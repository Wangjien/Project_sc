library(infercnv)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(AnnoProbe)
library()

# load data 
big.fastmnn = '/root/wangje/wangjien_new/liu-big/scRNA_FastMNN_new.RData'
load(big.fastmnn)

# 提取部分文件
scRNA_sub = scRNA[,scRNA$celltype %in% c('Fibroblasts','Endothelials','Epithelials')]
# 

