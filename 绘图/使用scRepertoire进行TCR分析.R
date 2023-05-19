library(scRepertoire)
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)


## 读入文件
setwd("/root/wangje/Project/刘老师/tcr/Data/rawData")
flist = list()
for(file in list.files('./',pattern='csv$')){
    name = stringr::str_split_fixed(file,'-',n=2)[1]
    # print(name)
    flist[[name]] = read.csv(file)
}

