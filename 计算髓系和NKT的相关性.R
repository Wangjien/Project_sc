
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)

# 读入NK & T cells
#########################################################################
################# 计算相关性，保留增值细胞
# load('/root/wangje/Project/刘老师/NK_T/Data/NKT_scRNA_seurat.RData')
NK_T_Fraction <- read.table('/root/wangje/Project/刘老师/NK_T/Data/合并TH17和Naive_结合_3HER2_以个样本大群的细胞数目为分母.txt', header=T, sep= '\t', fill =T, stringsAsFactors = F )
Myeloids_T_Fraction  <- read.table('/root/wangje/Project/刘老师/Myeloids/Data/以大群每个样本的细胞数目作为比例进行计算.txt', header=T ,sep ='\t', fill =T)
# 修改NK&T细胞亚群中的部分细胞名称
NK_T_Fraction$celltype <- gsub('Tdysp','Texp',NK_T_Fraction$celltype)
NK_T_Fraction$celltype <- gsub('TH1 cells','TH1-like',NK_T_Fraction$celltype)

# 提取文件中的Post数据
NKT_Sub <- NK_T_Fraction %>% select(Var1,celltype,Fraction,group) %>% 
    dplyr::rename_with(~'celltype1',2) %>% 
    filter(!Var1  %in% c("Three_CID3838--TNNaive","Three_CID3921--TNNaive","Three_CID45171--TNNaive"))

Myeloid_Sub <- Myeloids_T_Fraction %>% select(Var1,celltype,Fraction,group) %>% 
    dplyr::rename_with(~c("celltype2","Fraction2"),c(2,3))
Myeloid_Sub$Var1 <- gsub('CRPR-|SDPD-','',Myeloid_Sub$Var1)
merge_data <- cbind(NKT_Sub,Myeloid_Sub)

plist <- list()
for(i in unique(NKT_Sub$celltype1)){
    for(j in unique(Myeloid_Sub$celltype2)){
        tmp1 <- NKT_Sub %>% filter(celltype1 == i,group %in% c('R_Post','NR_Post'))
        tmp2 <- Myeloid_Sub %>% filter(celltype2 == j, group %in% c('R_Post','NR_Post'))
        tmp <- full_join(tmp1,tmp2,by=c('Var1','group'))
        print(head(tmp))
        p <- ggpubr::ggscatter(tmp,x='Fraction',y='Fraction2',
                            add = 'reg.line',conf.int = TRUE,
                            add.params = list(fill = 'lightgray'))+
                            stat_cor(method = 'pearson')+
                            labs(x = unique(tmp$celltype1), y = unique(tmp$celltype2))
        plist[[paste0(i,'_',j)]] <- p
    }
}

p <- wrap_plots(plist, ncol = 12)
png('/root/wangje/Project/刘老师/Myeloids/Fig/计算Myeloids和NKT的相关性_Rpost_NRpost.png',height=10000,width=12000,res=300)
p
dev.off()

######################################################################
############## 计算两者的相关性，保留增值细胞，去除p值大于0.05的组合
NKT_Sub <- NK_T_Fraction %>% select(Var1,celltype,Fraction,group) %>% 
    dplyr::rename_with(~'celltype1',2) %>% 
    filter(!Var1  %in% c("Three_CID3838--TNNaive","Three_CID3921--TNNaive","Three_CID45171--TNNaive"))

Myeloid_Sub <- Myeloids_T_Fraction %>% select(Var1,celltype,Fraction,group) %>% 
    dplyr::rename_with(~c("celltype2","Fraction2"),c(2,3))
Myeloid_Sub$Var1 <- gsub('CRPR-|SDPD-','',Myeloid_Sub$Var1)
merge_data <- cbind(NKT_Sub,Myeloid_Sub)

plist <- list()
for (i in unique(NKT_Sub$celltype1)) {
    for (j in unique(Myeloid_Sub$celltype2)) {
        tmp1 <- NKT_Sub %>% filter(celltype1 == i, group %in% c("R_Post", "NR_Post"))
        tmp2 <- Myeloid_Sub %>% filter(celltype2 == j, group %in% c("R_Post", "NR_Post"))
        tmp <- full_join(tmp1, tmp2, by = c("Var1", "group"))
        print(head(tmp))
        if ((cor.test(tmp$Fraction, tmp$Fraction2))$p.value <= 0.05) {
            p <- ggpubr::ggscatter(tmp,
                x = "Fraction", y = "Fraction2",
                add = "reg.line", conf.int = TRUE,
                add.params = list(fill = "lightgray")
            ) +
                stat_cor(method = "pearson") +
                labs(x = unique(tmp$celltype1), y = unique(tmp$celltype2))+
                theme(axis.title = element_text(size = 20,color = 'black'),
                    axis.text = element_text(size = 15,color = 'black'))
            plist[[paste0(i, "_", j)]] <- p
        }
    }
}

p <- wrap_plots(plist, ncol = 4)
png('/root/wangje/Project/刘老师/Myeloids/Fig/计算Myeloids和NKT的相关性_Rpost_NRpost_只保留p小于0.05的.png',height=4000,width=10000,res=300)
p & labs(title = 'NK & T cells  vs Myeloids R_Post and NR_Post') & theme(plot.title = element_text(size = 20 ,face = 'bold', hjust = 0.5))
dev.off()

#################################################################
###########计算R_Pre 和 NR_Pre之间的相关性，保留p-value小于0.05
plist <- list()
for (i in unique(NKT_Sub$celltype1)) {
    for (j in unique(Myeloid_Sub$celltype2)) {
        tmp1 <- NKT_Sub %>% filter(celltype1 == i, group %in% c("R_Pre", "NR_Pre"))
        tmp2 <- Myeloid_Sub %>% filter(celltype2 == j, group %in% c("R_Pre", "NR_Pre"))
        tmp <- full_join(tmp1, tmp2, by = c("Var1", "group"))
        print(head(tmp))
        if ((cor.test(tmp$Fraction, tmp$Fraction2))$p.value <= 0.05) {
            p <- ggpubr::ggscatter(tmp,
                x = "Fraction", y = "Fraction2",
                add = "reg.line", conf.int = TRUE,
                add.params = list(fill = "lightgray")
            ) +
                stat_cor(method = "pearson") +
                labs(x = unique(tmp$celltype1), y = unique(tmp$celltype2))+
                theme(axis.title = element_text(size = 20,color = 'black'),
                    axis.text = element_text(size = 15,color = 'black'))
            plist[[paste0(i, "_", j)]] <- p
        }
    }
}

p <- wrap_plots(plist, ncol = 4)
png('/root/wangje/Project/刘老师/Myeloids/Fig/计算Myeloids和NKT的相关性_Rpre_NRpre_只保留p小于0.05的.png',height=8000,width=10000,res=300)
p & labs(title = 'NK & T cells  vs Myeloids R_Pre and NR_Pre') & theme(plot.title = element_text(size = 20 ,face = 'bold', hjust = 0.5))
dev.off()


##################################################################
###########计算两者的相关性，去除增值细胞
Myeloid_Sub <- Myeloid_Sub %>% filter(celltype2 != 'Proliferating')
NKT_Sub <- NKT_Sub %>% filter(celltype1 != 'Proliferating cells')

plist <- list()
for(i in unique(NKT_Sub$celltype1)){
    for(j in unique(Myeloid_Sub$celltype2)){
        tmp1 <- NKT_Sub %>% filter(celltype1 == i,group %in% c('R_Post','NR_Post'))
        tmp2 <- Myeloid_Sub %>% filter(celltype2 == j, group %in% c('R_Post','NR_Post'))
        tmp <- full_join(tmp1,tmp2,by=c('Var1','group'))
        print(head(tmp))
        p <- ggpubr::ggscatter(tmp,x='Fraction',y='Fraction2',
                            add = 'reg.line',conf.int = TRUE,
                            add.params = list(fill = 'lightgray'))+
                            stat_cor(method = 'pearson')+
                            labs(x = unique(tmp$celltype1), y = unique(tmp$celltype2))+
                            theme(axis.title = element_text(size = 20),plot.title = element_text(size = 20))
        plist[[paste0(i,'_',j)]] <- p
    }
}
png('/root/wangje/Project/刘老师/Myeloids/Fig/计算Myeloids和NKT的相关性_Rpost_NRpost.png',height=10000,width=11000,res=300)
wrap_plots(plist,ncol=11)
dev.off()

####################################################################################
############ 所有细胞的相关性

