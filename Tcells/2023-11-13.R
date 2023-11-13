##--------------------------------------------------------------------
## TCR 数据分析
##--------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(qs)
library(patchwork)

##########################################################################################
#                     统计每种克隆的细胞数和克隆数，绘制点图
##########################################################################################
# 统计每种分组中每种TCR的Freqence
df_sub <- qread('C:/Users/Administrator/Desktop/IMG/TCR/NKT_metaData_02.qs')
result1 <- df_sub %>% dplyr::select(Treat_assess, TRBaa) %>%
    na.omit() %>%
  dplyr::group_by(Treat_assess,TRBaa) %>%
  summarise(n = n())
result1 <- table(result1$n, result1$Treat_assess) %>% as.data.frame()

library(forcats)
library(varhandle)
result1$Var1 <- unfactor(result1$Var1)

p6 <- ggplot(data = result1, aes(x = log2(Var1), y = log2(Freq+1),fill = Var2, shape = Var2))+
  geom_jitter(width = 0.1,size = 4,alpha = 0.8)+
  # geom_line()+
  scale_shape_manual(values = c(21,21,24,24))+
  scale_fill_manual(values = c('#99c7e5','#1e71ae','#b2df8a','#ff7f01','#fb9a99','#6a3d9a','#fbfb98','#b15928'))+
  scale_x_continuous(breaks = 0:15)+
  scale_y_continuous(breaks = 0:13)+
  geom_vline(xintercept = 1, lty="dashed", color = "red", linewidth = 0.5)+
  labs(x = 'log2(number of cells per clonotype)', y = 'log2(number of clonetypes)',fill = '', shape = '')+
  theme_test()+
  theme(axis.text = element_text(colour = 'black',size = 20),
        axis.title = element_text(size = 22, color = 'black'),
        legend.text = element_text(size = 18),
        legend.position = c(0.8,0.8),
        legend.key = element_rect(colour = '#dce1d8'),
        legend.spacing.x = unit(0.25, 'cm'),
        # legend.background = element_rect(color = 'black',linewidth = 0.8)
  )+
  guides(color=guide_legend(override.aes =list(size=5,alpha=1), shape = guide_legend(override.aes =list(size=5,alpha=1))))

pdf("C:/Users/Administrator/Desktop/IMG/TCR/New_TCR/T细胞TCR比例散点.pdf", height = 6,width = 8,family = 'ArialMT')
p6
dev.off()

##########################################################################################
#                     分别统计每种celltype的TCR克隆类型数，绘制柱形图
##########################################################################################
# 克隆型同届为 CloneSize =1 ，CloneSize = 2， 2< CloneSize <=5 和 CloneSize >5
# 分开CD4细胞和CD8细胞进行绘制

##--------------------------------------
##      CD8+ T cells
##--------------------------------------
# 设置clonetype类型
scRNA@meta.data <- scRNA@meta.data %>% mutate(
  cloneType = case_when(
    TRB_count == 1 ~ 'CloneSize = 1',
    TRB_count == 2 ~ 'CloneSize = 2',
    TRB_count > 2 & TRB_count <= 5 ~ '2 < CloneSize <=5',
    TRB_count > 5 ~ 'CloneSize > 5'
  )
)

# 筛选出CD4+ T cells和CD8+ T cells
cell1s <- c("CD8+ Tn","Texp","Tex","Tcyto")
cell2s <- c("CD4+ Tn","CD4+ Tcxcl13-ifng","Tfh","CD4+ Tem","SYNE2+ CD4+ T","Treg")

# 提取CD8 T cells进行统计
cd8_df <- scRNA@meta.data %>% dplyr::filter(celltype2 %in% cell1s) %>% dplyr::select(celltype2, Treat_assess, cloneType)
# 设置坐标轴
cd8_df <- cd8_df %>% mutate(
  x = case_when(
    celltype2 == 'CD8+ Tn' ~ 1,
    celltype2 == 'Texp' ~ 2,
    celltype2 == 'Tex' ~ 3,
    celltype2 == 'Tcyto' ~ 4
  )
)
test1 <- na.omit(cd8_df)
test1$cloneType <- factor(test1$cloneType, levels = c(
  rev(c("CloneSize = 1","CloneSize = 2","2 < CloneSize <=5","CloneSize > 5"))))
test1$group <- test1$cloneType

# 分组计算
test2 <- test1 %>% group_by(celltype2,Treat_assess,group,x) %>% summarise(n = n())
test1_R_Pre <- test2 %>% dplyr::filter(Treat_assess == 'R_Pre')
test1_R_Post <- test2 %>% dplyr::filter(Treat_assess == 'R_Post')
test1_NR_Pre <- test2 %>% dplyr::filter(Treat_assess == 'NR_Pre')
test1_NR_Post <- test2 %>% dplyr::filter(Treat_assess == 'NR_Post')
# 绘图
p3 <- ggplot() +
  geom_bar(data = test1_R_Pre,aes(x = x, y = n,fill = group),stat = 'identity',position = 'stack',width = 0.2)+
  geom_bar(data = test1_R_Post,aes(x = x + 0.22, y = n,fill = group),stat = 'identity',position = 'stack',width = 0.2)+
  geom_bar(data = test1_NR_Pre,aes(x = x + 0.44, y = n,fill = group),stat = 'identity',position = 'stack',width = 0.2)+
  geom_bar(data = test1_NR_Post,aes(x = x + 0.66, y = n,fill = group),stat = 'identity',position = 'stack',width = 0.2)+
  scale_fill_manual(values = rev(c('#673ca7','#b294c7','#fe992e','#fbfc94')))+
  theme_test(base_line_size = 1,base_rect_size = 1,base_size = 40)+
  scale_x_continuous(expand = c(0.01,0),breaks = c(1.3,2.3,3.3,4.3),
                     labels = c('CD8+ Tn','Texp','Tex','Tcyto'))+
  scale_y_continuous(expand = c(0.02,0))+
  labs(fill = '', y = 'Frequence')+
  theme(axis.text.x = element_text(colour = 'black',angle = 90,hjust = 1,vjust = 0.8),
        axis.text.y = element_text(colour = 'black'),
        axis.title.x = element_blank())+
  geom_vline(xintercept = 1.84, lty="dashed", color = "grey50", linewidth = 0.8)+
  geom_vline(xintercept = 2.84, lty="dashed", color = "grey50", linewidth = 0.8)+
  geom_vline(xintercept = 3.84, lty="dashed", color = "grey50", linewidth = 0.8)

pdf("/root/wangje/Project/刘老师/NK_T/Data/CCA新的分群/new_CD8/cd8_T细胞celltype堆叠柱形图_2.pdf", height = 5,width = 12,family = 'ArialMT')
p3
dev.off()

##--------------------------------------
##      CD4+ T cells
##--------------------------------------
# CD4 T cells
cd4_df <- scRNA@meta.data %>% dplyr::filter(celltype2 %in% cell2s) %>% select(celltype2, Treat_assess, cloneType)

# 绘图
cd4_df <- cd4_df %>% mutate(
  x = case_when(
    celltype2 == 'CD4+ Tn' ~ 1,
    celltype2 == 'CD4+ Tcxcl13-ifng' ~ 2,
    celltype2 == 'Tfh' ~ 3,
    celltype2 == 'CD4+ Tem' ~ 4,
    celltype2 == 'SYNE2+ CD4+ T' ~ 5,
    celltype2 == 'Treg' ~ 6
  )
)

# 分别提取数据
test1 <- na.omit(cd4_df)
test1$cloneType <- factor(test1$cloneType, levels = c(
  rev(c("CloneSize = 1","CloneSize = 2","2 < CloneSize <=5","CloneSize > 5")
)))
test1$group <- test1$cloneType

## 分组计算
test2 <- test1 %>% group_by(celltype2,Treat_assess,group,x) %>% summarise(n = n())


test1_R_Pre <- test2 %>% dplyr::filter(Treat_assess == 'R_Pre')
test1_R_Post <- test2 %>% dplyr::filter(Treat_assess == 'R_Post')
test1_NR_Pre <- test2 %>% dplyr::filter(Treat_assess == 'NR_Pre')
test1_NR_Post <- test2 %>% dplyr::filter(Treat_assess == 'NR_Post')

# 绘制柱形图
p4 <- ggplot() +
  geom_bar(data = test1_R_Pre,aes(x = x, y = n,fill = group),stat = 'identity',position = 'stack',width = 0.2)+
  geom_bar(data = test1_R_Post,aes(x = x + 0.22, y = n,fill = group),stat = 'identity',position = 'stack',width = 0.2)+
  geom_bar(data = test1_NR_Pre,aes(x = x + 0.44, y = n,fill = group),stat = 'identity',position = 'stack',width = 0.2)+
  geom_bar(data = test1_NR_Post,aes(x = x + 0.66, y = n,fill = group),stat = 'identity',position = 'stack',width = 0.2)+
  scale_fill_manual(values = rev(c('#673ca7','#b294c7','#fe992e','#fbfc94')))+
  theme_test(base_line_size = 1,base_rect_size = 1,base_size = 40)+
  scale_x_continuous(expand = c(0.01,0),breaks = c(1.3,2.3,3.3,4.3,5.3,6.3),
                     labels = c('CD4+ Tn','CD4+ Tcxcl13-ifng','Tfh','CD4+ Tem','SYNE2+ CD4+ T','Treg'))+
  scale_y_continuous(expand = c(0.02,0))+
  labs(fill = '', y = 'Frequence')+
  theme(axis.text.x = element_text(colour = 'black',angle = 90,hjust = 1,vjust = 0.8),
        axis.text.y = element_text(colour = 'black'),
        axis.title.x = element_blank())+
  geom_vline(xintercept = 1.84, lty="dashed", color = "grey50", linewidth = 0.8)+
  geom_vline(xintercept = 2.84, lty="dashed", color = "grey50", linewidth = 0.8)+
  geom_vline(xintercept = 3.84, lty="dashed", color = "grey50", linewidth = 0.8)+
  geom_vline(xintercept = 4.84, lty="dashed", color = "grey50", linewidth = 0.8)+
  geom_vline(xintercept = 5.84, lty="dashed", color = "grey50", linewidth = 0.8)

pdf("/root/wangje/Project/刘老师/NK_T/Data/CCA新的分群/new_CD4/cd4_T细胞celltype堆叠柱形图.pdf", height = 8,width = 18,family = 'ArialMT')
p4
dev.off()

##########################################################################################
#                    计算不同分组的TCR richness， 绘制箱线图
##########################################################################################
# 要求同时区分出CD4和CD8 Tcells，以及不同的celltype
library(qs)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
# 去除NK cells中的数据
scRNA <- scRNA[,!scRNA$celltype2 %in% c("CD16-CD56bright NK","CD16-CD56dim NK")]
# 区分出CD4或CD8 T cells
scRNA@meta.data <- scRNA@meta.data %>% mutate(
    celltype3 = case_when(
        celltype2 %in% c('"CD8+ Tn','Tex','Texp','Tcyto') ~ 'CD8+ T',
        TRUE ~ 'CD4+ T'
    )
)
# 分組进行统计
test1 <- scRNA@meta.data %>% 
  dplyr::select(Treat_assess,celltype3,Patient,TRBaa) %>% 
  na.omit() %>% 
  group_by(Treat_assess,celltype3, Patient) %>%  
  summarise(TRB_type = unique(TRBaa)) %>% summarise(n1 = n())
head(test1)
# # A tibble: 6 × 4
# # Groups:   Treat_assess, celltype3 [1]
#   Treat_assess celltype3 Patient    n1
#   <fct>        <chr>     <chr>   <int>
# 1 R_Pre        CD4+ T    CJME      151
# 2 R_Pre        CD4+ T    CMDI      664
# 3 R_Pre        CD4+ T    HDYU      178
# 4 R_Pre        CD4+ T    HXZH      699
# 5 R_Pre        CD4+ T    LCLI     1629
# 6 R_Pre        CD4+ T    WYJU       82  
test2 <-  scRNA@meta.data %>% 
  dplyr::select(Treat_assess,celltype3,Patient) %>%
  na.omit() %>% 
  group_by(Treat_assess,celltype3, Patient) %>% summarise(n = n())
# 合并数据
test3 <- left_join(test1, test2, by =c('Treat_assess','celltype3' ,'Patient'))
test3$Prop <- test3$n1/test3$n * 100
head(test3)
# # A tibble: 6 × 6
# # Groups:   Treat_assess, celltype3 [1]
#   Treat_assess celltype3 Patient    n1     n  Prop
#   <fct>        <chr>     <chr>   <int> <int> <dbl>
# 1 R_Pre        CD4+ T    CJME      151   316  47.8
# 2 R_Pre        CD4+ T    CMDI      664  3728  17.8
# 3 R_Pre        CD4+ T    HDYU      178   480  37.1
# 4 R_Pre        CD4+ T    HXZH      699  1621  43.1
# 5 R_Pre        CD4+ T    LCLI     1629  5083  32.0
# 6 R_Pre        CD4+ T    WYJU       82   283  29.0

# 计算显著性
stat.test <- test3 %>%
    group_by(celltype3) %>%
    wilcox_test(Prop ~ Treat_assess,comparisons = list(c('R_Pre','R_Post'),c('NR_Pre','NR_Post'),c('R_Pre','NR_Pre'),c('R_Post','NR_Post')))
stat.test <- stat.test %>%
    add_xy_position(x = "celltype3",dodge = 0.8)
stat.test
# head(stat.test)
# # A tibble: 6 × 15
#   celltype3 .y.   group1 group2     n1    n2 statistic     p p.adj
#   <chr>     <chr> <chr>  <chr>   <int> <int>     <dbl> <dbl> <dbl>
# 1 CD4+ T    Prop  R_Pre  R_Post      9     8        36 1     1    
# 2 CD4+ T    Prop  NR_Pre NR_Post     8     5         3 0.011 0.044
# 3 CD4+ T    Prop  R_Pre  NR_Pre      9     8        53 0.114 0.342
# 4 CD4+ T    Prop  R_Post NR_Post     8     5        12 0.284 0.568
# 5 CD8+ T    Prop  R_Pre  R_Post      9     8        37 0.963 1    
# 6 CD8+ T    Prop  NR_Pre NR_Post     8     5         7 0.065 0.237

# 绘图
p1 <- ggplot(test3)+
  geom_boxplot(aes(x = celltype3, y = Prop, fill = Treat_assess),outlier.shape = NA,position = position_dodge(width = 0.8))+
  geom_jitter(aes(x = celltype3, y = Prop, fill = Treat_assess),color = 'black',position = position_dodge(width = 0.8),pch = 21)+
  stat_pvalue_manual(
    stat.test, label = "p", tip.length = 0.00,size = 4,
      hide.ns = FALSE
  )+
  labs(y = '% TCR richness',x = '')+
  scale_fill_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#ff7f01','#fb9a99'))+
  theme_test(base_size = 22,base_line_size = 1, base_rect_size = 1)+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 18,color = 'black'),
        axis.text.y = element_text(size = 22, colour = 'black'),
        axis.title = element_text(size = 25, colour = 'black'),
        axis.text.x = element_text(hjust = 1,vjust = 0.5, colour = "black",size = 22,angle = 90))

##########################################################################################
#                    计算不同分组的TCR Clonality， 绘制箱线图
##########################################################################################
# Clonality: Clonality = 1 – Evenness = 1 - Diversity/log(Richness)









##########################################################################################
#                    Reactome GSEA
##########################################################################################
########### clusterProlifer; Reactome  gsePathway()













########### hypeR



