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
##      CD4+ T cells
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
##      CD8+ T cells
##--------------------------------------