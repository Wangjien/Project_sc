## 绘制每种celltype的TCR堆叠柱形图
### clonetypes类型分为了以下几种类型
'CloneSize = 1','CloneSize = 2','2 < CloneSize <=5','CloneSize > 5'
```R
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(qs)
# 读入数据
setwd('/root/wangje/Project/刘老师/NK_T/Data/CCA新的分群')
scRNA <- qread('new_CCA新的分群结果.qs')

# 数据中添加clone types
scRNA@meta.data <- scRNA@meta.data %>% mutate(
  cloneType = case_when(
    TRB_count == 1 ~ 'CloneSize = 1',
    TRB_count == 2 ~ 'CloneSize = 2',
    TRB_count > 2 & TRB_count <= 5 ~ '2 < CloneSize <=5',
    TRB_count > 5 ~ 'CloneSize > 5'
  )
)
# 筛选出CD4和CD8 T细胞
cell1s <- c("CD8+ Tn","Texp","Tex","Tcyto")
cell2s <- c("CD4+ Tn","CD4+ Tcxcl13-ifng","Tfh","CD4+ Tem","SYNE2+ CD4+ T","Treg","Cycling T")
# CD4 T cells
cd4_df <- scRNA@meta.data %>% dplyr::filter(celltype2 %in% cell1s) %>% select(celltype2, Treat_assess, cloneType)

# 绘图
cd4_df <- cd4_df %>% mutate(
  x = case_when(
    celltype2 == 'CD8+ Tn' ~ 1,
    celltype2 == 'Texp' ~ 2,
    celltype2 == 'Tex' ~ 3,
    celltype2 == 'Tcyto' ~ 4
  )
)
test1 <- na.omit(cd4_df)
test1$cloneType <- factor(test1$cloneType, levels = c(
  "CloneSize = 1","CloneSize = 2","2 < CloneSize <=5","CloneSize > 5"
))
test1$group <- test1$cloneType
```
[![piMx8bV.md.png](https://z1.ax1x.com/2023/11/04/piMx8bV.md.png)](https://imgse.com/i/piMx8bV)

```R
test2 <- test1 %>% group_by(celltype2,Treat_assess,group,x) %>% summarise(n = n())
head(test2)
```
[![piMx4KI.png](https://z1.ax1x.com/2023/11/04/piMx4KI.png)](https://imgse.com/i/piMx4KI)

```R
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
  scale_fill_manual(values = rev(c('#0d0887','#9c179e','#ed7953','#f0f921')))+
  theme_test(base_line_size = 1,base_rect_size = 1,base_size = 40)+
  scale_x_continuous(expand = c(0.01,0),breaks = c(1.3,2.3,3.3,4.3,5.3,6.3,7.3,8.3,9.3,10.3),
                     labels = c('CD4+ Tn','CD4+ Tcxcl13-ifng','Tfh','CD4+ Tem','Treg','SYNE2+ CD4+ T','CD8+ Tn',
                                'Tex','Texp','Cycling T'))+
  scale_y_continuous(expand = c(0.02,0))+
  labs(fill = '', y = 'Frequence')+
  theme(axis.text.x = element_text(colour = 'black',angle = 90,hjust = 1,vjust = 0.8),
        axis.text.y = element_text(colour = 'black'),
        axis.title.x = element_blank())+
  geom_vline(xintercept = 1.84, lty="dashed", color = "grey50", linewidth = 0.8)+
  geom_vline(xintercept = 2.84, lty="dashed", color = "grey50", linewidth = 0.8)+
  geom_vline(xintercept = 3.84, lty="dashed", color = "grey50", linewidth = 0.8)
p3
```
[![piMxxrq.png](https://z1.ax1x.com/2023/11/04/piMxxrq.png)](https://imgse.com/i/piMxxrq)

























