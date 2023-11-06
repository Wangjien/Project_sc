```R
library(Seurat)
library(qs)
library(ggplot2)
library(ggpubr)
library(patchwork)

## 按照文章中的条件进行复现
test <- scRNA@meta.data %>% select(cellSubType,expansion,clonotype,expansion,SAMPLE_ID) %>% na.omit

 test1 <- test %>%
      select(expansion,clonotype,SAMPLE_ID) %>%
      na.omit() %>%
      group_by(expansion, SAMPLE_ID) %>%
      distinct(clonotype) %>% summarise(n1 = n())

test2 <-  test %>% 
  select(expansion,SAMPLE_ID) %>%
  na.omit() %>% 
  group_by(expansion, SAMPLE_ID) %>% summarise(n = n())

test3 <- left_join(test1, test2, by =c('expansion','SAMPLE_ID'))
test3$Prop <- test3$n1/test3$n * 100

test3 <- test3 %>% dplyr::filter(expansion %in% c("E_On" ,"E_Pre" ,"NE_On" ,"NE_Pre"))

test3$expansion <- factor(test3$expansion, levels=c("NE_Pre",'E_Pre','NR_On','E_On'))
p1 <- ggplot(test3, aes(x = expansion, y = Prop,color = expansion,fill = expansion))+
    geom_boxplot(
        position = position_dodge(width = 0.8),
        outlier.shape = NA,color = 'black')+
    geom_jitter(aes(x = expansion, y = Prop,fill = expansion),
                color = 'black',position = position_dodge(width = 0.8),pch = 21)+
  ggpubr::stat_compare_means(list(c('ER_On','NE_On'),c('E_Pre','NE_Pre')))+
  labs(x = '', y = '% TCR richness')+
   scale_fill_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#ff7f01','#fb9a99'))+
    theme_test(base_size = 22,base_line_size = 1, base_rect_size = 1)+
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 18,color = 'black'),
          axis.text.y = element_text(size = 22, colour = 'black'),
          axis.title = element_text(size = 25, colour = 'black'),
          axis.text.x = element_text(hjust = 1,vjust = 0.5, colour = "black",size = 22,angle = 90))
## 以TRB的氨基酸序列作为指标进行分析
test <- scRNA@meta.data %>% select(cellSubType,expansion,TRBaa,expansion,SAMPLE_ID) %>% na.omit
 test1 <- test %>%
      select(expansion,TRBaa,SAMPLE_ID) %>%
      na.omit() %>%
      group_by(expansion, SAMPLE_ID) %>%
      distinct(TRBaa) %>% summarise(n1 = n())

test2 <-  test %>% 
  select(expansion,SAMPLE_ID) %>%
  na.omit() %>% 
  group_by(expansion, SAMPLE_ID) %>% summarise(n = n())

test3 <- left_join(test1, test2, by =c('expansion','SAMPLE_ID'))
test3$Prop <- test3$n1/test3$n * 100
test3 <- test3 %>% dplyr::filter(expansion %in% c("E_On" ,"E_Pre" ,"NE_On" ,"NE_Pre"))
# plot
test3$expansion <- factor(test3$expansion, levels=c("NE_Pre",'E_Pre','NE_On','E_On'))
p2 <- ggplot(test3, aes(x = expansion, y = Prop,color = expansion,fill = expansion))+
    geom_boxplot(
        position = position_dodge(width = 0.8),
        outlier.shape = NA,color = 'black')+
    geom_jitter(aes(x = expansion, y = Prop,fill = expansion),
                color = 'black',position = position_dodge(width = 0.8),pch = 21)+
  ggpubr::stat_compare_means(comparisons = list(c('E_On','NE_On'),c('E_Pre','NE_Pre')))+
  labs(x = '', y = '% TCR richness',title= 'PD1文章TRB氨基酸序列')+
   scale_fill_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#ff7f01','#fb9a99'))+
    theme_test(base_size = 22,base_line_size = 1, base_rect_size = 1)+
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 18,color = 'black'),
          axis.text.y = element_text(size = 22, colour = 'black'),
          axis.title = element_text(size = 25, colour = 'black'),
          axis.text.x = element_text(hjust = 1,vjust = 0.5, colour = "black",size = 22,angle = 90))

p1|p2
```
[![pilZmbF.png](https://z1.ax1x.com/2023/11/06/pilZmbF.png)](https://imgse.com/i/pilZmbF)
