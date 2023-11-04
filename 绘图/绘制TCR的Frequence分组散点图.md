## 绘制TCR分组Frequence散点图
```R
library(dplyr)
library(ggplot2)
library(qs)

# 读入文件
# 统计每种分组中每种TCR的Freqence
df_sub <- qread('C:/Users/Administrator/Desktop/IMG/TCR/NKT_metaData_02.qs')
result1 <- df_sub %>% dplyr::group_by(Treat_assess,TRBaa) %>% summarise(n = n())
result1 <- table(result1$n, result1$Treat_assess) %>% as.data.frame()

# 去除factor
library(forcats)
library(varhandle)
result1$Var1 <- unfactor(result1$Var1)

# 绘制图片
p6 <- ggplot(data = result1, aes(x = log2(Var1), y = log2(Freq+1),color = Var2, shape = Var2))+
  geom_jitter(width = 0.1,size = 3,alpha = 0.6)+
  scale_x_continuous(breaks = 0:15)+
  scale_y_continuous(breaks = 0:13)+
  geom_vline(xintercept = 1, lty="dashed", color = "red", linewidth = 0.5)+
  labs(x = 'log2(number of cells per clonotype)', y = 'log2(number of clonetypes)',color = '', shape = '')+
  scale_shape_manual(values = c(16,16,17,17))+
  scale_color_manual(values = c('#99c7e5','#1e71ae','#b2df8a','#ff7f01','#fb9a99','#6a3d9a','#fbfb98','#b15928'))+
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

# 保存图片
pdf("C:/Users/Administrator/Desktop/IMG/TCR/New_TCR/T细胞TCR比例散点.pdf", height = 6,width = 8,family = 'ArialMT')
p6
dev.off()
```
[![piQCp6J.png](https://z1.ax1x.com/2023/11/04/piQCp6J.png)](https://imgse.com/i/piQCp6J)
