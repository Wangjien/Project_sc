# 绘制Split_violin_plot 结合boxplot

## 安装需要的R包

```R
install.packages('rstatix')
if (!require(devtools)) {
    install.packages('devtools')
}
devtools::install_github('erocoar/gghalves')
devtools::install_github("psyteachr/introdataviz")
```

## 查看数据

```R
library(tidyverse)
library(ggplot2)
library(patchwork)
library(gghalves)
library(stringr)
library(rstatix)

set.seed(123)
data <- tibble(x = rep(c("A", "B"), each = 100),
               y = rep(c("C", "D"), times = 100),
               z = c(rnorm(100), rnorm(100, mean = 1, sd = 0.5)))
head(data)
```

<img src="https://s1.ax1x.com/2023/03/02/ppFOAG8.png" style="zoom: 150%;" />

```R
ggplot(data, aes(x = x, y = z, fill = y)) +
    geom_split_violin(trim = TRUE) +
    geom_half_boxplot(
        data = data %>% filter(y == "C"), aes(x = x, y = z),
        width = 0.15, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef = 0, color = "black", linetype = 1, fill = "white"
    ) +
    geom_half_boxplot(
        data = data %>% filter(y == 'D'), aes(x = x, y = z),
        width = 0.15, side = "r", notch = FALSE, notchwidth = .4, outlier.shape = NA, coef = 0, color = "black", linetype = 1, fill = "white"
    ) +
    labs(x = NULL, y = "Signature score", fill = "Group", title = 'Split violin plot and boxplot') +
    theme_classic() +
    theme(
        text = element_text(size = 20), 
        axis.title = element_text(size = 23, color = "black"),
        axis.text = element_text(size = 20, colour = "black"),
        plot.title = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = c("#3b9aac", "#a74947", "#4470a1", "#8da24f", "#6d5a87")) +
    stat_compare_means(method = "wilcox.test")
```

![image-20230302183958679](https://s1.ax1x.com/2023/03/02/ppFOERS.png)

















