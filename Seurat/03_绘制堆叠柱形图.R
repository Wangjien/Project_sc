library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(tidyverse)

# ---------------------------------------------------------------->>>
# 计算比例
# ---------------------------------------------------------------->>>
df.fraction <- as.data.frame(t(prop.table(table(scRNA_seurat$celltype,scRNA_seurat$new_Rename),margin=2))) %>% 
        dplyr::rename_with(~c("sample","celltype","fraction"),1:3)

# > head(df.fraction)
#    sample     celltype   fraction
# 1  P1-Pre NK & T cells 0.03229548
# 2 P1-Post NK & T cells 0.01845049
# 3  P2-Pre NK & T cells 0.47269400
# 4 P2-Post NK & T cells 0.41630263
# 5  P3-Pre NK & T cells 0.05994509
# 6 P3-Post NK & T cells 0.28127119

# ---------------------------------------------------------------->>>
# 绘图
# ---------------------------------------------------------------->>>
p1 <- ggplot(df.fraction, aes(sample, fraction, fill = celltype)) +
    geom_bar(stat = "identity", position = "fill", col = "black", ) +
    ggtitle("") +
    theme_classic() +
    theme(axis.ticks.length = unit(0.1, "cm")) +
    guides(fill = guide_legend(title = NULL)) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = "black", size = 18),
        axis.text.y = element_text(size = 15),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 18), axis.title.y = element_text(size = 24)
    ) +
    labs(x = " ", y = "Proportion") +
    guides(colour = guide_legend(override.aes = list(size = 10))) +
    scale_fill_manual(values = excel_color) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.key.height = unit(1.5, "line"))
