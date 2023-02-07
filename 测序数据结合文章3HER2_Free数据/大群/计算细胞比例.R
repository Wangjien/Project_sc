# library
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(cowplot)
library(stringr)
# load data
setwd("/root/wangje/Project/刘老师/大群/")
load("/root/wangje/Project/刘老师/大群/大群CCA.Rdata")## ---------------- 大群计算比例，在原有基础上增加3Her2数据 ---------------
scRNA_sub <- scRNA_seurat
scRNA_sub$patients <- paste0(scRNA_sub$Patient,"--",scRNA_sub$Treat_assess,"--",scRNA_sub$Tissue)
Proportion_analys <- function(scRNA_seurat, IDs = "patients", celltype = "celltype") { # nolint
    ids <- table(scRNA_seurat@meta.data[, IDs]) %>% as.data.frame()
    df.merge <- table(scRNA_seurat@meta.data[, c(IDs, celltype)]) %>% as.data.frame()
    for (i in 1:dim(ids)[1]) {
        prop.temp <- df.merge[which(df.merge$patients %in% ids[i, 1]), "Freq"] / ids[i,2]
        df.merge[which(df.merge$patients %in% ids[i, 1]),"SampleCellNum"] <- ids[i,2]
        df.merge[which(df.merge$patients %in% ids[i, 1]), "prop"] <- prop.temp
        df.merge <<- df.merge
    }
}
# 计算刘老师大群的比例
scRNA_sub %>% Proportion_analys()
df.merge <- df.merge %>%
    as_tibble() %>%
    mutate(
        sample = str_split_fixed(df.merge$patients, "--", n = 3)[, 1],
        group = str_split_fixed(df.merge$patients, "--", n = 3)[, 2],
        tissue <- str_split_fixed(df.merge$patients, "--", n = 3)[, 3]
    )
Liu_big <- df.merge
df.merge <- NULL

# ---------------------- 计算3HER2的大群比例 -------------------------
# 读入大群数据
scRNA_seurat@meta.data[which(scRNA_seurat$integrated_snn_res.0.4 %in% c(8)), "celltype"] <- "B & Plasma cells"
scRNA_seurat@meta.data[which(scRNA_seurat$integrated_snn_res.0.4 %in% c(10)), "celltype"] <- "Endothelial cells"
scRNA_seurat@meta.data[which(scRNA_seurat$integrated_snn_res.0.4 %in% c(3, 5, 11)), "celltype"] <- "Epithelial cells"
scRNA_seurat@meta.data[which(scRNA_seurat$integrated_snn_res.0.4 %in% c(9)), "celltype"] <- "Fibroblasts"
scRNA_seurat@meta.data[which(scRNA_seurat$integrated_snn_res.0.4 %in% c(6, 12)), "celltype"] <- "Myeloid cells"
scRNA_seurat@meta.data[which(scRNA_seurat$integrated_snn_res.0.4 %in% c(0, 1, 2, 4, 7)), "celltype"] <- "NK & T cells"
save(scRNA_seurat,file="/root/wangje/文献数据/5HER2/Wu_etal_2021_BRCA_scRNASeq/Wu/重新分析3HER2/bigcluster_3her2_cca.rdata")
load('/root/wangje/文献数据/5HER2/Wu_etal_2021_BRCA_scRNASeq/Wu/重新分析3HER2/bigcluster_3her2_cca.rdata')
## 计算3HER2的比例
scRNA_seurat$Patient <- "3HER2"
scRNA_seurat$Treat_assess <- "TRNaive"
scRNA_seurat$Tissue <- 'Breast'
scRNA_seurat$Patient <- scRNA_seurat$sample
scRNA_seurat$patients <- paste0(scRNA_seurat$Patient,"--",scRNA_seurat$Treat_assess,"--",scRNA_seurat$Tissue)
scRNA_seurat %>% Proportion_analys()
df.merge <- df.merge %>%
    as_tibble() %>%
    mutate(
        sample = str_split_fixed(df.merge$patients, "--", n = 3)[, 1],
        group = str_split_fixed(df.merge$patients, "--", n = 3)[, 2],
        tissue <- str_split_fixed(df.merge$patients, "--", n = 3)[, 3]
    )
# 3HER2数据
write.table(df.merge,file="/root/wangje/Project/刘老师/大群/结合3HER2数据进行分析/3HER2_大群CCA以每个样本的数目计算比例.txt",quote=F,sep="\t",row.names = F)
# Liu 大群数据
write.table(Liu_big,file="/root/wangje/Project/刘老师/大群/结合3HER2数据进行分析/Liu_大群CCA以每个样本的数目计算比例.txt",quote = F,sep = "\t",row.names = F)
# 合并两个数据的比例
df.merge <- rbind(Liu_big,df.merge)
write.table(df.merge,file="/root/wangje/Project/刘老师/大群/结合3HER2数据进行分析/合并Liu和3HER2_大群CCA以每个样本的数目计算比例.txt",quote = F,sep = "\t",row.names = F)

# 绘图
library(ggpubr)
library(ggpllot2)
library(patchwork)
library(cowplot)
library(tidyverse)
library(ggrepel)

df.merge$Patient <- str_split_fixed(df.merge$patients,"_",n=2)[,1]
df.merge$Patient <- str_split_fixed(df.merge$Patient,"--",n=2)[,1]
df.merge$Patient <- gsub("FHXHBS2","FHXHBS1",df.merge$Patient)
df.merge$Patient <- gsub("HEJX","HEJI",df.merge$Patient)
df.merge$group <- factor(df.merge$group, levels = c("TRNaive", "R_Pre", "R_Post", "NR_Pre", "NR_Post"))
compaired <- list(c("TRNaive", "R_Pre"), c("TRNaive", "R_Post"), c("TRNaive", "NR_Pre"), c("TRNaive", "NR_Post"), c("R_Pre", "R_Post"), c("NR_Pre", "NR_Post"), c("R_Pre", "NR_Pre"), c("R_Post", "NR_Post"))
excel_color <- c("#3b9aac", "#a74947", "#4470a1", "#8da24f", "#6d5a87")
plist <- list()
for (celltype in unique(df.merge$celltype)) {
    tmp <- df.merge[df.merge$celltype == celltype, ]
    p <- ggplot(tmp, aes(x = group, y = prop)) +
        geom_boxplot(aes(fill = group), show.legend = F, width = 0.6, outlier.colour = NA) +
        geom_jitter(width = 0.1,size = 2, fill = "black",color="black") + # 绘制散点
        geom_line(aes(group = Patient), color = "darkgrey", lwd = 0.5) + # 配对样本间连线
        theme(
            panel.grid = element_blank(),
            legend.position = "none",
            axis.line = element_line(colour = "black", size = 1.2),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5,size=23,face="bold"),
            plot.subtitle = element_text(size = 23, hjust = 0.5),
            axis.text = element_text(size = 23, color = "black"),
            axis.title = element_text(size = 23, color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        ) +
        labs(x = "", y = "Fraction", title = celltype) +
        geom_text_repel(aes(group, prop, label=sample),size=2)+
        # ggpubr::stat_compare_means(method = "wilcox.test", hide.ns = F, comparisons = compaired, label = "p.format", symnum.args = symnum.args, size = 6) +
        ggpubr::stat_compare_means(method="wilcox.test",hide.ns = F,comparisons = compaired,label="p.formot")+
        scale_fill_manual(values = c("#3b9aac","#a74947","#4470a1","#8da24f","#6d5a87"))
    plist[[celltype]] <- p
}

        
png("/root/wangje/Project/刘老师/大群/结合3HER2数据进行分析/合并Liu和3HER2_大群CCA以每个样本的数目计算比例.png", height = 4500, width = 6000, res = 300)
p_freq <- wrap_plots(plist, bycol = T, ncol = 3)
p_freq
dev.off()
p_freq
