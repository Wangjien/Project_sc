library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(VISION)
library(purrr)
library(tibble)
library(gghalves)
library(introdataviz)

# 读入数据
# -----------> Seurat 数据
setwd("/root/wangje/Project/刘老师/Myeloids/Vision_result/")
load("/root/wangje/Project/刘老师/Myeloids/Data/Myeloids_CCA.RData")
scRNA_seurat$Treat_assess <- gsub("CRPR-|SDPD-", "", scRNA_seurat$Treat_assess)
splitBy <- c("R_Pre", "R_Post", "NR_Pre", "NR_Post")

# ------------> Angiogenesis
# ! Angiogenesis
signatures <- c("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/Angiogenesis.gmt")

# signatures = c("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/Toll_like_receptor_signaling_pathway.gmt")
# singnaturs <- c("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/Antigen_Proc_and_Pres_via_MHC.gmt")

# -----------> 细胞凋亡
signatures <- c("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/细胞凋亡.gmt")

# -----------> IFN gamma
signatures <- c("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/IFN_GAMA.gmt")

# -----------> Antigen_Proc_and_Pres_via_MHC
signatures <- c("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/Antigen_Proc_and_Pres_via_MHC.gmt")

# -----------> Cytokine_cytokine_receptor
signatures <- c("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/Cytokine_cytokine_receptor.gmt")

object <- scRNA_seurat
flist <- list()
for (i in splitBy) {
    tmp <- object[, object$Treat_assess %in% i]
    print(dim(tmp))
    # 计算VISION
    Score <- Vision(tmp,
        signatures = signatures,
        projection_methods = NULL
    ) %>% analyze()
    flist[[i]] <- Score@SigScores
}

# 获取Myeloids的meta.data
load("/root/wangje/Project/刘老师/Myeloids/Data/Myeloids_CCA.RData")
meta.data <- scRNA_seurat@meta.data %>%
    select(c("orig.ident", "Patient", "Rename", "Treat_assess", "Tissue", "Rseponse", "celltype")) %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(Treat_assess = gsub("CRPR-|SDPD-", "", .$Treat_assess)) %>%
    write.table(file = "/root/wangje/Project/刘老师/Myeloids/Vision_result/Myeloids_metaData.txt", sep = "\t", row.names = F, quote = F)

# 合并flist
test <- purrr::map_df(
    flist, function(df) {
        as.data.frame(df) %>%
            tibble::rownames_to_column("cell") %>%
            dplyr::rename_with(~"values", 2)
    }
)

# 添加其他的信息
meta.data <- read.table("/root/wangje/Project/刘老师/Myeloids/Vision_result/Myeloids_metaData.txt", sep = "\t", header = T)
test <- test %>%
    dplyr::mutate(
        Patient = meta.data[match(test$cell, meta.data$cell), "Patient"],
        Rename = meta.data[match(test$cell, meta.data$cell), "Rename"],
        Treat_assess = meta.data[match(test$cell, meta.data$cell), "Treat_assess"],
        celltype = meta.data[match(test$cell, meta.data$cell), "celltype"]
    ) %>%
    dplyr::mutate(
        group1 = stringr::str_split_fixed(Treat_assess, "_", n = 2)[, 1],
        group2 = stringr::str_split_fixed(Treat_assess, "_", n = 2)[, 2]
    )

# 细胞凋亡
Apoptosis <- test
Apoptosis$Group <- "Cell Apoptosis"
write.table(Apoptosis, file = "/root/wangje/Project/刘老师/Myeloids/Vision_result/细胞凋亡.txt", row.names = F, sep = "\t", quote = F)

# 血管生成
Angiogenesis <- test
Angiogenesis$Group <- "Angiogenesis"
write.table(Angiogenesis, file = "/root/wangje/Project/刘老师/Myeloids/Vision_result/血管生成.txt", row.names = F, sep = "\t", quote = F)

# IFN Gamma
IFNGamma <- test
IFNGamma$Group <- "IFN Gamma"
write.table(IFNGamma, file = "/root/wangje/Project/刘老师/Myeloids/Vision_result/IFNGamma.txt", row.names = F, sep = "\t", quote = F)

# Cytokine
Cytokine <- test
Cytokine$Group <- "Cytokine"
write.table(Cytokine, file = "/root/wangje/Project/刘老师/Myeloids/Vision_result/Cytokine.txt", row.names = F, sep = "\t", quote = F)

# Antigen_Process
Antigen_Process <- test
Antigen_Process$Group <- "Antigen_Process"
write.table(Antigen_Process, file = "/root/wangje/Project/刘老师/Myeloids/Vision_result/Antigen_Process.txt", row.names = F, sep = "\t", quote = F)

# 绘图
df1 <- read.table("细胞凋亡.txt", header = T, sep = "\t", fill = T)
df2 <- read.table("血管生成.txt", header = T, sep = "\t", fill = T)
df3 <- read.table("Antigen_Process.txt", header = T, sep = "\t", fill = T)
df4 <- read.table("Cytokine.txt", header = T, sep = "\t", fill = T)
df5 <- read.table("IFNGamma.txt", header = T, sep = "\t", fill = T)
merge_df <- rbind(df1, df2, df3, df4, df5)

plist1 <- list()
for (celltype in unique(merge_df$Group)) {
    tmp <- merge_df[merge_df$Group == celltype, ]
    p <- ggplot(tmp, aes(x = group1, y = values, fill = group2)) +
        geom_split_violin(trim = TRUE) +
        geom_half_boxplot(
            data = tmp %>% filter(group2 == "Pre"), aes(x = group1, y = values),
            width = 0.15, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef = 0, color = "black", linetype = 1, fill = "white"
        ) +
        geom_half_boxplot(
            data = tmp %>% filter(group2 == "Post"), aes(x = group1, y = values),
            width = 0.15, side = "r", notch = FALSE, notchwidth = .4, outlier.shape = NA, coef = 0, color = "black", linetype = 1, fill = "white"
        ) +
        labs(x = NULL, y = "Signature Score", fill = "Group", title = celltype) +
        theme_classic() +
        theme(
            text = element_text(size = 20), axis.title = element_text(size = 23, color = "black"),
            axis.text = element_text(size = 20, colour = "black"),
            plot.title = element_text(hjust = 0.5)
        ) +
        scale_fill_manual(values = c("#3b9aac", "#a74947", "#4470a1", "#8da24f", "#6d5a87")) +
        stat_compare_means(method = "wilcox.test")
    plist1[[celltype]] <- p
    print(celltype)
}

p <- wrap_plots(plist1, ncol = 5)
ggsave(filename = "/root/wangje/Project/刘老师/Myeloids/Fig/signaturePlot.png", height = 8, width = 20, plot = p)


p <- ggviolin(merge_data, "ident", "value",
    fill = "ident",
    outlier.shape = NA, trim = F
)

p <- p + stat_compare_means(comparisons = list(c("R_Pre", "R_Post"), c("R_Pre", "NR_Pre"), c("R_Post", "NR_Pre"), c("NR_Pre", "NR_Post"))) +
    theme(
        axis.title.x = element_blank(),
        axis.line = element_line(linewidth = 0.8, colour = "black"),
        legend.position = "bottom",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_text(size = 18, colour = "black"),
        legend.text = element_text(size = 15)
    ) +
    labs(x = "", y = "GLI1", fill = "Group")

# 绘制不同分组的小提琴图
plist2 <- list()
for (des in unique(merge_df$Group)) {
    tmp <- merge_df[merge_df$Group == des, ]
    p <- ggviolin(tmp, "Treat_assess", "values", fill = "Treat_assess", outlier.shape = NA, trim = F, add = "mean_sd", error.plot = "crossbar") +
        stat_compare_means(comparisons = list(c("R_Pre", "R_Post"), c("R_Pre", "NR_Pre"), c("R_Post", "NR_Pre"), c("NR_Pre", "NR_Post"))) +
        theme(
            axis.title.x = element_blank(),
            axis.line = element_line(linewidth = 0.8, colour = "black"),
            legend.position = "none",
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 25),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.title = element_text(size = 18, colour = "black"),
            legend.text = element_text(size = 15),
            plot.title = element_text(hjust = 0.5)
        ) +
        labs(x = "", y = "Signature Score", fill = "Group", title = des)
    plist2[[des]] <- p
    print(des)
}

p1 <- wrap_plots(plist2, ncol = 5)
ggsave(filename = "/root/wangje/Project/刘老师/Myeloids/Fig/signaturePlot2.png", height = 8, width = 20, plot = p1)

p2 <- p / p1
ggsave(filename = "/root/wangje/Project/刘老师/Myeloids/Fig/signaturePlot3.png", height = 10, width = 30, plot = p2)
