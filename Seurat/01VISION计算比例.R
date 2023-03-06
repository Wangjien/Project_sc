library(Seurat)
library(ggplot2)
library(patchwork)
library(gghalves)
library(dplyr)
library(VISION)
library(stringr)
library(ggunchained)
library(ggsci)


# 读入vision的计算数据
load("/root/wangje/Project/刘老师/EPi/新结果/添加肝脏/new/new_CCA_EPI.RData")

id1 <- c("CJME", "CMDI", "HDYU", "HXZH", "LCLI", "WYJU", "WZLA", "ZXME", "ZJLI_0116")
id2 <- c("CZYI", "FHXHBS1", "FYYI", "HEJI", "LAWE", "ZEYI", "ZFXI", "LIPE")
id3 <- c("CJME_0707", "CMDI_0624", "HDYU_0720", "HXZH_0220", "LCLI_0623", "WYJU_0122", "ZXME_0223", "ZJLI_0312")
scRNA_CCA[["Patient"]] <- unlist(strsplit(rownames(scRNA_CCA@meta.data), split = "_[A|T|C|G]*$"))
scRNA_CCA@meta.data$Treat_assess <- ifelse(scRNA_CCA@meta.data$Patient %in% id1, "R_Pre",
    ifelse(scRNA_CCA@meta.data$Patient %in% id2, "NR_Pre",
        ifelse(scRNA_CCA@meta.data$Patient %in% id3, "R_Post", "NR_Post")
    )
)

# ! ------------------------------------------------------------------------> IFN gama <--------------------------------------------------------------------------------------------------------------------------------

signatures <- ("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/IFN_GAMA.gmt")
splitBy <- c("R_Pre", "R_Post", "NR_Pre", "NR_Post")
# signatures = c("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/Toll_like_receptor_signaling_pathway.gmt")
# singnaturs <- c("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/Antigen_Proc_and_Pres_via_MHC.gmt")
object <- scRNA_CCA
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

## 添加其他信息
test <- Reduce(rbind, flist) %>%
    as.data.frame()
test <- test %>% mutate(Patient = stringr::str_split_fixed(rownames(test), "_[A|T|G|C].*", n = 2)[, 1])
id1 <- c("CJME", "CMDI", "HDYU", "HXZH", "LCLI", "WYJU", "WZLA", "ZXME", "ZJLI_0116")
id2 <- c("CZYI", "FHXHBS1", "FYYI", "HEJI", "LAWE", "ZEYI", "ZFXI", "LIPE")
id3 <- c("CJME_0707", "CMDI_0624", "HDYU_0720", "HXZH_0220", "LCLI_0623", "WYJU_0122", "ZXME_0223", "ZJLI_0312")
test$group <- ifelse(test$Patient %in% id1, "R_Pre",
    ifelse(test$Patient %in% id2, "NR_Pre",
        ifelse(test$Patient %in% id3, "R_Post", "NR_Post")
    )
)

test$group1 <- str_split_fixed(test$group, "_", n = 2)[, 1]
test$group2 <- str_split_fixed(test$group, "_", n = 2)[, 2]

## 添加celltype
load("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/Liu_HER2_celltype.RData")
# data_merge <- rbind(HER2_Score, R_Pre_HER2, R_Post_HER2, NR_Pre_HER2, NR_Post_HER2)
test$celltype <- Liu_HER2_celltype[match(rownames(test), rownames(Liu_HER2_celltype)), "celltype"]
colnames(test)[1] <- "values"
test$new_celltype <- ifelse(test$celltype %in% c(0), "cancer cells", "Prolifer cells")

mytheme <- theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.line = element_line(linewidth = 1),
    panel.grid = element_blank(),
    panel.border =  element_blank()
    
)
test$group1 <- factor(test$group1, levels = c('R', 'NR'))
test$group2 <- factor(test$group2, levels = c('Pre', 'Post'))

png("split_violin.png", height = 2000, width = 4000, res = 300)
p <- ggplot(test, aes(x = group1, y = values, fill = group2)) +
    geom_split_violin(trim = T, color = "white") +
    # geom_boxplot(width = 0.1, position = position_dodge(0.6), outlier.shape = NA, height = 0.4) +
    # geom_boxplot(width = 0.1, position = position_dodge2(0.9), fill = "white")+
    geom_point(stat = "summary", fun = mean, position = position_dodge(0.4), size = 4) +
    scale_fill_aaas() +
    stat_summary(
        fun.min = function(x) {
            quantile(x)[2]
        },
        fun.max = function(x) {
            quantile(x)[4]
        },
        geom = "errorbar", color = "black",
        width = 0.1, size = 0.5,
        position = position_dodge(width = 0.4)
    ) +
    theme_classic() +
    mytheme +
    ggpubr::stat_compare_means(method = "wilcox.test") +
    facet_grid(. ~ new_celltype)
p
dev.off()
save(test,file="/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/IFN.RData")



## 绘图
library(ggplot2)
library(patchwork)
library(ggsci)
library(ggpubr)
library(gghalves)
library(ggunchained)

plotSplitViolin <- function(object,x = , y =  fill = split.by =  ){
    ggplot(data = object, mapping = aes())

}



# !------------------------------------------------------------------------------------------------------------------> Angiogenesis <--------------------------------------------------------------------------------------------------------------
signatures <- ("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/Angiogenesis.gmt")
splitBy <- c("R_Pre", "R_Post", "NR_Pre", "NR_Post")
object <- scRNA_CCA
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
test <- Reduce(rbind, flist) %>%
    as.data.frame()
test <- test %>% mutate(Patient = stringr::str_split_fixed(rownames(test), "_[A|T|G|C].*", n = 2)[, 1])
id1 <- c("CJME", "CMDI", "HDYU", "HXZH", "LCLI", "WYJU", "WZLA", "ZXME", "ZJLI_0116")
id2 <- c("CZYI", "FHXHBS1", "FYYI", "HEJI", "LAWE", "ZEYI", "ZFXI", "LIPE")
id3 <- c("CJME_0707", "CMDI_0624", "HDYU_0720", "HXZH_0220", "LCLI_0623", "WYJU_0122", "ZXME_0223", "ZJLI_0312")
test$group <- ifelse(test$Patient %in% id1, "R_Pre",
    ifelse(test$Patient %in% id2, "NR_Pre",
        ifelse(test$Patient %in% id3, "R_Post", "NR_Post")
    )
)

test$group1 <- str_split_fixed(test$group, "_", n = 2)[, 1]
test$group2 <- str_split_fixed(test$group, "_", n = 2)[, 2]

## 添加celltype
load("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/Liu_HER2_celltype.RData")
# data_merge <- rbind(HER2_Score, R_Pre_HER2, R_Post_HER2, NR_Pre_HER2, NR_Post_HER2)
test$celltype <- Liu_HER2_celltype[match(rownames(test), rownames(Liu_HER2_celltype)), "celltype"]
colnames(test)[1] <- "values"
test$new_celltype <- ifelse(test$celltype %in% c(0), "cancer cells", "Prolifer cells")
save(test,file="/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/Angiogenesis.RData")

# !-----------------------------------------------------------------------------------------------------------------> fat acid <------------------------------------------------------------------------------------------------------
signatures <- c("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/fatAcid.gmt")
splitBy <- c("R_Pre", "R_Post", "NR_Pre", "NR_Post")
object <- scRNA_CCA
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

test <- Reduce(rbind, flist) %>%
    as.data.frame()
test <- test %>% mutate(Patient = stringr::str_split_fixed(rownames(test), "_[A|T|G|C].*", n = 2)[, 1])
id1 <- c("CJME", "CMDI", "HDYU", "HXZH", "LCLI", "WYJU", "WZLA", "ZXME", "ZJLI_0116")
id2 <- c("CZYI", "FHXHBS1", "FYYI", "HEJI", "LAWE", "ZEYI", "ZFXI", "LIPE")
id3 <- c("CJME_0707", "CMDI_0624", "HDYU_0720", "HXZH_0220", "LCLI_0623", "WYJU_0122", "ZXME_0223", "ZJLI_0312")
test$group <- ifelse(test$Patient %in% id1, "R_Pre",
    ifelse(test$Patient %in% id2, "NR_Pre",
        ifelse(test$Patient %in% id3, "R_Post", "NR_Post")
    )
)

test$group1 <- str_split_fixed(test$group, "_", n = 2)[, 1]
test$group2 <- str_split_fixed(test$group, "_", n = 2)[, 2]

## 添加celltype
load("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/Liu_HER2_celltype.RData")
# data_merge <- rbind(HER2_Score, R_Pre_HER2, R_Post_HER2, NR_Pre_HER2, NR_Post_HER2)
test$celltype <- Liu_HER2_celltype[match(rownames(test), rownames(Liu_HER2_celltype)), "celltype"]
colnames(test)[1] <- "values"
test$new_celltype <- ifelse(test$celltype %in% c(0), "cancer cells", "Prolifer cells")
save(test,file="/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/fatAcid.RData")

# !---------------------------------------------------------------------------------------------------------------> 细胞凋亡 <-------------------------------------------------------------------------------------------------------
signatures <- c("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/细胞凋亡.gmt")
splitBy <- c("R_Pre", "R_Post", "NR_Pre", "NR_Post")
object <- scRNA_CCA
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

test <- Reduce(rbind, flist) %>%
    as.data.frame()
test <- test %>% mutate(Patient = stringr::str_split_fixed(rownames(test), "_[A|T|G|C].*", n = 2)[, 1])
id1 <- c("CJME", "CMDI", "HDYU", "HXZH", "LCLI", "WYJU", "WZLA", "ZXME", "ZJLI_0116")
id2 <- c("CZYI", "FHXHBS1", "FYYI", "HEJI", "LAWE", "ZEYI", "ZFXI", "LIPE")
id3 <- c("CJME_0707", "CMDI_0624", "HDYU_0720", "HXZH_0220", "LCLI_0623", "WYJU_0122", "ZXME_0223", "ZJLI_0312")
test$group <- ifelse(test$Patient %in% id1, "R_Pre",
    ifelse(test$Patient %in% id2, "NR_Pre",
        ifelse(test$Patient %in% id3, "R_Post", "NR_Post")
    )
)

test$group1 <- str_split_fixed(test$group, "_", n = 2)[, 1]
test$group2 <- str_split_fixed(test$group, "_", n = 2)[, 2]

## 添加celltype
load("/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/Liu_HER2_celltype.RData")
# data_merge <- rbind(HER2_Score, R_Pre_HER2, R_Post_HER2, NR_Pre_HER2, NR_Post_HER2)
test$celltype <- Liu_HER2_celltype[match(rownames(test), rownames(Liu_HER2_celltype)), "celltype"]
colnames(test)[1] <- "values"
test$new_celltype <- ifelse(test$celltype %in% c(0), "cancer cells", "Prolifer cells")
save(test,file="/root/wangje/Project/刘老师/合并3HER2_Epi和Liu_Epi/VISION_result/细胞凋亡.RData")

## !---------------------------------------------------------------------------------------------------------> 合并4组计算的signature scores <------------------------------------------------------------------------------------
