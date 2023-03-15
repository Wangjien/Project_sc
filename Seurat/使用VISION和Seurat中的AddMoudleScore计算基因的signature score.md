# 使用VISION和Seurat中的AddMoudleScore计算基因的signature score

## 使用VISION

* 计算不同基因集的score

```R
library(Seurat)
library(patchwork)
library(dplyr)
library(VISION)
library(Seurat)

############################################
### 读入Seurat文件
load("/root/wangje/Project/刘老师/Myeloids/Data/Myeloids_CCA.RData")

## 添加其他信息
scRNA_CCA <- scRNA_seurat
id1 <- c("CJME", "CMDI", "HDYU", "HXZH", "LCLI", "WYJU", "WZLA", "ZXME", "ZJLI_0116")
id2 <- c("CZYI", "FHXHBS1", "FYYI", "HEJI", "LAWE", "ZEYI", "ZFXI", "LIPE")
id3 <- c("CJME_0707", "CMDI_0624", "HDYU_0720", "HXZH_0220", "LCLI_0623", "WYJU_0122", "ZXME_0223", "ZJLI_0312")
scRNA_CCA[["Patient"]] <- unlist(strsplit(rownames(scRNA_CCA@meta.data), split = "_[A|T|C|G]*$"))
scRNA_CCA@meta.data$Treat_assess <- ifelse(scRNA_CCA@meta.data$Patient %in% id1, "R_Pre",
  ifelse(scRNA_CCA@meta.data$Patient %in% id2, "NR_Pre",
    ifelse(scRNA_CCA@meta.data$Patient %in% id3, "R_Post", "NR_Post")
  )
)
SplitBy <- c("R_Pre", "R_Post", "NR_Pre", "NR_Post")
## 计算 VISION Score
CaculateViosion <- function(data = data, signatures = signatures) {
  flist <- list()
  for (group in unique(SplitBy)) {
    print(group)
    tmp <- data[, data$Treat_assess == group]
    Score <- Vision(tmp,
      signatures = signatures,
      projection_methods = NULL
    ) %>% analyze()
    flist[[group]] <- Score@SigScores
  }
  return(flist)
}
## 合并数据
plist <- list()
for (i in paste0("res", 1:17)) {
  res_test <- Reduce(rbind, get(i))
  res_test <- as.data.frame(res_test) %>% tibble::rownames_to_column("cells")
  res_test <- res_test %>%
    mutate(
      Rename = anno[match(res_test$cell, rownames(anno)), "Renames"],
      Patient = anno[match(res_test$cell, row.names(anno)), "Patient"],
      celltype = anno[match(res_test$cell, row.names(anno)), "celltype"],
      Treat_assess = anno[match(res_test$cell, row.names(anno)), "Treat_assess"]
    )

plist_new <- purrr::map(plist, function(x) x + theme(plot.title = element_text(size = 10, hjust = 0.5)))
png("/root/wangje/Project/刘老师/Myeloids/Fig/03_Macro_Vicion.png", height = 4000, width = 8500, res = 300)
wrap_plots(plist_new, ncol = 5)
dev.off()



```

