# 计算并绘制Endothelials 和 Fibroblasts 的箱线图

### 查看DimPlot
```R
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(stringr)

load('/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/CCA/Data/seurat_cca.RData')
# 绘制分辨率为0.2的时候的
Idents(scRNA_CCA_res) <- 'integrated_snn_res.0.2'
DefaultAssay(scRNA_CCA_res) <- 'RNA'
plotFeature <- function(scRNA_CCA_res = scRNA_CCA_res) {
    features <- c("PDGFRA", "COL1A1", "ACTA2", "PDGFRB", "MCAM", "PECAM1", "CD34", "VWF","TOP2A","MKI67","STMN1","TUBA1B")
    p <- DimPlot(scRNA_CCA_res, label = T, repel = T) +
        theme(legend.position = "bottom") +
        FeaturePlot(scRNA_CCA_res, features = features, ncol = 6) +
        plot_layout(
            design = "AABBBBBB
                      AABBBBBB"
        )
    return(p)
}
png('/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/CCA/Fig/Seurat_CCA_O.2DimPlot结果.png', height = 2000, width = 8000, res= 300 )
plotFeature(scRNA_seurat = scRNA_CCA_res)
dev.off()

Idents(scRNA_CCA_res) <- 'integrated_snn_res.0.15'
DefaultAssay(scRNA_CCA_res) <- 'RNA'
plotFeature <- function(scRNA_CCA_res = scRNA_CCA_res) {
    features <- c("PDGFRA", "COL1A1", "ACTA2", "PDGFRB", "MCAM", "PECAM1", "CD34", "VWF","TOP2A","MKI67","STMN1","TUBA1B")
    p <- DimPlot(scRNA_CCA_res, label = T, repel = T) +
        theme(legend.position = "bottom") +
        FeaturePlot(scRNA_CCA_res, features = features, ncol = 6) +
        plot_layout(
            design = "AABBBBBB
                      AABBBBBB"
        )
    return(p)
}
png('/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/CCA/Fig/Seurat_CCA_O.15DimPlot结果.png', height = 2000, width = 8000, res= 300 )
plotFeature(scRNA_CCA_res = scRNA_CCA_res)
dev.off()


```
[![ppklTa9.md.png](https://s1.ax1x.com/2023/03/03/ppklTa9.md.png)](https://imgse.com/i/ppklTa9)

### 绘制分辨率为0.2时候的箱线图
#### 以所分析亚群的每个样本的全部细胞作为分母
```R
library(ggpubr)
library(patchwork)
library(ggplot2)

# 计算比例
## 添加Traet_assess
id1 <- c("CJME","CMDI","HDYU","HXZH","LCLI","WYJU","WZLA","ZXME","ZJLI_0116")
id2 <- c("CZYI","FHXHBS1","FYYI","HEJI","LAWE","ZEYI","ZFXI","LIPE")
id3 <- c("CJME_0707","CMDI_0624","HDYU_0720","HXZH_0220","LCLI_0623","WYJU_0122","ZXME_0223","ZJLI_0312")
scRNA_CCA_res[['Patient']] <- unlist(strsplit(rownames(scRNA_CCA_res@meta.data), split = '_[A|T|C|G]*$'))
scRNA_CCA_res@meta.data$Treat_assess <-ifelse(scRNA_CCA_res@meta.data$Patient %in% id1, "R_Pre", 
                                    ifelse(scRNA_CCA_res@meta.data$Patient %in% id2,"NR_Pre",
                                    ifelse(scRNA_CCA_res@meta.data$Patient %in% id3,"R_Post","NR_Post")))

ids4 <- c('CJME','CJME_0707','CMDI','CMDI_0624','HDYU','HDYU_0720','HXZH','HXZH_0220','LCLI','LCLI_0623','WYJU','WYJU_0122','ZXME','ZXME_0223','ZJLI_0116','ZJLI_0312',
'CZYI','CZYI_0702','FHXHBS1','FHXHBS2','HEJI','HEJX','LAWE','LAWE_0309','ZEYI','ZEYI_0204','WZLA','FYYI','ZFXI','LIPE')
ids5 <- c('RCJMEprePRLymph','RCJMEpostPRLymph','RCMDIprePRLymph','RCMDIpostPRLymph','RHDYUprePRBreast','RHDYUpostPRBrest','RHXZHprePRBreast','RHXZHpostPRBreast','RLCLIpreCRLymph','RLCLIpostChestWall','RWYJUprePRLymph','RWYJUpostPRLymph','RZXMEprePRBreast','RZXMEpostPRBreast','RZJLIprePRBreast','RZJLIpostPRBreast',
'NRCZYIprePDLiver','NRCZYIpostPDLiver','NRFHXHBSpreSDBreast','NRFHXHBSpostSDBreast','NRHEJIpreSDBreast','NRHEJXpostSDBreast','NRLAWEprePDLiver','NRLAWEpostPDLiver','NRZEYIprePDBreast','NRZEYIpostPDBreast',
'RRWZLApreChestWall','NRFYYIpreSDLiver','NRZFXIprePDBreast','NRLIPEprePDLiver')
ids6 <- c('Lymph','Lymph','Lymph','Lymph','Breast','Breast','Breast','Breast','Lymph','ChestWall','Lymph','Lymph','Breast','Breast','Breast','Breast',
'Liver','Liver','Breast','Breast','Breast','Breast','Liver','Liver','Breast','Breast','ChestWall','Liver','Breast','Liver')
ids7 <- c('R','R','R','R','R','R','R','R','R','R','R','R','R','R','R','R','NR','NR','NR','NR','NR','NR','NR','NR','NR','NR','R','NR','NR','NR')

for(i in 1:length(ids4)){
    print(ids4[i])
    print(ids5[i])
    print(ids6[i])
    print(ids7[i])
    scRNA_CCA_res@meta.data[which(scRNA_CCA_res$Patient == ids4[i]),"Rename"]  <- ids5[i]
    scRNA_CCA_res@meta.data[which(scRNA_CCA_res$Patient == ids4[i]),"Tissue"]  <- ids6[i]
    scRNA_CCA_res@meta.data[which(scRNA_CCA_res$Patient == ids4[i]),"Rseponse"]  <- ids7[i]
}

ids4 <- c(
    "CJME", "CJME_0707", "CMDI", "CMDI_0624", "HDYU", "HDYU_0720", "HXZH", "HXZH_0220", "LCLI", "LCLI_0623", "WYJU", "WYJU_0122",
    "ZXME", "ZXME_0223", "ZJLI_0116", "ZJLI_0312", "CZYI", "CZYI_0702", "FHXHBS1", "FHXHBS2", "HEJI", "HEJX", "LAWE", "LAWE_0309",
    "ZEYI", "ZEYI_0204", "WZLA", "FYYI", "ZFXI", "LIPE"
)
ids5 <- c(
    "P1-Pre", "P1-Post", "P2-Pre", "P2-Post", "P3-Pre", "P3-Post", "P4-Pre", "P4-Post", "P5-Pre", "P5-Post", "P6-Pre", "P6-Post",
    "P7-Pre", "P7-Post", "P8-Pre", "P8-Post", "P10-Pre", "P10-Post", "P11-Pre", "P11-Post", "P12-Pre", "P12-Post", "P13-Pre", "P13-Post",
    "P14-Pre", "P14-Post", "P9-Pre", "P15-Pre", "P16-Pre", "P17-Pre"
)

for (i in 1:length(ids4)) {
    print(ids4[i])
    print(ids5[i])
    scRNA_CCA_res@meta.data[which(scRNA_CCA_res$Patient == ids4[i]), "new_Rename"] <- ids5[i]
}
save(scRNA_CCA_res, file = '/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/CCA/Data/seurat_cca.RData')

## 计算比例
scRNA_CCA_res$new_Patient <- paste0(scRNA_CCA_res$Patient, '--', scRNA_CCA_res$Treat_assess)
Fraction_df  = prop.table(table(scRNA_CCA_res$new_Patient, scRNA_CCA_res$integrated_snn_res.0.15), margin = 1 ) %>% as.data.frame()
Fraction_df <- Fraction_df %>% 
    mutate(
        group = stringr::str_split_fixed(.$Var1, '--', n = 2)[,2],
        Patient = stringr::str_split_fixed(.$Var1, '--', n = 2)[,1]ss
        # sample <- stringr::str_split_fixed(.$Patient, "_", n = 1),
        # sample <- gsub('FHXHBS1|FHXHBS2','FHXHBS',.$sample)
    )
Fraction_df$sample <- stringr::str_split_fixed(Fraction_df$Patient, "_", n = 2)[,1]
Fraction_df$sample <- gsub('FHXHBS1|FHXHBS2','FHXHBS',Fraction_df$sample)
Fraction_df$group <- factor(Fraction_df$group, levels = c(' R_Pre', 'R_Post', 'NR_Pre', 'NR_Post'))
compaired <- list(c("R_Pre", "R_Post"), c("NR_Pre", "NR_Post"), c("R_Pre", "NR_Pre"), c("R_Post", "NR_Post"))

PlotFractionBoxplot <- function(df.fraction=df.fraction,compaired = compaired,ncols = ncols,text=NULL,...){
    # pacman::p_load("ggplot2","ggpubr","patchwork","cowplot","tidyverse","ggrepel")
    plist <- list()
    for (celltype in unique(df.fraction$celltype)) {
        tmp <- df.fraction[df.fraction$celltype == celltype, ]
        p <- ggplot(tmp, aes(x = group, y = Fraction)) +
            geom_boxplot(aes(fill = group), show.legend = F, width = 0.6, outlier.colour = NA) +
            geom_jitter(width = 0.1,size = 2, fill = "black",color="black") + # 绘制散点
            geom_line(aes(group = sample), color = "darkgrey", lwd = 0.5) + # 配对样本间连线
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
            geom_text_repel(aes(group, Fraction, label=sample),size=2)+
            # ggpubr::stat_compare_means(method = "wilcox.test", hide.ns = F, comparisons = compaired, label = "p.format", symnum.args = symnum.args, size = 6) +
            ggpubr::stat_compare_means(method="wilcox.test",hide.ns = F,comparisons = compaired,label="p.formot")+
            scale_fill_manual(values = c("#3b9aac","#a74947","#4470a1","#8da24f","#6d5a87"))
            plist[[celltype]] <- p
        }
	plist <<- plist
	return(plist)
}
result <- PlotFractionBoxplot(df.fraction = Fraction_df,compaired = compaired, ncols = 6)
png('/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/CCA/Fig/CCA_0.15_boxplot.png', height = 2000, width = 10500, res= 300)
wrap_plots(result, ncol =7)
dev.off()
```
[![ppk2HnU.md.png](https://s1.ax1x.com/2023/03/03/ppk2HnU.md.png)](https://imgse.com/i/ppk2HnU)

#### 以每个样本的大群细胞数目作为分母计算比例
```R
# 读入大群的细胞数目
df <- read.table('/root/wangje/Project/刘老师/大群/刘大群每个样本的细胞数目.txt', header=T, sep ='\t')
df_fib_endo <- table(scRNA_CCA_res$new_Patient, scRNA_CCA_res$integrated_snn_res.0.15) %>%
    as.data.frame() %>%
    mutate(
        Patient = str_split_fixed(.$Var1, "--", n = 2)[,1],
        group = str_split_fixed(.$Var1, "--", n = 2)[,2]
    )
df_fib_endo$sample <- stringr::str_split_fixed(df_fib_endo$Patient, "_", n = 2)[,1]
df_fib_endo$sample <- gsub('FHXHBS1|FHXHBS2','FHXHBS',df_fib_endo$sample)
colnames(df_fib_endo)[c(2,3)] <- c('celltype', 'Freq1')
colnames(df)[c(1:2)] <- c('Patient','Freq2')
merge_df <- left_join(df_fib_endo, df, by = 'Patient') %>% 
    mutate(Fraction = round(Freq1/Freq2, digits = 6))
# 绘图
result <- PlotFractionBoxplot(df.fraction = merge_df,compaired = compaired, ncols = 7)
png('/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/CCA/Fig/CCA_0.15_以大群样本细胞数目作为分母_boxplot.png', height = 2000, width = 10500, res= 300)
wrap_plots(result, ncol =7)
dev.off()
````
[![ppkXcRK.md.png](https://s1.ax1x.com/2023/03/03/ppkXcRK.md.png)](https://imgse.com/i/ppkXcRK)
