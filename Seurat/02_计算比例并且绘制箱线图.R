## 计算比例
CalculateFraction <-
    function(object = object,
             sampleInfo = sampleInfo,
             celltype = celltype,
             group = group,
             ...) {
        if (!require("pacman")) {
            install.packages("pacman")
        } else{
            pacman::p_load("Seurat", "tidyverse", "dplyr", "stringr")
        }
        # 计算每个人，每个celltype的细胞数目
        tmp1 <- object@meta.data %>% 
            group_by(!!!syms(sampleInfo),!!!syms(celltype),!!!syms(group)) %>%
            count(!!!syms(sampleInfo))
        tmp1 <- tmp1 %>% rename_with(~c("Patient","celltype","group","celltype_cells"))
        # print(tmp1)
        # 计算每个人的总细胞数
        tmp2 <- object@meta.data %>% 
            group_by(!!!syms(sampleInfo),!!!syms(group)) %>% 
            count(!!!syms(sampleInfo))
        tmp2 <- tmp2 %>% rename_with(~c("Patient","group","sample_cells"))
        # print(tmp2)
        # 合并两个文件，并且计算比例
        tmp <- tmp1 %>% left_join(tmp2, by = c("Patient","group")) %>% 
            mutate(Fraction = round(celltype_cells/sample_cells,digits = 4))
        print(tmp)
        return(tmp)
    }
# Test 
CalculateFraction(object = scRNA_seurat,sampleInfo = "Patient", celltype = "integrated_snn_res.0.1",group = "Treat_assess")

# Result
# ----------------------------->>>
# Patient      celltype group celltype_cells sample_cells Fraction
# BIOKEY_13_Pre	0	    Pre	   778		  3369	   0.2309
# BIOKEY_13_Pre	1	    Pre	   927		  3369	   0.2752
# BIOKEY_13_Pre	2	    Pre	   114		  3369	   0.0338
# BIOKEY_13_Pre	3	    Pre	   805		  3369	   0.2389
# BIOKEY_13_Pre	4	    Pre	   137		  3369	   0.0407
# BIOKEY_13_Pre	5 	  Pre	   196		  3369	   0.0582
# BIOKEY_13_Pre	6	    Pre	   139		  3369	   0.0413
# BIOKEY_13_Pre	7	    Pre	   241		  3369	   0.0715
# BIOKEY_13_Pre	8	    Pre	   12		    3369	   0.0036
# BIOKEY_13_Pre	9	    Pre	   20		    3369	   0.0059

# ----------------------------->>>
## 绘制boxplot
PlotFractionBoxplot <- function(df.fraction=df.fraction,compaired = compaired,ncols = ncols,text=NULL,...){
    pacman::p_load("ggplot2","ggpubr","patchwork","cowplot","tidyverse","ggrepel")
    plist <- list()
    for (celltype in unique(df.fraction$celltype)) {
        tmp <- df.fraction[df.fraction$celltype == celltype, ]
        p <- ggplot(tmp, aes(x = group, y = Fraction)) +
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
	return(plist)
}

## ------------------------------------------------------>>>
## 绘制boxplot
## ------------------------------------------------------>>>
# df.merge
# Patient          celltype group celltype_cells sample_cells Fraction
# 1    CJME  B & Plasma cells R_Pre             84        16380   0.0051
# 2    CJME Endothelial cells R_Pre           1098        16380   0.0670
# 3    CJME  Epithelial cells R_Pre           1811        16380   0.1106
# 4    CJME       Fibroblasts R_Pre           3367        16380   0.2056
# 5    CJME     Myeloid cells R_Pre           9491        16380   0.5794
# 6    CJME      NK & T cells R_Pre            529        16380   0.0323
### BoxPlot
PlotFractionBoxplot <- function(df.fraction=df.fraction,compaired = compaired,ncols = ncols,p.format = TRUE,text=NULL,...){
    pacman::p_load("ggplot2","ggpubr","patchwork","cowplot","tidyverse","ggrepel","patchwork")
    plist <- list()
    if(p.format){
    for (celltype in unique(df.fraction$celltype)) {
        tmp <- df.fraction[df.fraction$celltype == celltype, ]
        p <- ggplot(tmp, aes(x = group, y = Fraction)) +
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
            # geom_text_repel(aes(group, prop, label=sample),size=2)+
            # ggpubr::stat_compare_means(method = "wilcox.test", hide.ns = F, comparisons = compaired, label = "p.format", symnum.args = symnum.args, size = 6) +
            ggpubr::stat_compare_means(method="wilcox.test",hide.ns = F,comparisons = compaired,label="p.formot")+
            scale_fill_manual(values = c("#3b9aac","#a74947","#4470a1","#8da24f","#6d5a87"))
            plist[[celltype]] <- p
        }
    }else{
         for (celltype in unique(df.fraction$celltype)) {
        tmp <- df.fraction[df.fraction$celltype == celltype, ]
        p <- ggplot(tmp, aes(x = group, y = Fraction)) +
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
            # geom_text_repel(aes(group, prop, label=sample),size=2)+
            ggpubr::stat_compare_means(method = "wilcox.test", hide.ns = F, comparisons = compaired, label = "p.format", symnum.args = symnum.args, size = 6) +
            # ggpubr::stat_compare_means(method="wilcox.test",hide.ns = F,comparisons = compaired,label="p.formot")+
            scale_fill_manual(values = c("#3b9aac","#a74947","#4470a1","#8da24f","#6d5a87"))
            plist[[celltype]] <- p
      }
    }
    plist <<- plist
    return(patchwork::wrap_plots(plist,ncol = ncols))
}
# 举例
symnum.args <- list(cutpoints = c(0, 0.005, 0.01, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns"))
compaired = list(c('R_Pre','R_Post'),c('NR_Pre','NR_Post'),c('NR_Pre','R_Pre'),c('NR_Post','R_Post'))
png('./NKT_FastMNN_比例箱线图_单样本全部细胞.png', height = 5000, width= 6500, res=300)
PlotFractionBoxplot(df.fraction = df.merge,compaired = compaired,ncols = 4)
dev.off()







id1 <- c("CJME","CMDI","HDYU","HXZH","LCLI","WYJU","WZLA","ZXME","ZJLI_0116")
id2 <- c("CZYI","FHXHBS1","FYYI","HEJI","LAWE","ZEYI","ZFXI","LIPE")
id3 <- c("CJME_0707","CMDI_0624","HDYU_0720","HXZH_0220","LCLI_0623","WYJU_0122","ZXME_0223","ZJLI_0312")
test1$group <-ifelse(test1$sample %in% id1, "R_Pre", 
                                    ifelse(test1$sample %in% id2,"NR_Pre",
                                    ifelse(test1$sample %in% id3,"R_Post","NR_Post")))
