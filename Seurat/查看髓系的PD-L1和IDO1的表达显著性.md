# 查看髓系的PD-L1和IDO1的表达显著性

### 髓系ID01和PD-L1表达量比较

```R
.libPaths('/root/wangje/miniconda3/lib/R/library')
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(ggsci)
.libPaths(c('/root/wangje/miniconda3/lib/R/library','/root/wangje/R/x86_64-conda-linux-gnu-library/4.2'))
# 读入数据
load('/root/wangje/Project/刘老师/Myeloids/Data/Myeloids_CCA.RData')
#########################
#### 计算比例 PD-l1
scRNA_seurat$Treat_assess <- gsub("CRPR-|SDPD-", "", scRNA_seurat$Treat_assess)
scRNA_seurat$Treat_assess <- factor(scRNA_seurat$Treat_assess, levels = c("R_Pre", "R_Post", "NR_Pre", "NR_Post"))
scRNA_seurat$new_celltype <- paste0(scRNA_seurat$celltype, "__", scRNA_seurat$Treat_assess)

data1 <- VlnPlot(scRNA_seurat, group.by = "new_celltype", features = c("CD274"))$data
data1$group <- str_split_fixed(data1$ident, "__", n = 2)[, 2]
data1$celltype <- str_split_fixed(data1$ident, "__", n = 2)[, 1]
#########################
##### 计算IDO1
data1 <- VlnPlot(scRNA_seurat, group.by = "new_celltype", features = c("IDO1"))$data
data1$group <- str_split_fixed(data1$ident, "__", n = 2)[, 2]
data1$celltype <- str_split_fixed(data1$ident, "__", n = 2)[, 1]

# 绘图函数
PlotViolin <- function(data = data , x= x , y = y, fill = fill, ncol = ncol, y_label = y_label, point = NULL){
    plist <- list()
    if(is.null(point)){
        for (i in unique(data$celltype)) {
            tmp <- data %>% filter(celltype == i)
            p <- ggplot(tmp, aes(x = !!sym(x), y = !!sym(y), fill=!!sym(fill)))+
                geom_violin(width=0.5)+
                theme_classic()+
                labs(y = y_label, x = '', title = i)+
                theme(plot.background = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size = 20,family = 'Arial'),
                      plot.title = element_text(size = 25,hjust = 0.5,family = 'Arial'),
                      axis.text.y = element_text(size = 25,color='black',family = 'Arial'),
                      axis.text.x = element_text(size = 25,color = 'black',family = 'Arial',hjust = 1, vjust = 1, angle = 45))+
                scale_fill_manual(values = c("#3b9aac","#a74947","#4470a1","#8da24f","#6d5a87"))+
                stat_compare_means(comparisons = list(c('R_Pre','NR_Pre'),c('R_Post','NR_Post'),c('R_Pre','R_Post'),c('NR_Pre','NR_Post')), method = 'wilcox.test', label = "p.format") # 将p值以显著性标签的形式添加到图中   
            plist[[i]] <- p 
            print(i)
        }
        fig <- wrap_plots(plist, ncol = ncol)
        return(fig)
        if(is.null(ncol)){
            return(plist)
        }
    }else{
        for (i in unique(data$celltype)) {
            tmp <- data %>% filter(celltype == i)
            p <- ggplot(tmp, aes(x = !!sym(x), y = !!sym(y), fill=!!sym(fill)))+
                geom_violin(width=0.5)+
                geom_jitter(width=0.3,size=0.1)+
                theme_classic()+
                labs(y = y_label, x = '', title = i)+
                theme(plot.background = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size = 20,family = 'Arial'),
                      plot.title = element_text(size = 25,hjust = 0.5,family = 'Arial'),
                      axis.text.y = element_text(size = 25,color='black',family = 'Arial'),
                      axis.text.x = element_text(size = 25,color = 'black',family = 'Arial',hjust = 1, vjust = 1, angle = 45))+
                scale_fill_manual(values = c("#3b9aac","#a74947","#4470a1","#8da24f","#6d5a87"))+
                stat_compare_means(comparisons = list(c('R_Pre','NR_Pre'),c('R_Post','NR_Post'),c('R_Pre','R_Post'),c('NR_Pre','NR_Post')), method = 'wilcox.test', label = "p.format") # 将p值以显著性标签的形式添加到图中   
            plist[[i]] <- p 
            print(i)
        }   
        fig <- wrap_plots(plist, ncol = ncol)
        return(fig)
        if(is.null(ncol)){
            return(plist)
        }
    }
}

p1 <- PlotViolin(data = data1, x = 'group', y = 'CD274', fill = 'group', ncol = 6, y_label = 'PD-L1 Expression',point = NULL)
p2 <- PlotViolin(data = data2, x = 'group', y = 'IDO1', fill = 'group', ncol = 6, y_label = 'IDO1 Expression',point = NULL)

png('/root/wangje/Project/刘老师/Myeloids/Fig/01_ido1_Vlolin_去除point.png', height = 8000,width = 12000,res=300 )
p1/p2
dev.off()

p1 <- PlotViolin(data = data1, x = 'group', y = 'CD274', fill = 'group', ncol = 6, y_label = 'PD-L1 Expression',point = '')
p2 <- PlotViolin(data = data2, x = 'group', y = 'IDO1', fill = 'group', ncol = 6, y_label = 'IDO1 Expression',point = '')

png('/root/wangje/Project/刘老师/Myeloids/Fig/01_ido1_Vlolin.png', height = 8000,width = 12000,res=300 )
p1/p2
dev.off()
```

[![pp8GhW9.md.png](https://s1.ax1x.com/2023/03/16/pp8GhW9.md.png)](https://imgse.com/i/pp8GhW9)

[![pp8Gosx.md.png](https://s1.ax1x.com/2023/03/16/pp8Gosx.md.png)](https://imgse.com/i/pp8Gosx)

### 髓系新的分群，计算表达量

```R
scRNA_seurat$new_celltype <- ifelse(
    scRNA_seurat$celltype %in% grep("^Macro", unique(scRNA_seurat$celltype), value = T), "Macrophage", ifelse(
        scRNA_seurat$celltype %in% c("migDC", "cDC2_CD1C", "cDC1_CLEC9A"), "DC cells",
        ifelse(
            scRNA_seurat$celltype %in% c("Mono_CD14"), "Mono_CD14", scRNA_seurat$celltype
        )
    )
)
scRNA_seurat$Treat_assess <- gsub("CRPR-|SDPD-", "", scRNA_seurat$Treat_assess)
scRNA_seurat$Treat_assess <- factor(scRNA_seurat$Treat_assess, levels = c("R_Pre", "R_Post", "NR_Pre", "NR_Post"))
scRNA_seurat$new_celltype <- paste0(scRNA_seurat$new_celltype, "__", scRNA_seurat$Treat_assess)

data1 <- VlnPlot(scRNA_seurat, group.by = "new_celltype", features = c("CD274"))$data
data1$group <- str_split_fixed(data1$ident, "__", n = 2)[, 2]
data1$celltype <- str_split_fixed(data1$ident, "__", n = 2)[, 1]

data2 <- VlnPlot(scRNA_seurat, group.by = "new_celltype", features = c('IDO1'))$data
data2$group <- str_split_fixed(data1$ident, "__", n = 2)[, 2]
data2$celltype <- str_split_fixed(data1$ident, "__", n = 2)[, 1]

p3 <- PlotViolin(data = data1,x = 'group', y= 'CD274',fill = 'group',ncol = 5,y_label = 'PD-L1 Expression', point = '')
p4 <- PlotViolin(data = data2,x = 'group', y= 'IDO1',fill = 'group',ncol = 5,y_label = 'IDO1 Expression', point = '')

png('/root/wangje/Project/刘老师/Myeloids/Fig/髓系_大群_PDCD1和ID01表达量.png', height = 4000,width = 8000,res = 300)
p3/p4
dev.off()


p5 <- PlotViolin(data = data1,x = 'group', y= 'CD274',fill = 'group',ncol = 5,y_label = 'PD-L1 Expression', point = NULL)
p6 <- PlotViolin(data = data2,x = 'group', y= 'IDO1',fill = 'group',ncol = 5,y_label = 'IDO1 Expression', point = NULL)

png('/root/wangje/Project/刘老师/Myeloids/Fig/髓系_大群_PDCD1和ID01表达量_去除point.png', height = 4000,width = 8000,res = 300)
p5/p6
dev.off()

```

[![pp8GF2R.md.png](https://s1.ax1x.com/2023/03/16/pp8GF2R.md.png)](https://imgse.com/i/pp8GF2R)

[![pp8GnaD.md.png](https://s1.ax1x.com/2023/03/16/pp8GnaD.md.png)](https://imgse.com/i/pp8GnaD)

