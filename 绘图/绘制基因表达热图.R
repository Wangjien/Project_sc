library(Seurat)
library(Seurat)
library(patchwork)
library(tidyverse)
library(dplyr)
library(ggrepel)


# 绘制火山图
plot_volcanes = function(DEGs, title= 'R_Pre vs NR_Pre') {
    cat('DGEs的维度为:',dim(DEGs))
    # DEGs$log10padgs = ifelse(DEGs$p_val_adj==0, 1, log10(DEGs$p_val_adj))
    DEGs$log10padgs = -log10(DEGs$p_val_adj)
    DEGs$gene = rownames(DEGs)
    DEGs$group = ifelse(DEGs$avg_log2FC > 0.5 & DEGs$p_val_adj < 0.05, 'Up', 
                    ifelse(DEGs$avg_log2FC < -0.5 & DEGs$p_val_adj < 0.05, 'Down', 'Normal')) 
    # 绘制散点图
    p = ggplot(DEGs, aes(x = avg_log2FC, y = log10padgs))+ 
        geom_point(aes(color = group)) + 
        scale_color_manual(values = c('red','black', 'red'))+ 
        geom_text_repel(data = DEGs[which(DEGs$group %in% c('Up','Down')),], aes(label = gene))+ 
        labs(x = expression("log"[2]('fold change')), 
             y = expression("log"[10]("P Value adjusted")),
             title = title )+ 
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = 'None')

    return(p)
}

markers <- FindMarkers(scRNA_seurat,
                       ident.1 = "R_Pre",
                       ident.2 = 'NR_Pre',
                       group.by = 'Treat_assess',
                       subset.ident = c('R_Pre','NR_Pre'))

DEGs = function(data,ident1 = 'R_Pre', ident2 = 'NR_Pre', subset_ident = c('R_Pre','NR_Pre')){
    result = FindMarkers(
        data,
        ident.1 = ident1,
        ident.2 = ident2,
        group.by = 'Treat_assess',
        subset.ident = subset_ident
    )
    return(result)
}                       

# 添加SingleR注释
for(i in seq_len(nrow(celltype))){
    print(celltype[i,1:2])
    scRNA_seurat@meta.data[which(scRNA_seurat$integrated_snn_res.0.3 == celltype[i,1]),'singler'] = celltype[i,2]
}

# 绘制出图片
png('./res0.3.png', height = 2000, width = 6000, res=300)
DimPlot(scRNA_seurat, group.by='integrated_snn_res.0.3', label =T, raster=F) + theme(legend.position='bottom') -> p1
DimPlot(scRNA_seurat, group.by='celltype', label=T, repel=T, raster=F) + theme(legend.position='bottom') -> p2
DimPlot(scRNA_seurat, group.by='singler', label=T, repel=T, raster=F) + theme(legend.position='bottom') -> p3
p1|p2|p3
dev.off()
# 查看FeaturePlot
png('./FeaturePlot.png', height=4000,width=6000,res=300)
FeaturePlot(scRNA_seurat, features=c('IRF7','CXCR3','LILRA4','CD79A','MZB1','MS4A1'), ncol=3, raster=F)
dev.off()

# 查看不同分辨率中的表达分布
png('./ALL_featurePlot.png', height= 6000,width= 6000, res=300)
plotFeature(scRNA_data=scRNA_seurat,choose="Feature",col_num=6,marker.list=marker.list)
dev.off()


# 绘制热图
gene_exp = DotPlot(scRNA_seurat, group.by='new_group', features=c('CD274','PDCD1LG2','CTLA4'))
gene_data = gene_exp$data %>% as.data.frame() %>% mutate(
    group = str_split_fixed(id,'_', n=2)[,2],
    celltype = str_split_fixed(id, '_', n=2)[,1]
) %>% select(avg.exp.scaled,group,celltype,features.plot ) 
pivot_wider(id_cols = group, names_from = celltype, values_from = avg.exp.scaled)

## 绘制热图
test1 = gene_data %>% filter(group == 'R_Pre') %>% pivot_wider(names_from=celltype, values_from=avg.exp.scaled) %>% select(-group) %>% as.data.frame()
test2 = gene_data %>% filter(group == 'R_Post') %>% pivot_wider(names_from=celltype, values_from=avg.exp.scaled) %>% select(-group) %>% as.data.frame()
test3 = gene_data %>% filter(group == 'NR_Pre') %>% pivot_wider(names_from=celltype, values_from=avg.exp.scaled) %>% select(-group) %>% as.data.frame()
test4 = gene_data %>% filter(group == 'NR_Post') %>% pivot_wider(names_from=celltype, values_from=avg.exp.scaled) %>% select(-group) %>% as.data.frame()

rownames(test1) = test1$features.plot
rownames(test2) = test2$features.plot
rownames(test3) = test3$features.plot
rownames(test4) = test4$features.plot

test1 = test1 %>% select(-1)
test2 = test2 %>% select(-1)
test3 = test3 %>% select(-1)
test4 = test4 %>% select(-1)

## 绘制热图
h1 = ComplexHeatmap::pheatmap(test1,
        cluster_rows = FALSE,
        cluster_col = FALSE
    
    )

h2 = ComplexHeatmap::pheatmap(test2,
        cluster_rows = FALSE,
        cluster_col = FALSE
    )

h3 = ComplexHeatmap::pheatmap(test3,
        cluster_rows = FALSE,
        cluster_col = FALSE
    )

h4 = ComplexHeatmap::pheatmap(test4,
        cluster_rows = FALSE,
        cluster_col = FALSE
    )

ht_list = h1 %v% h2 %v% h3 %v% h4
png('/root/wangje/Project/刘老师/大群/Fig/大群基因热图.png',height = 2000,width = 2000,res=300)
draw(ht_list, column_km = 1)
dev.off()

## 使用solt=‘counts’
scRNA_seurat$new_group = paste0(scRNA_seurat$celltype, scRNA_seurat$Treat_assess)
scRNA_seurat$Treat_assess = factor(scRNA_seurat$Treat_assess, levels = c('NR_Post','NR_Pre', 'R_Post', 'R_Pre'))
test1 = AverageExpression(scRNA_seurat, features = c('CD274','PDCD1LG2','CTLA4'), group.by = 'new_group', slot = 'counts')
test2 = AverageExpression(scRNA_seurat, features = c('CD274','PDCD1LG2','CTLA4'), group.by = 'new_group', slot = 'data') 
test3 = AverageExpression(scRNA_seurat, features = c('CD274','PDCD1LG2','CTLA4'), group.by = 'new_group', slot = 'scale.data')

h1 = ComplexHeatmap::pheatmap(test1,
        cluster_rows = FALSE,
        cluster_col = FALSE
    )
png('/root/wangje/Project/刘老师/大群/Fig/大群基因热图_geneCount.png',height = 1000,width = 2000,res=300)
h1
dev.off()


h2 = ComplexHeatmap::pheatmap(test2,
        cluster_rows = FALSE,
        cluster_col = FALSE
    )
png('/root/wangje/Project/刘老师/大群/Fig/大群基因热图_geneData.png',height = 1000,width = 2000,res=300)
h2
dev.off()

##1111111111111111111111111111111 绘制基因比较小提琴图
plot_geneVlonlin = function(data, genes = NULL, group = NULL){
    if(is.null(genes) && is.null(group)){
        stop('请输入变量.')
    }
    if(all(genes %in% rownames(scRNA_seurat))){
        genes = genes
    }else{
        genes = genes[genes %in% rownames(scRNA_seurat)]
        cat(genes[!(genes %in% rownames(scRNA_seurat))],'不在seurat对象中'。)
    }
    gene_tmp = Seurat::FetchData(data, vars=c(genes, group))
    # 将宽数据转换为长数据
    long_df <<- tidyr::pivot_longer(
        gene_tmp,
        cols = 1:(ncol(gene_tmp)-1),
        names_to = 'gene',
        values_to = 'values'
    )
}

long_df$Treat_assess = factor(long_df$Treat_assess, levels = c('R_Pre','R_Post','NR_Pre','NR_Post'))
gene_tmp = Seurat::FetchData(DC, vars=c(rev(genes), 'Treat_assess'))
p = ggplot(data = long_df, aes(x = Treat_assess, y = values, fill = Treat_assess))+
        geom_violin() +
        ggpubr::stat_compare_means(method="wilcox.test",
            hide.ns = F,
            comparisons = list(c('R_Pre','R_Post'),c('NR_Pre','NR_Post'),c('R_Pre','NR_Pre'),c('R_Post','NR_Post')),
            label="p.formot") + 
            theme_bw()+
            theme(axis.title.x = element_blank(),
            axis.text.x = element_text(size = 15,colour = 'black', angle = 45, hjust = 1, vjust = 1),
            axis.title.y = element_text(size = 15, color = 'black'),
            strip.text.x = element_text(size = 25,colour = 'black'))+
            facet_wrap(.~gene,ncol = 5)
ggsave(filename = '/root/wangje/Project/刘老师/大群/Fig/DC_genes_小提琴图.png',height = 16,width = 20, plot=p, bg = 'white')

test = DotPlot(DC, features= rev(genes), group.by = 'new_group')
test = as.data.frame(test$data)

test1 = spread(test[3:5],key = id, value = avg.exp.scaled)
rownames(test1) = test1$features.plot
test1 = test1 %>% select(-1)
colnames(test1) =c('DC-NR_Post','DC-NR_Pre','DC-R_POst','DC-R_Pre')
test1 = test1 %>% select(4,3,2,1)
h1 = ComplexHeatmap::Heatmap(test1,
        name = 'avg.exp.score',
        cluster_rows = F,
        show_column_dend = FALSE,
        row_split = rep(c('IFN-responsiveness',
                        'Immune cell attraction',
                        'Co-stimulation',
                        'Antigen presentation & processing',
                        'Differentiation'),
                        times=c(11,2,1,5,1)),
        row_title_rot = 0,
        column_order = paste0('DC','-',c('R_Pre','R_Post','NR_Pre','NR_Post')),
        column_split = rep(LETTERS[1:4],1),
        column_title = NULL,
        row_names_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 20),
        column_names_gp = gpar(fontsize = 20)
    )
png('/root/wangje/Project/刘老师/大群/Fig/DC_genes_热图.png', height = 3000,width = 3000,res=300)
h1
dev.off()

# Myeloids
library(Seurat)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggplot2)
library(stringr)
library(ComplexHeatmap)
library(tidyverse)

setwd('/root/wangje/Project/刘老师/Myeloids/Data/')
load('/root/wangje/Project/刘老师/Myeloids/Data/Myeloids_CCA.RData')

# 查看分群
unique(scRNA_seurat$celltype)
myeloids = scRNA_seurat[,scRNA_seurat$celltype %in% c( "Macro_APOE",
                                                        "Macro_CCL2",
                                                        "Mono_CD14",
                                                        "Macro_SLC2A",
                                                        "Macro_LYVE1",
                                                        "Macro_MT1G",
                                                        "Macro_CCL5")]
myeloids$new_group = paste0(myeloids$celltype,'-', myeloids$Treat_assess)                                                
test = DotPlot(myeloids, features = c('CD274','PDCD1LG2','PDCD1'),group.by = 'new_group')
test= test$data
test$group = str_split_fixed(test$id, '-',n=2)[,2]
test$celltype = str_split_fixed(test$id, '-', n= 2)[,1]
test = test %>% select(3,4,5)

test_R_Pre = test %>% filter(grepl('-R_Pre$',id))
wide_R_Pre <- pivot_wider(test_R_Pre, names_from = id, values_from = avg.exp.scaled) %>% as.data.frame()
rownames(wide_R_Pre) = wide_R_Pre$features.plot
wide_R_Pre = wide_R_Pre%>% select(-1)
colnames(wide_R_Pre) = str_split_fixed(colnames(wide_R_Pre),'-', n= 2)[,1]
# 绘图
h1 = ComplexHeatmap::Heatmap(wide_R_Pre,
        name = 'avg.exp.score',
        cluster_rows = F,
        show_column_dend = FALSE,
        row_split = rep('R_Pre',
                        times=3),
        row_title_rot = 0,
        column_title = NULL,
        row_names_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 20),
        column_names_gp = gpar(fontsize = 20)
    )

tidy_data = function(pattern = '-NR_Pre'){
    tmp = test %>% filter(grepl(pattern = pattern ,id))
    wide_data <- pivot_wider(tmp, names_from = id, values_from = avg.exp.scaled) %>% as.data.frame()
    rownames(wide_data) = wide_data$features.plot
    wide_data = wide_data %>% select(-1)
    colnames(wide_data) = str_split_fixed(colnames(wide_data),'-', n= 2)[,1]
    return(wide_data)
}

plot_data = function(data, row_name='R_Pre'){
    h1 = ComplexHeatmap::Heatmap(data,
            name = 'avg.exp.score',
            cluster_rows = F,
            show_column_dend = FALSE,
            row_split = rep(row_name,
                            times=3),
            row_title_rot = 0,
            column_title = NULL,
            row_names_gp = gpar(fontsize = 20),
            row_title_gp = gpar(fontsize = 20),
            column_names_gp = gpar(fontsize = 20)
        )
    return(h1)
}

h1 = plot_data(data=wide_R_Pre, row_name = 'R_Pre')
NR_Pre = tidy_data(pattern = '-NR_Pre')
h2 = plot_data(data=NR_Pre, row_name = 'NR_Pre')
R_Post = tidy_data(pattern = '-R_Post')
h3 = plot_data(data=R_Post, row_name = 'R_Post')
NR_Post = tidy_data(pattern = '-NR_Post')
h4 = plot_data(data=NR_Post, row_name = 'NR_Post')

ht_list = h1 %v% h3 %v% h2 %v% h4
png('./基因热图.png',height = 2000,width = 2500,res=300)
draw(ht_list, column_km = 1)
dev.off()
