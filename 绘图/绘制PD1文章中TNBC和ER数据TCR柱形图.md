## 对TCR数据进行整理
### 数据整理1
```R
library(dplyr)
library(stringr)
library(qs)

setwd('/root/wangje/Artical_Data/PD1_ER_HER2_TNBC')
# 读入数据
df <- read.csv('./1881-BIOKEY_clonotypes_combined_cohort1.csv')
```
[![piQ2s8e.png](https://z1.ax1x.com/2023/11/05/piQ2s8e.png)](https://imgse.com/i/piQ2s8e)
```R
df <- df %>% mutate(
  TRAaa = str_split_fixed(cdr3s_aa,';', n = 2)[,1] %>% str_replace_all('TRA:',''),
  TRBaa = str_split_fixed(cdr3s_aa, ';', n = 2)[,2] %>% str_replace_all('TRB:',''),
  TRAnt = str_split_fixed(cdr3s_nt,';', n = 2)[,1] %>% str_replace_all('TRA:',''),
  TRBnt = str_split_fixed(cdr3s_nt , ';', n = 2)[,2] %>% str_replace_all('TRB:','')
)
head(df)
```
[![piQREa6.png](https://z1.ax1x.com/2023/11/05/piQREa6.png)](https://imgse.com/i/piQREa6)
### 数据整理2
```R
df1 <- read.csv('./1879-BIOKEY_barcodes_vdj_combined_cohort1.csv')
df1 <- df1 %>% mutate(
  TRAaa = str_split_fixed(cdr3_aa,';', n = 2)[,1] %>% str_replace_all('TRA:',''),
  TRBaa = str_split_fixed(cdr3_aa, ';', n = 2)[,2] %>% str_replace_all('TRB:',''),
  TRAnt = str_split_fixed(cdr3_nt,';', n = 2)[,1] %>% str_replace_all('TRA:',''),
  TRBnt = str_split_fixed(cdr3_nt , ';', n = 2)[,2] %>% str_replace_all('TRB:','')
)
head(df1)
df1$barcode <- paste0(df1$barcode,'-1')
# 写出数据
write.csv(df1, file = './1879-BIOKEY_barcodes_vdj_combined_cohort1_整理.csv',quote = F, row.names = F)
```
[![piQRMMd.png](https://z1.ax1x.com/2023/11/05/piQRMMd.png)](https://imgse.com/i/piQRMMd)
### 提取匹配ER数据的TCR
```R
ER <- qread('NKT_ER_seurat.qs')
ER@meta.data$barcode = rownames(ER@meta.data)
# 合并数据
ER@meta.data <- left_join(ER@meta.data, df1, by = 'barcode')
rownames(ER@meta.data) <- ER$barcode
write.table(ER@meta.data, file = './ER_TCR_metaData.csv',quote = F, row.names = F, sep= '\t')
```
### 提取TNBC数据的TCR
```R
TNBC <- qread('./NKT_TNBC_seurat.qs')
TNBC@meta.data$barcode = rownames(TNBC@meta.data)
# 合并数据
TNBC@meta.data <- left_join(TNBC@meta.data, df1, by = 'barcode')
rownames(TNBC@meta.data) <- TNBC$barcode
write.table(TNBC@meta.data, file = './TNBC_TCR_metaData.csv',quote = F, row.names = F, sep= '\t')
```
### ER TCR数据绘图
#### CD4 T cells
```R
# 筛选出CD4 T细胞
cell1s <- c("CD4+ Tn",'Tfh',"Treg","CD4+ Tem","CD4+ Tcxcl13-ifng","Cycling T")
ER_tcr <- ER_tcr %>% dplyr::select(celltype, TRBaa) %>%
    dplyr::filter(celltype %in% cell1s) %>% na.omit() %>%
    dplyr::filter(TRBaa != '')

# 计算重叠TCR数据
# 存放重叠TCR数目的数据框
tcrOverlap_CD4 <- data.frame(
    matrix(
        nrow = length(unique(ER_tcr$celltype)),
        ncol = length(unique(ER_tcr$celltype))
    )
)
colnames(tcrOverlap_CD4) <- unique(ER_tcr$celltype)
rownames(tcrOverlap_CD4) <- unique(ER_tcr$celltype)

# 使用交集计算不同celltype之间的TCR重叠数
for (i in unique(ER_tcr$celltype)) {
    for (j in unique(ER_tcr$celltype)) {
        tmp1 <- ER_tcr %>% filter(celltype == i) %>% pull(TRBaa)
        tmp2 <- ER_tcr %>% filter(celltype == j) %>% pull(TRBaa)
        tcrOverlap_CD4[i,j] <- length(intersect(tmp1,tmp2))
    }
}
# 将同一个celltype之间的TCR重叠数赋值为0
diag(tcrOverlap_CD4) <- 0

# 将宽数据转换为长数据
tcrOverlap_CD4 <- tcrOverlap_CD4 %>% tibble::rownames_to_column(var = 'cell1') %>%
    pivot_longer(cols = 2:last_col(), names_to = 'cell2', values_to = 'Freq1') %>%
    dplyr::filter(! cell1 == cell2) %>% na.omit()

# 计算每种细胞类型的unique TCR数据
test <- ER_tcr%>%
    group_by(celltype) %>%
    summarise(
        Freq2 = length(unique(TRBaa))
    )
colnames(test)[1] <- 'cell1'
# 计算重叠TCR数据在每种细胞类型的unique TCR数中的比例
tcrOverlap_CD4<- left_join(tcrOverlap_CD4,test, by='cell1')
tcrOverlap_CD4$prop <- tcrOverlap_CD4$Freq1/tcrOverlap_CD4$Freq2
head(tcrOverlap_CD4)

# 绘制柱形图
tcrOverlap_CD4$cell1 <- factor(
    tcrOverlap_CD4$cell1,
    levels = rev(c("CD4+ Tn","CD4+ Tem", "CD4+ Tcxcl13-ifng","Treg","Tfh","Cycling T"))
)

tcrOverlap_CD4$cell2 <- factor(
    tcrOverlap_CD4$cell2,
    levels = c("CD4+ Tn","CD4+ Tem", "CD4+ Tcxcl13-ifng","Treg","Tfh","Cycling T")
)

p1 <- ggplot(tcrOverlap_CD4,aes(x=prop,y=cell1,fill = cell2))+
    geom_bar(color="black", position = position_dodge2(),stat = "identity")+
    labs(y = '', fill = '',x= 'TCR Proportion')+
    coord_cartesian(xlim = c(0,0.5))+
    scale_x_continuous(expand = c(0,0),breaks = seq(0.1,0.4,0.1))+
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        # axis.title.x = element_blank(),
        axis.line = element_line(size = 0.8, colour = "black"),
        axis.text.y = element_text(size = 30 , color = "black"),
        axis.title.y = element_text(size = 35, colour = "black"),
        axis.title.x = element_text(size = 30, colour = "black"),
        axis.ticks.length.y =  unit(0.3,'cm'),
        legend.text = element_text(size = 25, colour = 'black'),
        axis.text.x = element_text(size = 30 ,color="black", angle = 90,hjust = 1, vjust = 0.5)
    )+ scale_fill_manual(values = c('#99c7e5','#1e71ae','#b2df8a','#ff7f01','#fb9a99','#6a3d9a','#fbfb98','#b15928'))

# 保存图片
pdf("C:/Users/Administrator/Desktop/IMG/TCR/New_TCR/ER_CD4_T细胞celltype重叠比例柱形图.pdf", height = 6,width = 14,family = 'ArialMT')
p1
dev.off()
```
#### CD8 T cells
```R


```











