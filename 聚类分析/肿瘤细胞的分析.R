# //----------------------------------------------------------------
# // 肿瘤细胞使用FastMNN进行分析
# // 分别使用FastMNN的结果和Seurat CCA的结果
# //
# //----------------------------------------------------------------

# // FastMNN的结果
library(Seurat)
library(ggplot2)
library(dplyr)

big.fastmnn = '/root/wangje/wangjien_new/liu-big/scRNA_FastMNN_new.RData'
load(big.fastmnn)



# //----------------------------------------------------------------
# // 使用rPCA进行分析
# //----------------------------------------------------------------


## RunUMAP函数中的dist参数
