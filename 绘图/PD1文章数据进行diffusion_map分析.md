## 读入数据,提取需要的数据，转化为H5对象
```R
library(Seurat)
library(qs)
library(dior)

# 读入数据
ER <- qread('/root/wangje/Artical_Data/PD1_ER_HER2_TNBC/NKT_ER_seurat.qs')
TNBC <- qread('/root/wangje/Artical_Data/PD1_ER_HER2_TNBC/NKT_TNBC_seurat.qs')

# 提取需要的clltype
ER_sub <- ER[, ER$celltype %in% c('CD4+ Tn','CD4+ Tcxcl13-ifng','Tfh')]
TNBC_sub <- TNBC[, TNBC$celltype %in% c('CD4+ Tn','CD4+ Tcxcl13-ifng','Tfh')]

# 转为H5数据
dior::write_h5(ER_sub, file = 'ER_CD4_tn_Th1-Tfh.h5')
dior::write_h5(TNBC_sub, file = 'TNBC_CD4_tn_Th1-Tfh.h5')
```
## 读入H5数据，进行diffusion分析
```python
import os
import scanpy as sc
import diopy
import matplotlib.pyplot as plt
import warnings
import numpy as np
import scvelo as scv

# ER
ER = diopy.input.read_h5('./ER_CD4_tn_Th1-Tfh.h5')
sc.pp.neighbors(ER,n_pcs=40)
sc.tl.diffmap(ER, n_comps=40)
# plot
sc.pl.diffmap(ER,color="celltype2",projection="3d")




# TNBC
TNBC = diopy.input.read_h5('./TNBC_CD4_tn_Th1-Tfh.h5')
sc.pp.neighbors(TNBC,n_pcs=40)
sc.tl.diffmap(TNBC, n_comps=40)

```
