import ktplotspy as kpy
import numpy as np
import pandas as pd
import anndata as ad
import diopy
import scanpy as sc


# 生成h5对象
# exp = pd.read_csv('./CJME_counts.txt', sep = '\t')
# exp = exp.set_index('gene')
# adata = ad.AnnData(exp.values.T, obs=exp.index, var=exp.columns)
'''
library(Seurat)
library(dior)
library(data.table)

'''


# 读入文件
adata = sc.read_h5ad('output/03.inferCNV/adata.h5')
adata_raw = adata.raw.to_adata()
means = pd.read_csv('means.csv')
pvals = pd.read_csv('pvalues.csv')
decon = pd.read_csv('deconvoluted.csv')

# 合并多个文件的结果
