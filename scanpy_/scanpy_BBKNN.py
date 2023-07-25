import scanpy as sc
import pandas as pd
import numpy as np 
import diopy
import scrublet as scr
import scipy
import matplotlib.pyplot as plt

data1 = diopy.input.read_h5("/root/wangje/Project/刘老师/大群/20230723/01.h5")
sc.pl.highest_expr_genes(data1, n_top=20, )
sc.pl.violin(data1, ['nCount_RNA', 'nFeature_RNA', 'percent.mt'],
             jitter=0, multi_panel=True)
sc.pp.normalize_total(data1, target_sum=1e4)
sc.pp.log1p(data1)
sc.pp.highly_variable_genes(data1, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes = 5000)
data1.raw = data1
data1 = data1[:, data1.var.highly_variable]
sc.pp.regress_out(data1, ['nCount_RNA', 'percent.mt'])
sc.pp.scale(data1, max_value=10)
sc.tl.pca(data1, svd_solver='arpack', n_comps = 50)
sc.pp.neighbors(data1, n_neighbors=50, n_pcs=50, metric='cosine')
sc.tl.umap(data1)
sc.pl.umap(data1, color=['Patient'])
sc.pp.neighbors(data1, n_neighbors=50, n_pcs=50, metric='cosine')
sc.external.pp.bbknn(data1, batch_key='Patient')
sc.tl.umap(data1)
sc.pl.umap(data1, color=['Patient'])
