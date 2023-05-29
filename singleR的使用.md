## singleR的使用

### 1 下载SingleR和需要的数据库

```R
library(BiocManager)
library(Seurat)
library(dplyr)
library(BiocParallel)
BiocManager::install('SingleR') # 安装软件
```

SingleR有7个数据库，其中5个是人类的数据库，其中有两个是小鼠的数据库

```R
hpca.se=HumanPrimaryCellAtlasData() ##第一次载入会下载数据集，可能会慢一些，后面在用时就不用下载了
Blue.se=BlueprintEncodeData() 
Immune.se=DatabaseImmuneCellExpressionData()
Nover.se=NovershternHematopoieticData()
MonacoIm.se=MonacoImmuneData()
ImmGen.se=ImmGenData() #(鼠)
Mouse.se=MouseRNAseqData() #(鼠)
```

### 2 使用SingleR对seurat对象进行注释

#### 1) 使用SingleR中的一个数据库进行注释

```R
meta=pbmc@meta.data #pbmc的meta文件，包含了seurat的聚类结果
pbmc_for_SingleR <- GetAssayData(pbmc, slot="data")
pbmc.hesc <- SingleR(test = pbmc_for_SingleR, 
                     ref = hpca.se, 
                     labels = hpca.se$label.main) # 使用HumanPrimaryCellAtlasData参考数据集，main大类注释，也可使用fine小类注释，不过小类注释准确性不好确定
table(pbmc.hesc[[i]]$labels,meta$seurat_clusters) ##查看新注释的标签与seurat分类的结果的交叠情况

```

#### 2) 使用SingleR中的多个数据库进行注释

##### 人类

```R
scTCR_list.singleR <- GetAssayData(scRNA_seurat, solt = 'data')
clusters <- scRNA_seurat$integrated_snn_res.0.2
scTCR_list.singleRAll <- SingleR(
    test = scTCR_list.singleR,
    ref = list(encode, hema, hpca, immune, monaImmune),
    labels = list(
    encode$label.main,
    hema$label.main,
    hpca$label.main,
    immune$label.main,
    monaImmune$label.main),
    clusters = clusters,
    BPPARAM = BiocParallel::MulticoreParam(6))
celltype = data.frame(ClusterID = rownames(scTCR_list.singleRAll), celltype=scTCR_list.singleRAll$labels, stringsAsFactors = F)
```

#####  小鼠

```R
scTCR_list.singleR <- GetAssayData(scRNA_seurat, solt = 'data')
clusters <- scRNA_seurat$integrated_snn_res.0.2
scTCR_list.singleRAll <- SingleR(
    test = scTCR_list.singleR,
    ref = list(ImmGen.se,Mouse.se),
    labels = list(
        ImmGen.se$label.main,
        Mouse.se$label.main
    ),clusters = clusters,
    BPPARAM = BiocParallel::MulticoreParam(6))
celltype = data.frame(ClusterID = rownames(scTCR_list.singleRAll), celltype=scTCR_list.singleRAll$labels, stringsAsFactors = F)
```



