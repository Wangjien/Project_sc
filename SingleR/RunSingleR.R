# ------------------------------------ singleR 函数 ---------------------------------------------
RunSingleR <- function(
    scRNA_data = scRNA_seurat,
    cluster = "integrated_snn_res.0.4",
    refData = "/root/wangje/singleR.RData",
    ...) {
    load(refData)
    pacman::p_load("SingleR", "BiocParallel")
    scTCR_list.singleR <- GetAssayData(scRNA_data)
    clusters <- scRNA_data@meta.data[cluster] # 定义分辨率
    scTCR_list.singleRAll <- SingleR( # nolint
        test = scTCR_list.singleR,
        ref = list(encode, hema, hpca, immune, monaImmune),
        labels = list(
            encode$label.main,
            hema$label.main,
            hpca$label.main,
            immune$label.main,
            monaImmune$label.main
        ),
        clusters = clusters,
        BPPARAM = MulticoreParam(10),
        ...
    )
    celltype <- data.frame(ClusterID = rownames(scTCR_list.singleRAll), celltype = scTCR_list.singleRAll$labels, stringsAsFactors = F)
    return(celltype)
}
RunSingleR(scRNA_data=scRNA_seurat,cluster="integrated_snn_res.0.4")
# 单独运行
pbmc.hesc <-SingleR(test pbmc_for_SingleR,ref list(encode,hema,hpca,immune,monaImmune),
        labels list(
        encode$label.main,
        hema$label.main,
        hpca$label.main,
        immune$label.main,
        monaImmune$label.main),clusters=clusters)
