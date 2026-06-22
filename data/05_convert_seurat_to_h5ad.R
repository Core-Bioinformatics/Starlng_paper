library(SeuratDisk)
library(anndataR)
library(Seurat)
library(dplyr)
library(ggplot2)

if (basename(getwd()) != "data") {
    setwd("data")
}
n_tfs <- c(5, 10, 15, 20)

for (n_tf in n_tfs) {
    for (count_type in c("raw", "noisy")) {
        for (generation_type in c("", "dynamic_ds12_")) {
            print(paste0("Processing ", n_tf, " TFs", generation_type, " ", count_type, " counts"))
            n_soft_clusters <- ceiling(n_tf * 0.6)

            input_file <- paste0("sergio_", generation_type, count_type, "_", n_tf, "_tfs_seurat.qs2")
            output_file_h5seurat <- paste0("sergio_", generation_type, count_type, "_", n_tf, "_tfs_seurat.h5ad")

            so <- qs2::qs_read(input_file, nthreads = 30)
            if (file.exists(output_file_h5seurat)) {
                file.remove(output_file_h5seurat)
            }
            
            expr_mat <- GetAssayData(so, assay = "RNA", slot = "counts")
            mask_cells <- colSums(expr_mat) == 0
            mask_genes <- rowSums(expr_mat) == 0
            if (any(mask_cells)) {
                so <- so[, !mask_cells]
            }
            if (any(mask_genes)) {
                so <- so[!mask_genes, ]
            }
            so@assays$RNA@layers$scale.data <- NULL

            so <- UpdateSeuratObject(so)
            so[["RNA"]] <- as(so[["RNA"]], "Assay")

            MuDataSeurat::WriteH5AD(so, output_file_h5seurat)

        }
    }
}

# organs <- c("Liver", "Pancreas", "Lung")
organs <- c("Liver_subset")

for (organ in organs) {
    print(paste0("Processing ", organ))
    input_file <- paste0("cao_", organ, "_filtered_normalized.qs2")
    output_file_h5seurat <- paste0("cao_", organ, ".h5ad")

    so <- qs2::qs_read(input_file, nthreads = 30)
    gc()
    if (file.exists(output_file_h5seurat)) {
        file.remove(output_file_h5seurat)
    }
    expr_mat <- GetAssayData(so, assay = "RNA", slot = "counts")
    mask_cells <- colSums(expr_mat) == 0
    mask_genes <- rowSums(expr_mat) == 0
    if (any(mask_cells)) {
        so <- so[, !mask_cells]
    }
    if (any(mask_genes)) {
        so <- so[!mask_genes, ]
    }
    so@assays$RNA@layers$scale.data <- NULL
    write_h5ad(so, output_file_h5seurat)
}


dts_path <- "../data/immuneCellsSCTransformed.rds"
dts_prefix <- "masld_immune"

so <- qs2::qs_read("masld_immune_filtered_normalized.qs2", nthreads = 30)

# so <- readRDS(dts_path)
# ca <- readRDS(file.path("..", "starlng_run", paste0(dts_prefix, "_clustassess_output.rds")))
# so <- so[rownames(so@assays$SCT), ]
# DefaultAssay(so) <- "RNA"
# so@assays$SCT <- NULL
# so$stable_24_clusters <- factor(ClustAssess::get_clusters_from_clustassess_object(
#     clustassess_object = ca,
#     feature_type = "HV",
#     feature_size = 2000,
#     clustering_method = "SLM",
#     nclusters = 24
# )$"24"$partitions[[1]]$mb, levels = 1:24)
# so <- NormalizeData(so)
# so <- FindVariableFeatures(so)
# so <- ScaleData(so, features = rownames(so))
# set.seed(42)
# so <- RunPCA(so, features = VariableFeatures(so))
# so@reductions$pca@cell.embeddings <- ca$"HV"$"2000"$pca
# so <- RunUMAP(so, reduction = "pca", dims = 1:30)
# so@reductions$umap@cell.embeddings <- ca$"HV"$"2000"$umap
# colnames(so@reductions$umap@cell.embeddings) <- c("umap_1", "umap_2")
# so@graphs$SCT_nn <- NULL
# so@graphs$SCT_snn <- NULL
# so <- FindNeighbors(so, reduction = "pca", dims = 1:30)
# so <- FindClusters(so, resolution = 0.5, random.seed = 42)
# expr_mat <- GetAssayData(so, assay = "RNA", layers = "counts")
# mask_cells <- colSums(expr_mat) == 0
# mask_genes <- rowSums(expr_mat) == 0
# if (any(mask_cells)) {
#     so <- so[, !mask_cells]
# }
# if (any(mask_genes)) {
#     so <- so[!mask_genes, ]
# }
# delete scale./data
LayerData(so, assay = "RNA", layer = "scale.data") <- NULL
write_h5ad(so, paste0(dts_prefix, ".h5ad"))

