library(Seurat)
library(SCENIC)
library(SingleCellExperiment)
library(SCopeLoomR)

add_cell_annotation <- function(loom, cellAnnotation) {
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}


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
            output_file_loom <- paste0("sergio_", generation_type, count_type, "_", n_tf, "_tfs_seurat.loom")

            so <- qs2::qs_read(input_file, nthreads = 30)
            if (file.exists(output_file_loom)) {
                next
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
            expr_matrix <- GetAssayData(so, assay = "RNA", layer = "data")
            meta_data <- so@meta.data

            loom <- build_loom(file.path(output_file_loom), dgem = expr_matrix)
            loom <- add_cell_annotation(loom, meta_data)
            close_loom(loom)
        }
    }
}

# organs <- c("Liver", "Pancreas", "Lung")
organs <- c("Liver_subset")

for (organ in organs) {
    print(paste0("Processing ", organ))
    input_file <- paste0("cao_", organ, "_filtered_normalized.qs2")
    output_file_loom <- paste0("cao_", organ, ".loom")

    so <- qs2::qs_read(input_file, nthreads = 30)
    gc()
    if (file.exists(output_file_loom)) {
        next
    }
    expr_matrix <- GetAssayData(so, assay = "RNA", layer = "data")
    meta_data <- so@meta.data

    loom <- build_loom(file.path(output_file_loom), dgem = expr_matrix)
    loom <- add_cell_annotation(loom, meta_data)
    close_loom(loom)
}


dts_path <- "../data/immuneCellsSCTransformed.rds"
dts_prefix <- "masld_immune"

so <- readRDS(dts_path)
ca <- readRDS(file.path("..", "starlng_run", paste0(dts_prefix, "_clustassess_output.rds")))
so <- so[rownames(so@assays$SCT), ]
DefaultAssay(so) <- "RNA"
so@assays$SCT <- NULL
so$stable_24_clusters <- factor(ClustAssess::get_clusters_from_clustassess_object(
    clustassess_object = ca,
    feature_type = "HV",
    feature_size = 2000,
    clustering_method = "SLM",
    nclusters = 24
)$"24"$partitions[[1]]$mb, levels = 1:24)
so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so, features = rownames(so))
set.seed(42)
so <- RunPCA(so, features = VariableFeatures(so))
so@reductions$pca@cell.embeddings <- ca$"HV"$"2000"$pca
so <- RunUMAP(so, reduction = "pca", dims = 1:30)
so@reductions$umap@cell.embeddings <- ca$"HV"$"2000"$umap
colnames(so@reductions$umap@cell.embeddings) <- c("umap_1", "umap_2")
so@graphs$SCT_nn <- NULL
so@graphs$SCT_snn <- NULL
so <- FindNeighbors(so, reduction = "pca", dims = 1:30)
so <- FindClusters(so, resolution = 0.5, random.seed = 42)

output_file_loom <- paste0(dts_prefix, ".loom")
expr_matrix <- GetAssayData(so, assay = "RNA", layer = "data")
meta_data <- so@meta.data

loom <- build_loom(file.path(output_file_loom), dgem = expr_matrix)
loom <- add_cell_annotation(loom, meta_data)
close_loom(loom)


# used in the original run of Scenic for the immune cells
# so <- readRDS("immuneCellsSCTransformed.rds")
# so <- SCTransform(so, variable.features.n = 3000, return.only.var.genes = FALSE, verbose = FALSE)
# so <- RunPCA(so, npcs = 50, verbose = FALSE)
# so <- RunUMAP(so, reduction = "pca", dims = 1:50, verbose = FALSE)

# cc_genes <- cc.genes.updated.2019
# so <- CellCycleScoring(so, s.features = cc_genes$s.genes, g2m.features = cc_genes$g2m.genes)
# qs2::qs_save(so, "immuneCellsSCTransformed.qs2")

# expr_matrix <- so@assays$SCT@data
# meta_data <- so@meta.data

# loom <- build_loom(file.path("seurat_info.loom"), dgem = expr_matrix)
# loom <- add_cell_annotation(loom, meta_data)
# close_loom(loom)
