library(Seurat)
library(rhdf5)
library(dplyr)
library(Matrix)
library(ggplot2)

if (basename(getwd()) != "data") {
    setwd("data")
}

target_celltypes <- list(
    # "Liver" = c("Hepatoblasts"),
    # "Pancreas" = c("Ductal cells", "Acinar cells", "Islet endocrine cells")
    "Lung" = c("Bronchiolar and alveolar epithelial cells")
)
target_organs <- names(target_celltypes)
obj_names <- paste0("cao_", target_organs, "_gene_count.RDS")

count_matrices <- lapply(obj_names, function(obj_name) {
    obj_full_path <- file.path(obj_name)
    return(readRDS(obj_full_path))
})
names(count_matrices) <- target_organs

metadata_genes <- readRDS("cao_df_gene.RDS")
metadata_genes$gene_short_name <- as.character(metadata_genes$gene_short_name)
metadata_pc_genes <- metadata_genes %>% filter(gene_type == "protein_coding")
metadata_cells <- readRDS("cao_df_cell.RDS")
rownames(metadata_cells) <- metadata_cells$sample


so_objects <- lapply(names(count_matrices), function(name) {
    count_matrix <- count_matrices[[name]]
    common_genes <- intersect(rownames(count_matrix), rownames(metadata_pc_genes))
    keep_genes <- metadata_pc_genes[common_genes, ]
    count_matrix <- count_matrix[common_genes, ]
    print(sum(count_matrix))
    comparison_colsum <- Matrix::colSums(count_matrix)

    i <- match(keep_genes$gene_short_name, unique(keep_genes$gene_short_name))
    j <- match(rownames(keep_genes), rownames(count_matrix))
    M <- sparseMatrix(
        i = i,
        j = j,
        x = 1,
        dims = c(length(unique(keep_genes$gene_short_name)), nrow(count_matrix)),
    )

    aggregate_matrix <- M %*% count_matrix
    rownames(aggregate_matrix) <- unique(keep_genes$gene_short_name)
    print(sum(aggregate_matrix))
    obtained_colsum <- Matrix::colSums(aggregate_matrix)
    print(all.equal(comparison_colsum, obtained_colsum))

    so <- CreateSeuratObject(
        counts = aggregate_matrix
    ) 
    
    so@meta.data <- metadata_cells[colnames(so), ]

    subset_celltypes <- target_celltypes[[name]]
    subset_celltypes <- paste0(name, "-", subset_celltypes)
    print(subset_celltypes)
    print(unique(so$Organ_cell_lineage))
    so <- subset(so, subset = Organ_cell_lineage %in% subset_celltypes)

    return(so)
})
names(so_objects) <- names(count_matrices)

for (organ in target_organs) {
    so <- so_objects[[organ]]
    mt_genes <- sort(grep("^MT-", rownames(so), value = FALSE, ignore.case = TRUE))

    rp_genes <- sort(c(
        grep("^RP[SL]\\d+", rownames(so), value = FALSE, ignore.case = TRUE),
        grep("^MRP[SL]\\d+", rownames(so), value = FALSE, ignore.case = TRUE)
    ))

    print(rownames(so)[mt_genes])
    print(rownames(so)[rp_genes])

    so$nFeature_RNA <- Matrix::colSums(GetAssayData(so, assay = "RNA", layer = "counts") > 0)
    so$nCount_RNA <- Matrix::colSums(GetAssayData(so, assay = "RNA", layer = "counts"))
    so$percent_mt <- PercentageFeatureSet(so, features = rownames(so)[mt_genes])
    so$percent_rp <- PercentageFeatureSet(so, features = rownames(so)[rp_genes])

    print("Total sum of reads incident to mt_genes")
    print(sum(GetAssayData(so, assay = "RNA", layer = "counts")[mt_genes, ]))
    print("Total sum of reads incident to rp_genes")
    print(sum(GetAssayData(so, assay = "RNA", layer = "counts")[rp_genes, ]))
    so <- subset(so, features = -c(mt_genes, rp_genes))
    so$nFeature_RNA <- Matrix::colSums(GetAssayData(so, assay = "RNA", layer = "counts") > 0)
    so$nCount_RNA <- Matrix::colSums(GetAssayData(so, assay = "RNA", layer = "counts"))

    saveRDS(so, paste0("cao_", organ, "_filtered.RDS"))

    so_objects[[organ]] <- so
}

# filtering
for (organ in target_organs) {
    so <- so_objects[[organ]]
    so <- subset(so, nFeature_RNA > 500 & nCount_RNA < 5e3 & percent_rp < 10)
    print(ncol(so))
    Idents(so) <- "Organ"
    saveRDS(so, paste0("cao_", organ, "_filtered.RDS"))
    so_objects[[organ]] <- so
}

# QC plots
mtd_rna <- c("nFeature_RNA", "nCount_RNA", "percent_rp")
for (organ in target_organs) {
    so <- so_objects[[organ]]
    so <- subset(so, nFeature_RNA > 500 & nCount_RNA < 5e3 & percent_rp < 10)
    print(ncol(so))
    Idents(so) <- "Organ"
    plots <- lapply(mtd_rna, function(mtd) {
        VlnPlot(so, features = mtd, ncol = 3, pt.size = 0, log = FALSE) +
        ggtitle(paste0(organ)) +
        theme(legend.position = "none")
    })
    print(cowplot::plot_grid(plotlist = plots, nrow = 1))
}

# normalisation
for (organ in target_organs) {
    so <- so_objects[[organ]]
    so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
    so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
    so <- ScaleData(so, features = rownames(so), verbose = FALSE)
    so <- RunPCA(so)
    set.seed(42)
    so <- RunUMAP(so, dims = 1:30, reduction = "pca")

    print(DimPlot(so, group.by = "RT_group"))
    print(DimPlot(so, group.by = "Organ_cell_lineage"))
    qs2::qs_save(so, paste0("cao_", organ, "_filtered_normalized.qs2"), nthreads = 30)
}
