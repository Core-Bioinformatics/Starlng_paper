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

            input_file <- paste0("sergio_", generation_type, count_type, "_counts_soft_", n_soft_clusters, "_clusters_", n_tf, "_tfs.csv")
            output_file <- paste0("sergio_", generation_type, count_type, "_", n_tf, "_tfs_seurat.qs2")


            expr_matrix <- read.table(input_file, sep = ",", header = FALSE)
            gene_names <- readRDS(paste0("sergio_ground_truth_", n_tf, "_tfs.rds"))$index
            colnames(expr_matrix) <- paste0("cell", seq_len(ncol(expr_matrix)))
            rownames(expr_matrix) <- names(sort(gene_names))

            so <- CreateSeuratObject(counts = expr_matrix)
            ncells_per_clusters <- ncol(expr_matrix) / n_soft_clusters
            so$gt_clusters <- rep(c(paste0("cluster", seq_len(n_soft_clusters))), each = ncells_per_clusters)

            so <- NormalizeData(so)
            so <- FindVariableFeatures(so, nfeatures = 2000)
            so <- ScaleData(so)
            so <- RunPCA(so, npcs = 30)
            so <- RunUMAP(so, dims = 1:30)

            print(DimPlot(so, group.by = "gt_clusters") + ggtitle(paste0("Sergio - ", n_tf, " TFs ", generation_type, " ", count_type, " counts")))

            # qs2::qs_save(so, output_file, nthreads = 30)

        }
    }

}


