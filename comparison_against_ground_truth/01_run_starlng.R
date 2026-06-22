library(Seurat)
library(dplyr)
library(Starlng)

if (basename(getwd()) != "comparison_against_ground_truth") {
    setwd("comparison_against_ground_truth")
}
n_tfs <- c(5, 10, 15, 20)

for (n_tf in n_tfs) {
    for (count_type in c("raw", "noisy")) {
        for (generation_type in c("", "dynamic_ds12_")) {
            print(paste0("Processing ", n_tf, " TFs", generation_type, " ", count_type, " counts"))
            n_soft_clusters <- ceiling(n_tf * 0.6)


            input_file <- file.path("..", "data", paste0("sergio_", generation_type, count_type, "_", n_tf, "_tfs_seurat.qs2"))

            so <- qs2::qs_read(input_file)
            expr_matrix <- GetAssayData(so, assay = "RNA", layer = "data")
            colnames(expr_matrix) <- as.character(colnames(expr_matrix))
            rownames(expr_matrix) <- as.character(rownames(expr_matrix))
            cell_mask <- apply(expr_matrix, 2, function(x) sum(as.numeric(x)) > 0)
            expr_matrix <- expr_matrix[, cell_mask, drop = FALSE]

            mtd_df <- so@meta.data[cell_mask, , drop = FALSE]
            pca_emb <- Embeddings(so, reduction = "pca")[cell_mask, 1:30]
            umap_emb <- Embeddings(so, reduction = "umap")[cell_mask, ]
            rm(so)
            gc()

            app_folder <- paste0("starlng/sergio_", generation_type, count_type, "_", n_tf, "_tfs_starlng_app")
            if (file.exists(file.path(app_folder, "objects", "module_summaries.h5"))) {
                next
            }

            test_run <- starlng_write_app_default(
                folder_path = app_folder,
                app_title_name = paste0("Sergio - ", n_tf, " TFs ", generation_type, " ", count_type, " counts"),
                expression_matrix = expr_matrix,
                metadata_df = mtd_df,
                pca_embedding = pca_emb,
                umap_embedding = umap_emb,
                gene_filtering_function = function(df) {
                    df %>% dplyr::filter(morans_I > 0.01, q_value < 0.05) %>% rownames
                },
                verbose = TRUE,
                nthreads = 10,
                skip_tf_identification = TRUE,
                learn_graph_parameters = list(
                    nodes_per_log10_cells = 30,
                    learn_graph_controls = list(
                        eps = 1e-5,
                        maxiter = 10,
                        prune_graph = TRUE,
                        minimal_branch_len = 5
                    )
                ),
                clustering_parameters = list(
                    "n_neighbours" = seq(from = 5, to = 50, by = 5),
                    "graph_type" = "snn",
                    "prune_value" = -1,
                    "resolutions" = list(
                        "RBConfigurationVertexPartition" = seq(from = 0.1, to = 2.5, by = 0.01)
                    ),
                    "number_iterations" = 5,
                    "number_repetitions" = 100
                ),
                ecc_threshold = 0.85
            )
        }
    }
}
