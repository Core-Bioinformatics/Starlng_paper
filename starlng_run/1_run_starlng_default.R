library(ClustAssess)
library(Seurat)
library(Starlng)

library(dplyr)

if (basename(getwd()) != "starlng_run") {
    setwd("starlng_run")
}
# devtools::load_all("/mnt/d/Starlng")
cao_organ <- c("Liver_subset")

for (organ in cao_organ) {
    print(organ)
    so <- qs2::qs_read(paste0("../data/cao_", organ, "_filtered_normalized.qs2"), nthreads = 30)
    expr_matrix <- GetAssayData(so, assay = "RNA", layer = "data")
    genes_near_zero_var <- apply(expr_matrix, 1, function(x) var(as.numeric(x)) < 1e-5)
    expr_matrix <- expr_matrix[!genes_near_zero_var, ]
    colnames(expr_matrix) <- as.character(colnames(expr_matrix))
    rownames(expr_matrix) <- as.character(rownames(expr_matrix))

    mtd_df <- so@meta.data
    pca_emb <- Embeddings(so, reduction = "pca")[, 1:30]
    umap_emb <- Embeddings(so, reduction = "umap")
    rm(so)
    gc()

    test_run <- starlng_write_app_default(
        folder_path = paste0("cao_", tolower(organ), "_starlng_app"),
        app_title_name = paste0("Cao et al. - ", organ),
        expression_matrix = expr_matrix,
        metadata_df = mtd_df,
        pca_embedding = pca_emb,
        umap_embedding = umap_emb,
        gene_filtering_function = function(df) {
            df %>% dplyr::filter(morans_I > 0.1, q_value < 0.05) %>% rownames
        },
        verbose = TRUE,
        nthreads = 5,
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
                "RBConfigurationVertexPartition" = seq(from = 0.1, to = 2.5, by = 0.1)
            ),
            "number_iterations" = 5,
            "number_repetitions" = 100

        )
    )
}

