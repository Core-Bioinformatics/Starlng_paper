library(Seurat)
library(dplyr)
library(monocle3)
library(Starlng)

if (basename(getwd()) != "robustness_analysis") {
    setwd("robustness_analysis")
}

app_path <- file.path("..", "starlng_run", "masld_immune_high_moran_starlng_app")

monocle_obj <- qs2::qs_read(file.path(app_path, "objects", "monocle_object.qs2"))
n_cores <- 8

set.seed(42)
monocle_obj <- preprocess_cds(
    cds = monocle_obj,
    method = "PCA",
    norm_method = "none",
    num_dim = 50
)

set.seed(42)
monocle_obj <- reduce_dimension(
    cds = monocle_obj,
    reduction_method = "UMAP",
    preprocess_method = "PCA",
    umap.metric = "cosine",
    umap.n_neighbors = 30,
    umap.min_dist = 0.3
)

monocle_obj <- custom_learn_graph(
    mon_obj = monocle_obj,
    use_partition = FALSE,
    nodes_per_log10_cells = 50,
    learn_graph_control = list(
        eps = 1e-5,
        maxiter = 100
    )
)

autocorr_test <- graph_test(
    cds = monocle_obj,
    neighbor_graph = "principal_graph",
    cores = n_cores / 2
)

near_zero_var_genes <- apply(exprs(monocle_obj), 1, function(x) var(as.numeric(x)) < 1e-5)
autocorr_test <- autocorr_test[!near_zero_var_genes, ]
monocle_obj <- monocle_obj[!near_zero_var_genes, ]

qs2::qs_save(monocle_obj, "baseline_monocle_object.qs2", nthreads = 4)
write.csv(autocorr_test, "baseline_autocorr_test.csv", row.names = TRUE)

filtered_genes <- autocorr_test %>% dplyr::filter(morans_I > 0.1, q_value < 0.05) %>% rownames

gene_loadings <- monocle_obj@reduce_dim_aux$"PCA"$model$svd_v[filtered_genes, ]
if (n_cores > 1) {
    cl <- parallel::makePSOCKcluster(n_cores)
    doParallel::registerDoParallel(cl)
}

clustering_assessment <- clustering_pipeline(
    embedding = gene_loadings,
    n_neighbours = seq(25, 50, by = 5),
    graph_type = "snn",
    prune_value = -1,
    resolutions = list(
        "RBConfigurationVertexPartition" = seq(from = 0.1, to = 1.5, by = 0.05)
    ),
    number_iterations = 5,
    number_repetitions = 100,
    merge_identical_partitions = TRUE
)

if (n_cores > 1) {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
}

qs2::qs_save(clustering_assessment, "baseline_clustering_assessment.qs2", nthreads = 4)