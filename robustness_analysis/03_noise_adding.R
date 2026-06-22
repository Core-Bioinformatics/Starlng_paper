library(Seurat)
library(dplyr)
library(monocle3)
library(Starlng)

if (basename(getwd()) != "robustness_analysis") {
    setwd("robustness_analysis")
}

n_runs <- 10
gaussian_noise_sd <- c(0.1, 0.25, 0.5, 1, 2)
n_cores <- 8

autocorr_test <- read.csv("baseline_autocorr_test.csv", row.names = 1)
filtered_genes <- autocorr_test %>% dplyr::filter(morans_I > 0.1, q_value < 0.05) %>% rownames

output_mon_objs <- "noise_matrix.qs2"
output_autocorr_tests <- "noise_autocorr_tests.qs2"
output_clusters <- "noise_clustering_assessments.qs2"

if (file.exists(output_mon_objs)) {
    sub_mon_objs <- qs2::qs_read(output_mon_objs, nthreads = 4)
} else {
    sub_mon_objs <- list()
}
if (file.exists(output_autocorr_tests)) {
    sub_autocorr <- qs2::qs_read(output_autocorr_tests, nthreads = 4)
} else {
    sub_autocorr <- list()
}
if (file.exists(output_clusters)) {
    sub_clusters <- qs2::qs_read(output_clusters, nthreads = 4)
} else {
    sub_clusters <- list()
}

for (run in seq_len(n_runs)) {
    if (isFALSE(paste0("run_", run) %in% names(sub_mon_objs))) {
        sub_mon_objs[[paste0("run_", run)]] <- list()
    }
    if (isFALSE(paste0("run_", run) %in% names(sub_autocorr))) {
        sub_autocorr[[paste0("run_", run)]] <- list()
    }
    if (isFALSE(paste0("run_", run) %in% names(sub_clusters))) {
        sub_clusters[[paste0("run_", run)]] <- list()
    }
    for (noise_sd in gaussian_noise_sd) {
        print(paste0("Adding Gaussian noise with SD ", noise_sd, ", run ", run, "/", n_runs))
        if (isFALSE(as.character(noise_sd) %in% names(sub_clusters[[paste0("run_", run)]]))) {
            set.seed(run * 42)
            mon_obj <- qs2::qs_read("baseline_monocle_object.qs2", nthreads = 4)
            noise_sd_matrix <- matrix(rnorm(n = nrow(mon_obj) * ncol(mon_obj), mean = 0, sd = noise_sd), nrow = nrow(mon_obj), ncol = ncol(mon_obj))
            mon_obj@assays@data$normalized_data <- mon_obj@assays@data$normalized_data + noise_sd_matrix
            mon_obj@assays@data$counts <- mon_obj@assays@data$counts + noise_sd_matrix
            # convert to sparse
            mon_obj@assays@data$counts <- Matrix::Matrix(mon_obj@assays@data$counts, sparse = TRUE)
            # rm(noise_sd_matrix)
            # gc()

            RhpcBLASctl::blas_set_num_threads(n_cores)
            set.seed(run * 42)
            mon_obj <- preprocess_cds(
                cds = mon_obj,
                method = "PCA",
                norm_method = "none",
                num_dim = 50
            )
            RhpcBLASctl::blas_set_num_threads(1)

            set.seed(run * 42)
            mon_obj <- reduce_dimension(
                cds = mon_obj,
                reduction_method = "UMAP",
                preprocess_method = "PCA",
                umap.metric = "cosine",
                umap.n_neighbors = 30,
                umap.min_dist = 0.3
            )

            mon_obj <- custom_learn_graph(
                mon_obj = mon_obj,
                use_partition = FALSE,
                nodes_per_log10_cells = 50,
                learn_graph_control = list(
                    eps = 1e-5,
                    maxiter = 100
                )
            )
            sub_mon_objs[[paste0("run_", run)]][[as.character(noise_sd)]] <- noise_sd_matrix
            qs2::qs_save(sub_mon_objs, output_mon_objs, nthreads = 4)
        } 

        if (as.character(noise_sd) %in% names(sub_autocorr[[paste0("run_", run)]])) {
            sub_autocorr_test <- sub_autocorr[[paste0("run_", run)]][[as.character(noise_sd)]]
        } else {
            sub_autocorr_test <- graph_test(
                cds = mon_obj,
                neighbor_graph = "principal_graph",
                cores = n_cores / 2
            )
            sub_autocorr[[paste0("run_", run)]][[as.character(noise_sd)]] <- sub_autocorr_test
            qs2::qs_save(sub_autocorr, output_autocorr_tests, nthreads = 4)
        }

        if (as.character(noise_sd) %in% names(sub_clusters[[paste0("run_", run)]])) {
            clustering_assessment <- sub_clusters[[paste0("run_", run)]][[as.character(noise_sd)]]
        } else {
            filtered_genes <- autocorr_test %>% dplyr::filter(morans_I > 0.1, q_value < 0.05) %>% rownames
            gene_loadings <- mon_obj@reduce_dim_aux$"PCA"$model$svd_v[filtered_genes, ]
            rm(mon_obj)
            gc()
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
                number_repetitions = 100
            )
            if (n_cores > 1) {
                parallel::stopCluster(cl)
                foreach::registerDoSEQ()
            }

            sub_clusters[[paste0("run_", run)]][[as.character(noise_sd)]] <- clustering_assessment
            qs2::qs_save(sub_clusters, output_clusters, nthreads = 4)
        }
    }
}

