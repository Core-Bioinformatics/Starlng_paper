library(Seurat)
library(dplyr)
library(monocle3)
library(Starlng)

if (basename(getwd()) != "robustness_analysis") {
    setwd("robustness_analysis")
}

n_runs <- 10
sub_percentage <- c(0.9, 0.75, 0.5, 0.25)
n_cores <- 8

autocorr_test <- read.csv("baseline_autocorr_test.csv", row.names = 1)
filtered_genes <- autocorr_test %>% dplyr::filter(morans_I > 0.1, q_value < 0.05) %>% rownames

output_mon_objs <- "subsampling_monocle_objects.qs2"
output_autocorr_tests <- "subsampling_autocorr_tests.qs2"
output_clusters <- "subsampling_clustering_assessments.qs2"

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
    for (sub_perc in sub_percentage) {
        print(paste0("Subsampling ", sub_perc * 100, "% of cells, run ", run, "/", n_runs))

        if (as.character(sub_perc) %in% names(sub_mon_objs[[paste0("run_", run)]])) {
            sub_mon_obj <- sub_mon_objs[[paste0("run_", run)]][[as.character(sub_perc)]]
        } else {
            set.seed(run * 42)
            mon_obj <- qs2::qs_read("baseline_monocle_object.qs2", nthreads = 4)
            sampled_cells <- sample(colnames(mon_obj), size = floor(sub_perc * ncol(mon_obj)), replace = FALSE)

            sub_mon_obj <- mon_obj[, sampled_cells]
            set.seed(run * 42)
            sub_mon_obj <- preprocess_cds(
                cds = sub_mon_obj,
                method = "PCA",
                norm_method = "none",
                num_dim = 50
            )

            set.seed(run * 42)
            sub_mon_obj <- reduce_dimension(
                cds = sub_mon_obj,
                reduction_method = "UMAP",
                preprocess_method = "PCA",
                umap.metric = "cosine",
                umap.n_neighbors = 30,
                umap.min_dist = 0.3
            )

            sub_mon_obj <- custom_learn_graph(
                mon_obj = sub_mon_obj,
                use_partition = FALSE,
                nodes_per_log10_cells = 50,
                learn_graph_control = list(
                    eps = 1e-5,
                    maxiter = 100
                )
            )
            sub_mon_objs[[paste0("run_", run)]][[as.character(sub_perc)]] <- sub_mon_obj
            qs2::qs_save(sub_mon_objs, output_mon_objs, nthreads = 4)
        }

        if (as.character(sub_perc) %in% names(sub_autocorr[[paste0("run_", run)]])) {
            sub_autocorr_test <- sub_autocorr[[paste0("run_", run)]][[as.character(sub_perc)]]
        } else {
            sub_autocorr_test <- graph_test(
                cds = sub_mon_obj,
                neighbor_graph = "principal_graph",
                cores = n_cores / 2
            )
            sub_autocorr[[paste0("run_", run)]][[as.character(sub_perc)]] <- sub_autocorr_test
            qs2::qs_save(sub_autocorr, output_autocorr_tests, nthreads = 4)
        }

        if (as.character(sub_perc) %in% names(sub_clusters[[paste0("run_", run)]])) {
            clustering_assessment <- sub_clusters[[paste0("run_", run)]][[as.character(sub_perc)]]
        } else {
            filtered_genes <- autocorr_test %>% dplyr::filter(morans_I > 0.1, q_value < 0.05) %>% rownames
            gene_loadings <- sub_mon_obj@reduce_dim_aux$"PCA"$model$svd_v
            missing_genes <- setdiff(filtered_genes, rownames(gene_loadings))
            if (length(missing_genes) > 0) {
                gene_loadings <- rbind(gene_loadings, matrix(0, nrow = length(missing_genes), ncol = ncol(gene_loadings)))
                rownames(gene_loadings)[(nrow(gene_loadings) - length(missing_genes) + 1):nrow(gene_loadings)] <- missing_genes
            }
            gene_loadings <- gene_loadings[filtered_genes, ]
            rm(sub_mon_obj)
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

            sub_clusters[[paste0("run_", run)]][[as.character(sub_perc)]] <- clustering_assessment
            qs2::qs_save(sub_clusters, output_clusters, nthreads = 4)
        }
    }
}
