library(dplyr)
library(ggplot2)
library(Starlng)
library(ClustAssess)

if (basename(getwd()) != "robustness_analysis") {
    setwd("robustness_analysis")
}

baseline_clustering <- qs2::qs_read("baseline_clustering_assessment.qs2")
baseline_best_config <- select_best_configuration(baseline_clustering$clusters_list)
ecc_threshold <- 0.9
freq_threshold <- 30
considered_k_app_threshold <- 5
baseline_stable_clusters <- choose_stable_clusters(
    baseline_clustering$clusters_list[[baseline_best_config[[1]]]][[baseline_best_config[[2]]]]$k,
    ecc_threshold = ecc_threshold,
    freq = freq_threshold
)

str(baseline_clustering, max.level = 2)

#### SUBSAMPLE ####
sub_clustering <- qs2::qs_read("subsampling_clustering_assessments.qs2")
sub_stable_clusters <- list()
for (run in names(sub_clustering)) {
    for (sub_perc in names(sub_clustering[[run]])) {
        sub_best_config <- select_best_configuration(sub_clustering[[run]][[sub_perc]]$clusters_list)
        current_stable_clust <- choose_stable_clusters(
            sub_clustering[[run]][[sub_perc]]$clusters_list[[sub_best_config[[1]]]][[sub_best_config[[2]]]]$k,
            ecc_threshold = 0, #ecc_threshold,
            freq = 0 # freq_threshold
        )

        if (isFALSE(sub_perc %in% names(sub_stable_clusters))) {
            sub_stable_clusters[[sub_perc]] <- list()
        }

        sub_stable_clusters[[sub_perc]][[run]] <- current_stable_clust
    }
}

sub_per_run_ecc <- list()
sub_per_run_df <- NULL
for (sub_perc in names(sub_stable_clusters)) {
    sub_per_run_ecc[[sub_perc]] <- list()
    k_list_per_run <- lapply(sub_stable_clusters[[sub_perc]], function(x) names(x))
    print(paste0("Subsampling ", sub_perc, "% of cells"))

    # select the k that appears at least 5 time
    k_freq <- table(unlist(k_list_per_run))
    selected_k <- names(k_freq)[k_freq >= considered_k_app_threshold]
    selected_k <- intersect(selected_k, names(baseline_stable_clusters))

    for (k in selected_k) {
        available_clusters <- lapply(sub_stable_clusters[[sub_perc]], function(x) {
            if (k %in% names(x)) {
                return(x[[k]]$partitions[[1]]$mb)
            } else {
                return(NULL)
            }
        })
        available_clusters <- available_clusters[!sapply(available_clusters, is.null)]
        sub_per_run_ecc[[sub_perc]][[k]] <- element_consistency(
            available_clusters
        )
        sub_per_run_df <- rbind(
            sub_per_run_df,
            data.frame(
                sub_percentage = sub_perc,
                k = k,
                ecc = sub_per_run_ecc[[sub_perc]][[k]]
            )
        )
    }
}
sub_per_run_df$k <- factor(sub_per_run_df$k, levels = stringr::str_sort(unique(sub_per_run_df$k), numeric = TRUE))
pdf("subsampling_consistency_per_k.pdf", width = 5, height = 7.5)
ggplot(sub_per_run_df, aes(x = k, y = ecc)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~sub_percentage, ncol = 2) +
    theme_bw() +
    labs(
        title = "Element consistency of stable clusters across subsampling runs",
        x = "Subsampling percentage",
        y = "Element consistency"
    )
dev.off()

pdf("subsampling_overall_consistency.pdf", width = 6, height = 6)
ggplot(sub_per_run_df, aes(x = sub_percentage, y = ecc)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    labs(
        title = "Element consistency of stable clusters across subsampling runs",
        x = "Subsampling percentage",
        y = "Element consistency"
    )
dev.off()

sub_against_gt_ecc <- list()
sub_gt_df <- NULL
for (sub_perc in names(sub_stable_clusters)) {
    sub_per_run_ecc[[sub_perc]] <- list()
    k_list_per_run <- lapply(sub_stable_clusters[[sub_perc]], function(x) names(x))
    print(paste0("Subsampling ", sub_perc, "% of cells"))
    # print(k_list_per_run)


    # select the k that appears at least 5 time
    k_freq <- table(unlist(k_list_per_run))
    selected_k <- names(k_freq)[k_freq >= considered_k_app_threshold]
    selected_k <- intersect(selected_k, names(baseline_stable_clusters))
    print(length(selected_k))

    for (k in selected_k) {
        available_clusters <- lapply(sub_stable_clusters[[sub_perc]], function(x) {
            if (k %in% names(x)) {
                return(x[[k]]$partitions[[1]]$mb)
            } else {
                return(NULL)
            }
        })
        available_clusters <- available_clusters[!sapply(available_clusters, is.null)]
        sub_against_gt_ecc[[sub_perc]][[k]] <- element_agreement(
            reference_clustering = baseline_stable_clusters[[k]]$partitions[[1]]$mb,
            clustering_list = available_clusters
        )

        sub_gt_df <- rbind(
            sub_gt_df,
            data.frame(
                sub_percentage = sub_perc,
                k = k,
                ecc = sub_against_gt_ecc[[sub_perc]][[k]]
            )
        )
    }
}
sub_gt_df$k <- factor(sub_gt_df$k, levels = stringr::str_sort(unique(sub_gt_df$k), numeric = TRUE))
print(sub_gt_df %>% group_by(sub_percentage) %>% summarise(avg = mean(ecc), med = median(ecc)))

pdf("subsampling_comparison_to_ground_truth_per_k.pdf", width = 5, height = 7.5)
ggplot(sub_gt_df, aes(x = k, y = ecc)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~sub_percentage, ncol = 2) +
    theme_bw() +
    labs(
        title = "Element agreement of stable clusters against ground truth across subsampling runs",
        x = "Subsampling percentage",
        y = "Element agreement"
    )
dev.off()

pdf("subsampling_overall_comparison_to_ground_truth.pdf", width = 6, height = 6)
ggplot(sub_gt_df, aes(x = sub_percentage, y = ecc)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    labs(
        title = "Element agreement of stable clusters against ground truth across subsampling runs",
        x = "Subsampling percentage",
        y = "Element agreement"
    )
dev.off()

#### NOISE ####
noise_clustering <- qs2::qs_read("noise_clustering_assessments.qs2")
noise_stable_clusters <- list()
for (run in names(noise_clustering)) {
    for (noise_sd in names(noise_clustering[[run]])) {
        noise_best_config <- select_best_configuration(noise_clustering[[run]][[noise_sd]]$clusters_list)
        current_stable_clust <- choose_stable_clusters(
            noise_clustering[[run]][[noise_sd]]$clusters_list[[noise_best_config[[1]]]][[noise_best_config[[2]]]]$k,
            ecc_threshold = 0, # ecc_threshold,
            freq = 0 #freq_threshold
        )

        if (isFALSE(noise_sd %in% names(noise_stable_clusters))) {
            noise_stable_clusters[[noise_sd]] <- list()
        }

        noise_stable_clusters[[noise_sd]][[run]] <- current_stable_clust
    }
}

noise_per_run_ecc <- list()
noise_per_run_df <- NULL
for (noise_sd in names(noise_stable_clusters)) {
    noise_per_run_ecc[[noise_sd]] <- list()
    k_list_per_run <- lapply(noise_stable_clusters[[noise_sd]], function(x) names(x))
    print(paste0("Noise SD: ", noise_sd))

    k_freq <- table(unlist(k_list_per_run))
    selected_k <- names(k_freq)[k_freq >= min(considered_k_app_threshold, max(k_freq) %/% 2)]
    selected_k <- intersect(selected_k, names(baseline_stable_clusters))

    for (k in selected_k) {
        available_clusters <- lapply(noise_stable_clusters[[noise_sd]], function(x) {
            if (k %in% names(x)) {
                return(x[[k]]$partitions[[1]]$mb)
            } else {
                return(NULL)
            }
        })
        available_clusters <- available_clusters[!sapply(available_clusters, is.null)]
        noise_per_run_ecc[[noise_sd]][[k]] <- element_consistency(
            available_clusters
        )
        noise_per_run_df <- rbind(
            noise_per_run_df,
            data.frame(
                noise_sd = noise_sd,
                k = k,
                ecc = noise_per_run_ecc[[noise_sd]][[k]]
            )
        )
    }
}
noise_per_run_df$k <- factor(noise_per_run_df$k, levels = stringr::str_sort(unique(noise_per_run_df$k), numeric = TRUE))
pdf("noise_consistency_per_k.pdf", width = 5, height = 7.5)
ggplot(noise_per_run_df, aes(x = k, y = ecc)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~noise_sd, ncol = 2) +
    theme_bw() +
    labs(
        title = "Element consistency of stable clusters across noise runs",
        x = "Noise SD",
        y = "Element consistency"
    )
dev.off()

pdf("noise_overall_consistency.pdf", width = 6, height = 6)
ggplot(noise_per_run_df, aes(x = noise_sd, y = ecc)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    labs(
        title = "Element consistency of stable clusters across noise runs",
        x = "Noise SD",
        y = "Element consistency"
    )
dev.off()

noise_against_gt_ecc <- list()
noise_gt_df <- NULL
for (noise_sd in names(noise_stable_clusters)) {
    noise_against_gt_ecc[[noise_sd]] <- list()
    k_list_per_run <- lapply(noise_stable_clusters[[noise_sd]], function(x) names(x))
    print(paste0("Noise SD: ", noise_sd))

    k_freq <- table(unlist(k_list_per_run))
    selected_k <- names(k_freq)[k_freq >= min(considered_k_app_threshold, max(k_freq) %/% 2)]
    selected_k <- intersect(selected_k, names(baseline_stable_clusters))
    print(length(selected_k))

    for (k in selected_k) {
        available_clusters <- lapply(noise_stable_clusters[[noise_sd]], function(x) {
            if (k %in% names(x)) {
                return(x[[k]]$partitions[[1]]$mb)
            } else {
                return(NULL)
            }
        })
        available_clusters <- available_clusters[!sapply(available_clusters, is.null)]
        noise_against_gt_ecc[[noise_sd]][[k]] <- element_agreement(
            reference_clustering = baseline_stable_clusters[[k]]$partitions[[1]]$mb,
            clustering_list = available_clusters
        )
        noise_gt_df <- rbind(
            noise_gt_df,
            data.frame(
                noise_sd = noise_sd,
                k = k,
                ecc = noise_against_gt_ecc[[noise_sd]][[k]]
            )
        )
    }
}
noise_gt_df$k <- factor(noise_gt_df$k, levels = stringr::str_sort(unique(noise_gt_df$k), numeric = TRUE))
print(noise_gt_df %>% group_by(noise_sd) %>% summarise(avg = mean(ecc), med = median(ecc)))

pdf("noise_comparison_to_ground_truth_per_k.pdf", width = 5, height = 7.5)
ggplot(noise_gt_df, aes(x = k, y = ecc)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~noise_sd, ncol = 2) +
    theme_bw() +
    labs(
        title = "Element agreement of stable clusters against ground truth across noise runs",
        x = "Noise SD",
        y = "Element agreement"
    )
dev.off()

pdf("noise_overall_comparison_to_ground_truth.pdf", width = 6, height = 6)
ggplot(noise_gt_df, aes(x = noise_sd, y = ecc)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    labs(
        title = "Element agreement of stable clusters against ground truth across noise runs",
        x = "Noise SD",
        y = "Element agreement"
    )
dev.off()

# reference
so <- qs2::qs_read("../data/masld_immune_filtered_normalized.qs2", nthreads = 4)
per_gene_sd <- apply(Seurat::GetAssayData(so, assay = "RNA", layer = "data"), 1, sd)
quantile(per_gene_sd, probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99))
