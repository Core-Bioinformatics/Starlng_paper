library(ggplot2)
library(dplyr)
library(Seurat)
library(Starlng)
devtools::load_all("/mnt/d/Starlng")

if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}

# dts_names <- c(paste0("cao_", c("Liver", "Pancreas", "Lung")), "masld_immune")
dts_names <- c("cao_Liver_subset")
method_names <- c("starlng", "hotspot", "hdwgcna", "paga_cellrank", "scenic")
output_dir <- "comparison_files"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

output_file_avg_summaries <- file.path(output_dir, "3_avg_summaries.qs2")
output_file_module_stats_summary <- file.path(output_dir, "3_module_stats_summary.qs2")
output_file_outlier_coverage_summary <- file.path(output_dir, "3_outlier_coverage_summary.qs2")

if (file.exists(output_file_avg_summaries)) {
    avg_summaries_list <- qs2::qs_read(output_file_avg_summaries, nthreads = 30)
} else {
    avg_summaries_list <- list()
}

if (file.exists(output_file_module_stats_summary)) {
    module_stats_summary_list <- qs2::qs_read(output_file_module_stats_summary, nthreads = 30)
} else {
    module_stats_summary_list <- list()
}
if (file.exists(output_file_outlier_coverage_summary)) {
    outlier_coverage_summary_list <- qs2::qs_read(output_file_outlier_coverage_summary, nthreads = 30)
} else {
    outlier_coverage_summary_list <- list()
}

for (dts_name in dts_names) {
    if (!dts_name %in% names(avg_summaries_list)) {
        avg_summaries_list[[dts_name]] <- list()
        module_stats_summary_list[[dts_name]] <- list()
        outlier_coverage_summary_list[[dts_name]] <- list()
    }

    so <- qs2::qs_read(file.path("..", "data", paste0(dts_name, "_filtered_normalized.qs2")), nthreads = 30)
    expr_matrix <- GetAssayData(so, assay = "RNA", layer = "data")
    umap_emb <- Embeddings(so, reduction = "umap")
    rm(so)
    gc()
    pseudotime_details <- read.csv(file.path("..", "data", paste0(dts_name, "_recommended_pseudotime.csv")))
    umap_median_dist <- pseudotime_details$umap_median_dist[1]
    pseudotime_value <- pseudotime_details$recommended_pseudotime
    names(pseudotime_value) <- pseudotime_details$X
    psd_mask <- !is.na(pseudotime_value)

    if (dts_name == "masld_immune") {
        actual_methods <- c(method_names[-1], "starlng_high_moran", "starlng_low_moran")
    } else {
        actual_methods <- method_names
    }
    for (method_name in actual_methods) {
        if (startsWith(method_name, "starlng") && dts_name == "masld_immune") {
            prefix <- gsub("starlng_", "", method_name)
            input_file <- file.path(paste0("starlng", "_results"), paste0(dts_name, "_", prefix, "_processed.csv"))
        } else {
            input_file <- file.path(paste0(method_name, "_results"), paste0(dts_name, "_processed.csv"))
        }

        if (!file.exists(input_file)) {
            next
        }
        print(paste0("Processing ", dts_name, " with method ", method_name))
        # if (method_name %in% names(avg_summaries_list[[dts_name]])) {
            # next
        # }

        modules_df <- read.csv(input_file)
        modules_df$module <- as.character(modules_df$module)
        module_names <- stringr::str_sort(unique(modules_df$module), numeric = TRUE)
        module_avg_scaled_expr <- lapply(module_names, function(module) {
            module_genes <- modules_df$gene[modules_df$module == module]
            module_expr <- expr_matrix[module_genes, , drop = FALSE]

            expr_value <- voting_scheme(
                expression_matrix = module_expr,
                thresh_percentile = 0,
                thresh_value = 0,
                n_coexpressed_thresh = 1,
                summary_function = mean
            )
            (expr_value - min(expr_value)) / (max(expr_value) - min(expr_value))
        })
        names(module_avg_scaled_expr) <- module_names
        avg_summaries_list[[dts_name]][[method_name]] <- module_avg_scaled_expr

        module_mask <- build_module_masks(
            module_summ = module_avg_scaled_expr,
            psd_mask = psd_mask,
            scale_threshold = 0.5,
            top_cells_percent = 100
        )
        rownames(module_mask) <- names(module_avg_scaled_expr)

        modules_stats <- get_module_stats(
            module_summ = module_avg_scaled_expr,
            module_mask = module_mask,
            psd_value = pseudotime_value,
            umap_df = umap_emb,
            centroid = TRUE
        )

        modules_stats_summary <- summarise_module_stats(
            modules_stats = modules_stats,
            gene_modules = split(modules_df$gene, modules_df$module)
        )

        outlier_assessment <- detect_outlier(
            modules_stats = modules_stats_summary,
            cell_masks = module_mask,
            psd_value = pseudotime_value,
            umap_dist_threshold = umap_median_dist * 0.85
        )

        
        modules_stats_summary$is_outlier <- outlier_assessment$outlier_output[rownames(modules_stats_summary)]
        module_stats_summary_list[[dts_name]][[method_name]] <- modules_stats_summary
        outlier_coverage_summary_list[[dts_name]][[method_name]] <- outlier_assessment$coverage_evolution_df


        filtered_module_mask <- module_mask[modules_stats_summary %>% filter(is_outlier == "no") %>% pull(module), , drop = FALSE]
        qs2::qs_save(avg_summaries_list, output_file_avg_summaries, nthreads = 5)
        qs2::qs_save(module_stats_summary_list, output_file_module_stats_summary, nthreads = 5)
        qs2::qs_save(outlier_coverage_summary_list, output_file_outlier_coverage_summary, nthreads = 5)
    }
}

qs2::qs_save(avg_summaries_list, output_file_avg_summaries, nthreads = 5)
qs2::qs_save(module_stats_summary_list, output_file_module_stats_summary, nthreads = 5)
qs2::qs_save(outlier_coverage_summary_list, output_file_outlier_coverage_summary, nthreads = 5)
