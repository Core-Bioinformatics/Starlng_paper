library(dplyr)
library(ggplot2)
library(ComplexUpset)
library(tidyr)
library(patchwork)

if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}

input_dir <- "comparison_files"
output_dir <- "panels"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

enrichment_all <- qs2::qs_read(file.path(input_dir, "4_enrichment_analysis.qs2"), nthreads = 30)
enrichment_top <- qs2::qs_read(file.path(input_dir, "4_enrichment_top_10.qs2"), nthreads = 30)
outlier_assessment <- qs2::qs_read(file.path(input_dir, "3_module_stats_summary.qs2"), nthreads = 30)
top_n_modules <- 5
term_size_threshold <- 500

comparison_table_file <- file.path("..", "comparison_table.csv")
if (!file.exists(comparison_table_file)) {
    comparison_table <- data.frame(
        method_name = character(),
        dataset_name = character(),
        comparison_type = character(),
        value = numeric(),
        rank = numeric()
    )
} else {
    comparison_table <- read.csv(comparison_table_file, stringsAsFactors = FALSE)
    # remove rows with comparison_type starting with "enrichment_all_modules_" or "enrichment_top10hubs_modules_"
    comparison_table <- comparison_table %>%
        filter(!(comparison_type %in% c(
            "enrichment_all_modules_GO:BP",
            "enrichment_all_modules_TF",
            "enrichment_top10hubs_modules_GO:BP",
            "enrichment_top10hubs_modules_TF",
            "selected_module_proportion",
            "selected_module_total_coverage",
            "selected_module_max_iqr"
        )))
}


method_colours <- qualpalr::qualpal(7)$hex
names(method_colours) <- c("Starlng", "Hotspot", "hdWGCNA", "SCENIC", "PAGA+\nCellRank", "Starlng\n(low Moran's I)", "Starlng\n(high Moran's I)")

build_upset_plot <- function(current_terms_df, dts_name, source_name, terms_label, n_intersections = 5) {
    upset_input <- current_terms_df %>%
        group_by(.data$method, .data$term_id) %>%
        summarise(in_set = 1L, .groups = "drop") %>%
        pivot_wider(names_from = "method", values_from = "in_set", values_fill = 0) %>%
        as.data.frame()

    if (!"term_id" %in% colnames(upset_input)) {
        message(sprintf("Skipping %s: not enough methods with %s for UpSet plot.", dts_name, terms_label))
        return(NULL)
    }

    method_cols <- setdiff(colnames(upset_input), "term_id")
    non_empty_method_cols <- method_cols[colSums(upset_input[, method_cols, drop = FALSE] > 0) > 0]
    if (length(non_empty_method_cols) < 2) {
        message(sprintf("Skipping %s: fewer than two methods with %s for UpSet plot.", dts_name, terms_label))
        return(NULL)
    }

    upset_input <- upset_input[, c("term_id", non_empty_method_cols), drop = FALSE]
    upset_input[, non_empty_method_cols] <- (upset_input[, non_empty_method_cols, drop = FALSE] > 0) * 1L

    rownames(upset_input) <- upset_input$term_id
    upset_input <- upset_input[, -1]

    if (ncol(upset_input) == 0 || nrow(upset_input) == 0) {
        return(list(
            plot = ggplot() + theme_void() + ggtitle(paste(dts_name, source_name)),
            include_source = FALSE
        ))
    }

    list(
        plot = ComplexUpset::upset(
            upset_input,
            intersect = colnames(upset_input),
            name = "Method",
            n_intersections = n_intersections
        ) + ggtitle(paste(dts_name, source_name)),
        include_source = TRUE
    )
}

upset_all <- list()
upset_top <- list()

upset_all_sources <- list()
upset_top_sources <- list()

for (dts_name in names(enrichment_all)) {
    if (dts_name %in% c("cao_Liver", "cao_Lung")) {
        next
    }
    upset_all[[dts_name]] <- list(
        ggplot() + theme_void() + ggtitle(paste(dts_name, "GO:BP")),
        ggplot() + theme_void() + ggtitle(paste(dts_name, "TF"))
    )
    upset_top[[dts_name]] <- list(
        ggplot() + theme_void() + ggtitle(paste(dts_name, "GO:BP")),
        ggplot() + theme_void() + ggtitle(paste(dts_name, "TF"))
    )
    upset_all_sources[[dts_name]] <- character(0)
    upset_top_sources[[dts_name]] <- character(0)

    unique_methods <- recode(names(enrichment_all[[dts_name]]), "starlng" = "Starlng", "hotspot" = "Hotspot", "hdwgcna" = "hdWGCNA", "scenic" = "SCENIC", "paga_cellrank" = "PAGA+\nCellRank", "starlng_low_moran" = "Starlng\n(low Moran's I)", "starlng_high_moran" = "Starlng\n(high Moran's I)")

    all_terms_df <- NULL
    top_terms_df <- NULL
    for (method_name in names(enrichment_all[[dts_name]])) {
        print(method_name)
        method_all_df <- NULL
        considered_modules <- outlier_assessment[[dts_name]][[method_name]] %>%
            filter(is_outlier == "no") %>%
            pull(module)
        enrich_method <- enrichment_all[[dts_name]][[method_name]]
        searched_modules <- intersect(considered_modules, names(enrich_method))
        results_per_module <- sapply(searched_modules, function(module_name) {
            module_enrich <- enrich_method[[module_name]]
            if (is.null(module_enrich)) {
                return(0)
            }
            module_enrich <- module_enrich$result
            if (is.null(module_enrich) || nrow(module_enrich) == 0) {
                return(0)
            }
            return(module_enrich %>% filter(source == "GO:BP", term_size <= term_size_threshold) %>% nrow())
        })

        if (length(searched_modules) == 0) {
            next
        }
        n_selected_modules <- min(top_n_modules, length(results_per_module))
        if (n_selected_modules == 0) {
            selected_modules <- character(0)
        } else {
            selected_modules <- names(sort(results_per_module, decreasing = TRUE))[seq_len(n_selected_modules)]
        }
        for (module_name in searched_modules) {
            print(paste("Module:", module_name))
            module_enrich <- enrich_method[[module_name]]
            if (is.null(module_enrich)) {
                next
            }
            module_enrich <- module_enrich$result

            if (is.null(module_enrich) || nrow(module_enrich) == 0) {
                next
            }
            print(nrow(module_enrich))
            module_enrich <- module_enrich %>% filter(term_size <= term_size_threshold)
            print(nrow(module_enrich))
            if (nrow(module_enrich) == 0) {
                next
            }
            module_enrich <- module_enrich[, c("term_id", "source")]
            module_enrich$method <- method_name
            module_enrich$module <- module_name
            method_all_df <- rbind(method_all_df, module_enrich)
        }

        # make unique term_id - method combinations
        method_all_df <- method_all_df %>%
            distinct(term_id, method, .keep_all = TRUE)
        all_terms_df <- rbind(all_terms_df, method_all_df)
    }
    all_terms_df$method <- recode(all_terms_df$method, "starlng" = "Starlng", "hotspot" = "Hotspot", "hdwgcna" = "hdWGCNA", "scenic" = "SCENIC", "paga_cellrank" = "PAGA+\nCellRank", "starlng_low_moran" = "Starlng\n(low Moran's I)", "starlng_high_moran" = "Starlng\n(high Moran's I)")

    if (isFALSE(is.null(all_terms_df) || nrow(all_terms_df) == 0)) {
        unique_sources <- unique(all_terms_df$source)
        index <- 1
        for (source_name in unique_sources) {
            current_terms_df <- all_terms_df %>%
                filter(source == source_name) %>%
                group_by(method, term_id) %>%
                ungroup()

            upset_result <- build_upset_plot(
                current_terms_df = current_terms_df,
                dts_name = dts_name,
                source_name = source_name,
                terms_label = "enrichment terms"
            )
            if (is.null(upset_result)) {
                next
            }

            upset_all[[dts_name]][[index]] <- upset_result$plot
            if (upset_result$include_source) {
                upset_all_sources[[dts_name]][index] <- source_name
            }
            index <- index + 1

            method_ranking <- current_terms_df %>%
                group_by(method) %>%
                summarise(value = n_distinct(term_id), .groups = "drop")
            colnames(method_ranking) <- c("method_name", "value")
            method_ranking$rank <- rank(-method_ranking$value, ties.method = "min")
            for (missing_method in setdiff(unique_methods, method_ranking$method_name)) {
                method_ranking <- rbind(method_ranking, data.frame(method_name = missing_method, value = 0, rank = length(unique_methods)))
            }
            method_ranking$dataset_name <- dts_name
            method_ranking$comparison_type <- paste0("enrichment_all_modules_", source_name)
            comparison_table <- rbind(comparison_table, method_ranking[, colnames(comparison_table)])
        }
    }




    # for (method_name in names(enrichment_top[[dts_name]])) {
    #     method_df <- NULL
    #     top_method <- enrichment_top[[dts_name]][[method_name]]
    #     if (is.null(top_method) || length(top_method) == 0) {
    #         next
    #     }
    #     considered_modules <- outlier_assessment[[dts_name]][[method_name]] %>%
    #         filter(is_outlier == "no") %>%
    #         pull(module)
    #     searched_modules <- intersect(considered_modules, names(top_method))
    #     if (length(searched_modules) == 0) {
    #         next
    #     }
    #     for (module_name in names(top_method)) {
    #         module_enrich <- top_method[[module_name]]
    #         if (is.null(module_enrich)) {
    #             next
    #         }
    #         module_enrich <- module_enrich$result

    #         if (is.null(module_enrich) || nrow(module_enrich) == 0) {
    #             next
    #         }
    #         module_enrich <- module_enrich[, c("term_id", "source")]
    #         module_enrich$method <- method_name
    #         module_enrich$module <- module_name
    #         method_df <- rbind(method_df, module_enrich)
    #     }
    #     method_df <- method_df %>%
    #         distinct(term_id, method, .keep_all = TRUE)
    #     top_terms_df <- rbind(top_terms_df, method_df)
    # }

    # top_terms_df$method <- recode(top_terms_df$method, "starlng" = "Starlng", "hotspot" = "Hotspot", "hdwgcna" = "hdWGCNA", "scenic" = "SCENIC", "paga_cellrank" = "PAGA+\nCellRank", "starlng_low_moran" = "Starlng\n(low Moran's I)", "starlng_high_moran" = "Starlng\n(high Moran's I)")

    # if (is.null(top_terms_df) || nrow(top_terms_df) == 0) {
    #     message(sprintf("Skipping %s: no top enrichment terms found.", dts_name))
    #     next
    # }

    # unique_sources_top <- unique(top_terms_df$source)

    # # pdf(file.path(output_dir, paste0(dts_name, "_enrichment_top_upset.pdf")), width = 7, height = 5)
    # for (source_name in unique_sources_top) {
    #     current_terms_df <- top_terms_df %>%
    #         filter(source == source_name) %>%
    #         group_by(method, term_id) %>%
    #         ungroup()

    #     upset_input <- current_terms_df %>%
    #         group_by(method, term_id) %>%
    #         summarise(in_set = 1L, .groups = "drop") %>%
    #         pivot_wider(names_from = method, values_from = in_set, values_fill = 0) %>%
    #         as.data.frame()

    #     if (!"term_id" %in% colnames(upset_input)) {
    #         message(sprintf("Skipping %s: not enough methods with top enrichment terms for UpSet plot.", dts_name))
    #         next
    #     }

    #     method_cols <- setdiff(colnames(upset_input), "term_id")
    #     non_empty_method_cols <- method_cols[colSums(upset_input[, method_cols, drop = FALSE] > 0) > 0]
    #     if (length(non_empty_method_cols) < 2) {
    #         message(sprintf("Skipping %s: fewer than two methods with top enrichment terms for UpSet plot.", dts_name))
    #         next
    #     }

    #     upset_input <- upset_input[, c("term_id", non_empty_method_cols), drop = FALSE]
    #     upset_input[, non_empty_method_cols] <- (upset_input[, non_empty_method_cols, drop = FALSE] > 0) * 1L

    #     rownames(upset_input) <- upset_input$term_id
    #     upset_input <- upset_input[, -1]

    #     if (ncol(upset_input) == 0 || nrow(upset_input) == 0) {
    #         message(sprintf("Skipping %s: empty UpSet input after preprocessing for top terms.", dts_name))
    #         next
    #     }

    #     print(
    #         upset(upset_input, order.by = "freq", nsets = length(enrichment_top[[dts_name]])) 
    #     )
    # }
    # dev.off()


    top5_modules_terms_df <- NULL
    for (method_name in names(enrichment_top[[dts_name]])) {
        method_df <- NULL
        top_method <- enrichment_top[[dts_name]][[method_name]]
        if (is.null(top_method) || length(top_method) == 0) {
            next
        }

        considered_modules <- outlier_assessment[[dts_name]][[method_name]] %>%
            filter(is_outlier == "no") %>%
            pull(module)
        searched_modules <- intersect(considered_modules, names(top_method))
        results_per_module <- sapply(searched_modules, function(module_name) {
            module_enrich <- top_method[[module_name]]
            if (is.null(module_enrich)) {
                return(0)
            }
            module_enrich <- module_enrich$result
            if (is.null(module_enrich) || nrow(module_enrich) == 0) {
                return(0)
            }
            return(module_enrich %>% filter(source == "GO:BP", term_size <= term_size_threshold) %>% nrow())
        })
        if (length(searched_modules) == 0) {
            next
        }

        n_selected_modules <- min(top_n_modules, length(results_per_module))
        if (n_selected_modules == 0) {
            selected_modules <- character(0)
        } else {
            selected_modules <- names(sort(results_per_module, decreasing = TRUE))[seq_len(n_selected_modules)]
        }
        for (module_name in selected_modules) {
            module_enrich <- top_method[[module_name]]
            if (is.null(module_enrich)) {
                next
            }
            module_enrich <- module_enrich$result

            if (is.null(module_enrich) || nrow(module_enrich) == 0) {
                next
            }
            module_enrich <- module_enrich %>% filter(term_size <= term_size_threshold)
            if (nrow(module_enrich) == 0) {
                next
            }
            module_enrich <- module_enrich[, c("term_id", "source")]
            module_enrich$method <- method_name
            module_enrich$module <- module_name
            method_df <- rbind(method_df, module_enrich)
        }

        method_df <- method_df %>%
            distinct(term_id, method, .keep_all = TRUE)
        top5_modules_terms_df <- rbind(top5_modules_terms_df, method_df)
    }

    if (is.null(top5_modules_terms_df) || nrow(top5_modules_terms_df) == 0) {
        message(sprintf("Skipping %s: no top enrichment terms found for top 5 modules per method.", dts_name))
        next
    }

    top5_modules_terms_df$method <- recode(top5_modules_terms_df$method, "starlng" = "Starlng", "hotspot" = "Hotspot", "hdwgcna" = "hdWGCNA", "scenic" = "SCENIC", "paga_cellrank" = "PAGA+\nCellRank", "starlng_low_moran" = "Starlng\n(low Moran's I)", "starlng_high_moran" = "Starlng\n(high Moran's I)")

    unique_sources_top5 <- unique(top5_modules_terms_df$source)

    # pdf(file.path(output_dir, paste0(dts_name, "_enrichment_top_upset_top5_modules.pdf")), width = 7, height = 5)
    index <- 1
    for (source_name in unique_sources_top5) {
        current_terms_df <- top5_modules_terms_df %>%
            filter(source == source_name) %>%
            group_by(method, term_id) %>%
            ungroup()

        upset_result <- build_upset_plot(
            current_terms_df = current_terms_df,
            dts_name = dts_name,
            source_name = source_name,
            terms_label = "top enrichment terms for top-5-modules"
        )
        if (is.null(upset_result)) {
            next
        }

        upset_top[[dts_name]][[index]] <- upset_result$plot
        if (upset_result$include_source) {
            upset_top_sources[[dts_name]][index] <- source_name
        }
        index <- index + 1

        # create a ranking of methods based on the number of unique terms in the top 5 modules for this source
        method_ranking <- current_terms_df %>%
            group_by(method) %>%
            summarise(value = n_distinct(term_id), .groups = "drop")
        colnames(method_ranking) <- c("method_name", "value")
        method_ranking$rank <- rank(-method_ranking$value, ties.method = "min")
        for (missing_method in setdiff(unique_methods, method_ranking$method_name)) {
            method_ranking <- rbind(method_ranking, data.frame(method_name = missing_method, value = 0, rank = length(unique_methods)))
        }
        method_ranking$dataset_name <- dts_name
        method_ranking$comparison_type <- paste0("enrichment_top10hubs_modules_", source_name)
        comparison_table <- rbind(comparison_table, method_ranking[, colnames(comparison_table)])
    }

}

rows_all <- lapply(upset_all, function(x) {
    wrap_plots(x, nrow = 1)
})
columns_all <- wrap_plots(rows_all, ncol = 1)

rows_top <- lapply(upset_top, function(x) {
    wrap_plots(x, nrow = 1)
})
columns_top <- wrap_plots(rows_top, ncol = 1)

row_enrich <- wrap_plots(columns_all, columns_top, nrow = 1)



outlier_coverage_data <- qs2::qs_read(file.path(input_dir, "3_outlier_coverage_summary.qs2"), nthreads = 30)
module_stats <- qs2::qs_read(file.path(input_dir, "3_module_stats_summary.qs2"), nthreads = 30)

bar_colours <- qualpalr::qualpal(3)$hex
names(bar_colours) <- c("selected", "outlier", "redundant")

method_colours <- qualpalr::qualpal(7)$hex
names(method_colours) <- c("Starlng", "Hotspot", "hdWGCNA", "SCENIC", "PAGA+\nCellRank", "Starlng\n(low Moran's I)", "Starlng\n(high Moran's I)")

outlier_perc <- list()
outlier_coverage <- list()
for (dts_name in names(outlier_coverage_data)) {
    if (dts_name %in% c("cao_Liver", "cao_Lung")) {
        next
    }
    perc_output_file <- file.path(output_dir, paste0(dts_name, "_barplot_outliers.pdf"))
    coverage_output_file <- file.path(output_dir, paste0(dts_name, "_cumulative_coverage.pdf"))
    perc_df <- NULL
    coverage_df <- NULL
    for (method_name in names(outlier_coverage_data[[dts_name]])) {
        perc_df <- rbind(perc_df, data.frame(
            "method" = method_name,
            "module_type" = module_stats[[dts_name]][[method_name]]$is_outlier,
            "n_genes" = module_stats[[dts_name]][[method_name]]$n_genes
        ))
        coverage_df <- rbind(coverage_df, data.frame(
            "method" = method_name,
            "coverage" = outlier_coverage_data[[dts_name]][[method_name]]$coverage_added,
            "module" = outlier_coverage_data[[dts_name]][[method_name]]$added_module,
            "iqr" = outlier_coverage_data[[dts_name]][[method_name]]$module_iqr
        ))
    }
    # remap no to - selected, yes to - outlier, redundant to - redundant
    perc_df$method <- recode(perc_df$method, "starlng" = "Starlng", "hotspot" = "Hotspot", "hdwgcna" = "hdWGCNA", "scenic" = "SCENIC", "paga_cellrank" = "PAGA+\nCellRank", "starlng_low_moran" = "Starlng\n(low Moran's I)", "starlng_high_moran" = "Starlng\n(high Moran's I)")
    coverage_df$method <- recode(coverage_df$method, "starlng" = "Starlng", "hotspot" = "Hotspot", "hdwgcna" = "hdWGCNA", "scenic" = "SCENIC", "paga_cellrank" = "PAGA+\nCellRank", "starlng_low_moran" = "Starlng\n(low Moran's I)", "starlng_high_moran" = "Starlng\n(high Moran's I)")
    perc_df$module_type <- recode(perc_df$module_type, "no" = "selected", "yes" = "outlier", "redundant" = "redundant")
    perc_df$module_type <- factor(perc_df$module_type, levels = c("outlier", "redundant", "selected"))

    summarised_perc_df <- perc_df %>%
        group_by(method) %>%
        summarise(
            n_modules = n(),
            .groups = "drop"
        )

    module_plot_df <- perc_df %>%
        group_by(method, module_type) %>%
        summarise(
            n_modules_type = n(),
            genes_covered = sum(n_genes, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        group_by(method) %>%
        arrange(module_type, .by_group = TRUE) %>%
        mutate(
            proportion = n_modules_type / sum(n_modules_type),
            label = paste0("n=", format(genes_covered, big.mark = ",", scientific = FALSE, trim = TRUE))
        ) %>%
        ungroup() %>%
        as.data.frame()
    module_plot_df$module_type <- factor(module_plot_df$module_type, levels = c("outlier", "redundant", "selected"))
    module_plot_df <- module_plot_df %>%
        group_by(method) %>%
        arrange(desc(module_type), .by_group = TRUE) %>%
        mutate(
            x_mid = cumsum(proportion) - 0.5 * proportion
        ) %>%
        ungroup() %>%
        as.data.frame()
    print(module_plot_df)

    outlier_perc[[dts_name]] <-  ggplot(module_plot_df, aes(x = proportion, y = method, fill = module_type)) +
        geom_col() +
        ggrepel::geom_label_repel(
            data = module_plot_df,
            aes(x = x_mid, y = method, label = label),
            inherit.aes = FALSE,
            size = 3.1,
            fontface = "italic",
            nudge_y = 0.28
        ) +
        scale_fill_manual(values = bar_colours ) +
        labs(x = "Proportion of Modules", y = "Method", fill = "Module Type") +
        theme_classic() +
        theme(
            axis.text.x = element_text(),
            legend.position = "bottom"
        ) +
        # add the number of modules at the right of the bars
        geom_text(data = summarised_perc_df, aes(x = 1.05, y = method, label = n_modules), size = 5, inherit.aes = FALSE, fontface = "bold") +
        scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1.1)) +
        labs(
            title = paste0(dts_name, "\nProportion of Outlier Modules by Method")
        ) +
        theme(
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            title = element_text(size = 16),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14)
        )
    
    # add to overall results
    method_ranking <- module_plot_df %>%
        group_by(method) %>%
        summarise(percent_selected = proportion[module_type == "selected"], .groups = "drop")
    colnames(method_ranking) <- c("method_name", "value")
    method_ranking$rank <- rank(-method_ranking$value, ties.method = "min")
    for (missing_method in setdiff(unique_methods, method_ranking$method_name)) {
        method_ranking <- rbind(method_ranking, data.frame(method_name = missing_method, value = 0, rank = length(unique_methods)))
    }
    method_ranking$dataset_name <- dts_name
    method_ranking$comparison_type <- "selected_module_proportion"
    comparison_table <- rbind(comparison_table, method_ranking[, colnames(comparison_table)])
    


    # convert coverage to cumulative, ordered by iqr
    coverage_df <- coverage_df %>%
        group_by(method) %>%
        arrange(iqr) %>%
        mutate(cumulative_coverage = cumsum(coverage))

    outlier_coverage[[dts_name]] <- ggplot(coverage_df, aes(x = iqr, y = cumulative_coverage, color = method)) +
        geom_line(size = 1) +
        geom_point(size = 2) +
        labs(x = "Module IQR", y = "Cumulative Coverage", color = "Method") +
        theme_classic() +
        guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
        scale_color_manual(values = method_colours) +
        labs(
            title = paste0(dts_name, "\nCumulative coverage of outlier modules by method")
        ) +
        theme(
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            title = element_text(size = 16),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            # position legend inside in the right side
            legend.position = c(0.8, 0.2),
        )
    # add to overall results - total coverage
    method_ranking <- coverage_df %>%
        group_by(method) %>%
        summarise(total_coverage = max(cumulative_coverage), .groups = "drop")
    colnames(method_ranking) <- c("method_name", "value")
    method_ranking$rank <- rank(-method_ranking$value, ties.method = "min")
    for (missing_method in setdiff(unique_methods, method_ranking$method_name)) {
        method_ranking <- rbind(method_ranking, data.frame(method_name = missing_method, value = 0, rank = length(unique_methods)))
    }
    method_ranking$dataset_name <- dts_name
    method_ranking$comparison_type <- "selected_module_total_coverage"

    comparison_table <- rbind(comparison_table, method_ranking[, colnames(comparison_table)])

    # add to overall results - min max IQR
    method_ranking <- coverage_df %>%
        group_by(method) %>%
        summarise(max_iqr = max(iqr), .groups = "drop")
    colnames(method_ranking) <- c("method_name", "value")
    method_ranking$rank <- rank(method_ranking$value, ties.method = "min")
    for (missing_method in setdiff(unique_methods, method_ranking$method_name)) {
        method_ranking <- rbind(method_ranking, data.frame(method_name = missing_method, value = NA, rank = length(unique_methods)))
    }
    method_ranking$dataset_name <- dts_name
    method_ranking$comparison_type <- "selected_module_max_iqr"

    comparison_table <- rbind(comparison_table, method_ranking[, colnames(comparison_table)])
}




row_perc <- wrap_plots(outlier_perc, nrow = 1)
row_coverage <- wrap_plots(outlier_coverage, nrow = 1)


final_plot <- wrap_plots(row_perc, row_coverage, row_enrich, nrow = 3, heights = c(1, 1, 4))




pdf(file.path(output_dir, "sup_enrichment_comparison_upset_outlier.pdf"), width = 15, height = 21)
print(final_plot)
dev.off()

# replace \n with "" in method_name
comparison_table$method_name <- gsub("\n", "", comparison_table$method_name)
write.csv(comparison_table, comparison_table_file, row.names = FALSE)
