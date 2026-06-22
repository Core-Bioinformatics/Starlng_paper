library(dplyr)
library(ggplot2)

if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}

input_dir <- "comparison_files"
output_dir <- "panels"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

outlier_coverage <- qs2::qs_read(file.path(input_dir, "3_outlier_coverage_summary.qs2"), nthreads = 30)
module_stats <- qs2::qs_read(file.path(input_dir, "3_module_stats_summary.qs2"), nthreads = 30)

bar_colours <- qualpalr::qualpal(3)$hex
names(bar_colours) <- c("selected", "outlier", "redundant")

method_colours <- qualpalr::qualpal(7)$hex
names(method_colours) <- c("Starlng", "Hotspot", "hdWGCNA", "SCENIC", "PAGA+\nCellRank", "Starlng\n(low Moran's I)", "Starlng\n(high Moran's I)")

for (dts_name in names(outlier_coverage)) {
    perc_output_file <- file.path(output_dir, paste0(dts_name, "_barplot_outliers.pdf"))
    coverage_output_file <- file.path(output_dir, paste0(dts_name, "_cumulative_coverage.pdf"))
    perc_df <- NULL
    coverage_df <- NULL
    for (method_name in names(outlier_coverage[[dts_name]])) {
        perc_df <- rbind(perc_df, data.frame(
            "method" = method_name,
            "module_type" = module_stats[[dts_name]][[method_name]]$is_outlier,
            "n_genes" = module_stats[[dts_name]][[method_name]]$n_genes
        ))
        coverage_df <- rbind(coverage_df, data.frame(
            "method" = method_name,
            "coverage" = outlier_coverage[[dts_name]][[method_name]]$coverage_added,
            "module" = outlier_coverage[[dts_name]][[method_name]]$added_module,
            "iqr" = outlier_coverage[[dts_name]][[method_name]]$module_iqr
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

    pdf(perc_output_file, width = 7, height = 5)
    print(ggplot(module_plot_df, aes(x = proportion, y = method, fill = module_type)) +
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
    )
    dev.off()

    # convert coverage to cumulative, ordered by iqr
    coverage_df <- coverage_df %>%
        group_by(method) %>%
        arrange(iqr) %>%
        mutate(cumulative_coverage = cumsum(coverage))

    pdf(coverage_output_file, width = 7, height = 5)
    print(ggplot(coverage_df, aes(x = iqr, y = cumulative_coverage, color = method)) +
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
    )
    dev.off()
}
