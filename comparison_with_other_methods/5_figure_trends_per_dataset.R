library(dplyr)
library(ggplot2)
library(Starlng)
library(patchwork)
library(Seurat)
library(grid)
devtools::load_all("/mnt/d/Starlng")

if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}

input_dir <- "comparison_files"
output_dir <- "panels"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

module_stats <- qs2::qs_read(file.path(input_dir, "3_module_stats_summary.qs2"), nthreads = 30)
avg_summs <- qs2::qs_read(file.path(input_dir, "3_avg_summaries.qs2"), nthreads = 30)


for (dts_name in names(module_stats)) {
    if (dts_name != "cao_Liver_subset") {
        next
    }
    trend_plot_list <- list()
    heatmap_plot_list <- list()

    so <- qs2::qs_read(file.path("..", "data", paste0(dts_name, "_filtered_normalized.qs2")), nthreads = 30)
    expr_matrix <- GetAssayData(so, assay = "RNA", layer = "data")
    pseudotime_vals <- read.csv(file.path("..", "data", paste0(dts_name, "_recommended_pseudotime.csv")))
    pseudotime_vals <- pseudotime_vals$recommended_pseudotime
    for (method_name in names(module_stats[[dts_name]])) {
        print(paste0("Processing ", dts_name, " - ", method_name))
        if (method_name %in% trend_plot_list[[dts_name]]) {
            next
        }
        if (startsWith(method_name, "starlng") && dts_name == "masld_immune") {
            prefix <- gsub("starlng_", "", method_name)
            input_file <- file.path(paste0("starlng", "_results"), paste0(dts_name, "_", prefix, "_processed.csv"))
        } else {
            input_file <- file.path(paste0(method_name, "_results"), paste0(dts_name, "_processed.csv"))
        }
        method_stats <- module_stats[[dts_name]][[method_name]]
        selected_modules <- method_stats %>%
            filter(is_outlier == "no") %>%
            arrange(median_pseudotime) %>%
            pull(module)
        method_avg <- avg_summs[[dts_name]][[method_name]][selected_modules]
        method_clust_df <- read.csv(input_file)
        method_clust_df <- method_clust_df %>% filter(module %in% selected_modules)


        trend_plot_list[[method_name]] <- plot_module_trends_over_pseudotime(
            expression_list = method_avg,
            pseudotime = pseudotime_vals,
            show_labels = FALSE,
            smooth_threshold_max = 0.2
        ) + ggtitle(paste0(dts_name, " - ", method_name, ": Expression trends over pseudotime")) +
        theme(
            legend.position = "none",
            plot.title = element_text(size = 20),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15)
        )

        blue_red_palette <- colorRampPalette(c("blue", "white", "red"))

        used_metadata <- data.frame(
            pseudotime = pseudotime_vals,
            row.names = colnames(expr_matrix)
        )
        if (dts_name == "masld_immune") {
            used_metadata$cell_clusters <- so$"stable_24_clusters"
        } else {
            used_metadata$cell_clusters <- so$"seurat_clusters"
        }
        heatmap_plot_list[[method_name]] <- invisible(ComplexHeatmap::draw(generate_cell_heatmap(
            expression_matrix = expr_matrix[method_clust_df$gene, ],
            gene_family_list = split(method_clust_df$gene, method_clust_df$module)[selected_modules],
            metadata_df = used_metadata,
            metadata_name = "cell_clusters",
            cap = 2,
            continuous_colors = blue_red_palette(100)
        ), merge_legends = TRUE))
    }

    a <- wrap_plots(trend_plot_list) + plot_layout(nrow = length(trend_plot_list))
    heatmap_grob <- lapply(heatmap_plot_list, function(p) {
        if (inherits(p, "gg") || inherits(p, "ggplot") || inherits(p, "patchwork")) {
            return(p)
        }
        if (inherits(p, "Heatmap") || inherits(p, "HeatmapList")) {
            heatmap_grob <- grid.grabExpr(ComplexHeatmap::draw(p))
            return(wrap_elements(full = heatmap_grob))
        }
        if (!is.null(p$gtable)) {
            return(wrap_elements(full = p$gtable))
        }
        if (inherits(p, "grob") || inherits(p, "gTree")) {
            return(wrap_elements(full = p))
        }
        stop("Unsupported heatmap plot object class: ", paste(class(p), collapse = ", "))
    })
    b <- wrap_plots(heatmap_grob) + plot_layout(nrow = length(heatmap_grob))

    first_row <- wrap_plots(a, b, nrow = 1, widths = c(1, 1.5))

    pdf(paste0("panels/sup_trends_", dts_name, ".pdf"), width = 13, height = 20)
    print(first_row)
    dev.off()

}


