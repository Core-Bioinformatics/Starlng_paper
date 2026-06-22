library(dplyr)
library(ggplot2)
library(Starlng)
library(patchwork)
library(Seurat)
library(grid)

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

trend_plot_list <- list()
heatmap_plot_list <- list()

as_patchwork_plot <- function(p) {
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
}

for (dts_name in names(module_stats)) {
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


        # trend_plot_list[[dts_name]][[method_name]] <- plot_module_trends_over_pseudotime(
        #     expression_list = method_avg,
        #     pseudotime = pseudotime_vals
        # ) + theme(legend.position = "none") + ggtitle(paste0(dts_name, " - ", method_name))


        used_metadata <- data.frame(
            pseudotime = pseudotime_vals,
            row.names = colnames(expr_matrix)
        )
        if (dts_name == "masld_immune") {
            used_metadata$cell_clusters <- so$"stable_24_clusters"
        } else {
            used_metadata$cell_clusters <- so$"seurat_clusters"
        }
        heatmap_plot_list[[dts_name]][[method_name]] <- generate_cell_heatmap(
            expression_matrix = expr_matrix[method_clust_df$gene, ],
            gene_family_list = split(method_clust_df$gene, method_clust_df$module)[selected_modules],
            metadata_df = used_metadata,
            metadata_name = "cell_clusters",
            cap = 2
        )
    }
}

# pdf(file.path(output_dir, "module_trends_over_pseudotime.pdf"), width = 20, height = 15)
# # each dts should be a column
# column_list <- list()
# for (dts_name in names(trend_plot_list)) {
#     p <- trend_plot_list[[dts_name]][[1]]
#     for (method_name in names(trend_plot_list[[dts_name]])[-1]) {
#         if (method_name == "starlng_low_moran") {
#             next
#         }
#         p <- p / trend_plot_list[[dts_name]][[method_name]]
#     }
#     column_list[[dts_name]] <- p
# }
# wrap_plots(column_list, ncol = length(column_list))
# dev.off()

pdf(file.path(output_dir, "module_trends_heatmaps.pdf"), width = 20, height = 15)
column_list <- list()
for (dts_name in names(heatmap_plot_list)) {
    method_plot_list <- heatmap_plot_list[[dts_name]]
    if (dts_name == "masld_immune") {
        kept_names <- setdiff(names(heatmap_plot_list[[dts_name]]), "starlng_low_moran")
        method_plot_list <- method_plot_list[kept_names]
    }

    converted_method_plots <- lapply(method_plot_list, function(p) {
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
    column_list[[dts_name]] <- wrap_plots(converted_method_plots, ncol = 1)
}
wrap_plots(column_list, ncol = length(column_list))
dev.off()
