library(dplyr)
library(ggplot2)
library(tidyr)
library(Starlng)
library(patchwork)

if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}

input_dir <- "comparison_files"
output_dir <- "panels"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}


method_colours <- qualpalr::qualpal(7)$hex
dts_names <- c("masld_immune", "cao_Liver_subset", "cao_Pancreas")
tf_bubbleplot_list <- list()
tf_graph_list <- list()

for (dts_name in dts_names) {
    if (dts_name == dts_names[1]) {
        app_dir <- file.path("..", "starlng_run", "masld_immune_high_moran_starlng_app")
    } else {
        app_dir <- file.path("..", "starlng_run", paste0(dts_name, "_starlng_app"))
    }
    app_dir <- file.path(app_dir, "objects")
    sum_path <- file.path(app_dir, "module_summaries.h5")
    emb_path <- file.path(app_dir, "gene_embedding.qs2")
    trajectory_object <- qs2::qs_read(file.path(app_dir, "trajectory_object.qs2"), nthreads = 30)
    cell_umap <- qs2::qs_read(file.path(app_dir, "metadata.qs2"), nthreads = 30)
    cell_umap <- cell_umap[, c(ncol(cell_umap) - 1, ncol(cell_umap))]

    emb <- qs2::qs_read(emb_path, nthreads = 30)
    gene_umap <- emb$umap
    gene_adj_matrix <- emb$adj_matrix

    module_summs <- rhdf5::h5read(sum_path, "/")
    module_names <- rhdf5::h5read(sum_path, "all_modules") %>% as.character()
    selected_module <- as.character(max(as.integer(module_names)))
    module_list <- split(rhdf5::h5read(sum_path, "genes") %>% as.character(), rhdf5::h5read(sum_path, paste0(selected_module, "/clustering")))
    expr_module <- rhdf5::h5read(sum_path, paste0(selected_module, "/expression_summaries"))
    expr_module <- split(expr_module, col(expr_module))
    module_medians <- get_module_centroid(
        module_expr = expr_module,
        cell_umap = cell_umap,
        expression_threshold = 0.5,
        scale = TRUE
    )
    closest_nodes_to_module <- setNames(
        rhdf5::h5read(sum_path, paste0(selected_module, "/closest_nodes_to_module")),
        rhdf5::h5read(sum_path, paste0(selected_module, "/modules"))
    )
    module_stats <- rhdf5::h5read(sum_path, paste0(selected_module, "/modules_stats_summary")) %>%
        filter(is_outlier == "no") %>%
        arrange(median_pseudotime)
    filtered_modules <- module_stats$module
    gene_hub_stats <- rhdf5::h5read(sum_path, paste0(selected_module, "/gene_hub_stats"))
    gene_hub_stats$gene <- rhdf5::h5read(sum_path, "genes") %>% as.character()
    tfs <- rhdf5::h5read(sum_path, paste0(selected_module, "/tfs"))

    ## TF
    tf_modules <- intersect(filtered_modules, unique(tfs$module))
    if (length(tf_modules) > 3) {
        n_top_hubs <- 3
        n_top_tfs <- 2
    } else {
        n_top_hubs <- 10
        n_top_tfs <- 4
    }

    gene_hub_stats <- rhdf5::h5read(sum_path, paste0(selected_module, "/gene_hub_stats"))
    gene_hub_stats$gene <- rhdf5::h5read(sum_path, "genes") %>% as.character()
    gene_hub_stats <- gene_hub_stats %>%
        group_by(module) %>%
        slice_max(combined_score, n = n_top_hubs) %>%
        filter(module %in% filtered_modules) %>%
        ungroup()


    tf_module_adj <- get_module_transitions(
        trajectory_object = trajectory_object,
        closest_module = closest_nodes_to_module[tf_modules],
        similarity_values = module_medians[tf_modules, , drop = FALSE]
    )

    tf_gene_graph <- get_tf_gene_network(
        tf_stats = tfs,
        module_adjacency = tf_module_adj,
        top_n_factors = n_top_tfs,
        hub_genes = gene_hub_stats
    )

    tf_graph_list[[dts_name]] <- plot_module_tfs_ggraph(
        module_tf_g = tf_gene_graph,
        exclude_non_hub_genes = TRUE
    ) + ggtitle(paste(dts_name, " - top ", n_top_tfs, " TFs per module")) +
        theme(
            plot.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            legend.position = "bottom",
            legend.direction = "horizontal",
            # stack one over the other
            legend.box = "vertical"
        )


    gene_hub_stats <- rhdf5::h5read(sum_path, paste0(selected_module, "/gene_hub_stats"))
    gene_hub_stats$gene <- rhdf5::h5read(sum_path, "genes") %>% as.character()
    gene_hub_stats <- gene_hub_stats %>%
        group_by(module) %>%
        slice_max(combined_score, n = 15) %>%
        filter(module %in% filtered_modules) %>%
        ungroup()
        tfs <- add_tf_hub_stats(tfs, gene_hub_stats)

    tf_bubbleplot_list[[dts_name]] <- plot_module_tfs_bubbleplot(
        tf_stats = tfs,
        point_size_range = c(1, 7)
    ) +
    coord_flip() +
    ggtitle(paste(dts_name, " - top 10 TFs per module")) +
    theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)
    ) +
    guides(size = guide_legend(nrow = 1))


}

row_bubble <- wrap_plots(tf_bubbleplot_list, ncol = 1)
row_tf <- wrap_plots(tf_graph_list, nrow = 1, widths = c(2, 1, 2))
final_fig <- wrap_plots(row_bubble, row_tf, ncol = 1, heights = c(3, 1))

pdf(file.path(output_dir, "sup_starlng_hub_genes_tf.pdf"), width = 14, height = 20)
print(final_fig)
dev.off()
