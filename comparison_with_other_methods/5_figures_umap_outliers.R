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
gradient <- c("grey85", "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84",
            "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000")

input_dir <- "comparison_files"
output_dir <- "panels"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

module_stats <- qs2::qs_read(file.path(input_dir, "3_module_stats_summary.qs2"), nthreads = 30)
avg_summs <- qs2::qs_read(file.path(input_dir, "3_avg_summaries.qs2"), nthreads = 30)


for (dts_name in names(module_stats)) {
    if (dts_name %in% c("cao_Liver", "cao_Lung")) {
        next
    }
    umap_plot_list <- list()

    so <- qs2::qs_read(file.path("..", "data", paste0(dts_name, "_filtered_normalized.qs2")), nthreads = 30)
    expr_matrix <- GetAssayData(so, assay = "RNA", layer = "data")
    umap_emb <- data.frame(Embeddings(so, reduction = "umap"))
    pseudotime_vals <- read.csv(file.path("..", "data", paste0(dts_name, "_recommended_pseudotime.csv")))
    pseudotime_vals <- pseudotime_vals$recommended_pseudotime
    for (method_name in names(module_stats[[dts_name]])) {
        print(paste0("Processing ", dts_name, " - ", method_name))
        if (method_name %in% umap_plot_list[[dts_name]]) {
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
            filter(is_outlier != "no") %>%
            arrange(median_pseudotime) %>%
            pull(module)
        method_avg <- avg_summs[[dts_name]][[method_name]][selected_modules]
        method_clust_df <- read.csv(input_file)
        method_clust_df <- method_clust_df %>% filter(module %in% selected_modules)

        selected_outlier <- method_stats %>%
            filter(is_outlier == "yes") %>%
            arrange(desc(iqr_pseudotime)) %>%
            pull(module) %>%
            head(1)
        if (length(selected_outlier) == 0) {
            plot_outlier <- ggplot() + theme_void() 
            selected_outlier <- ""
        } else {
            umap_emb$colour <- avg_summs[[dts_name]][[method_name]][[selected_outlier]]
            umap_emb$colour[umap_emb$colour < 0.5] <- 0

            plot_outlier <- ggplot() +
                geom_point(data = umap_emb %>% filter(colour == 0), aes(x = umap_1, y = umap_2), colour = "grey85", size = 0.5) +
                geom_point(data = umap_emb %>% filter(colour > 0), aes(x = umap_1, y = umap_2, colour = colour), size = 1) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                scale_colour_gradientn(colours = gradient, limits = c(0, 1)) +
                theme_void() +
                theme(legend.position = "right", plot.title = element_text(size = 20), legend.title = element_text(size = 18), legend.text = element_text(size = 16)) +
                ggtitle(paste0(method_name, ": Outlier module ", selected_outlier))
        }

        selected_redundant <- method_stats %>%
            filter(is_outlier %in% c("redundant", "yes"), module != selected_outlier) %>%
            arrange(is_outlier, desc(iqr_pseudotime)) %>%
            head(1)

        if (nrow(selected_redundant) == 0) {
            plot_redundant <- ggplot() + theme_void()
        } else {
            module_type <- selected_redundant$is_outlier
            selected_redundant <- selected_redundant$module
            umap_emb$colour <- avg_summs[[dts_name]][[method_name]][[selected_redundant]]
            umap_emb$colour[umap_emb$colour < 0.5] <- 0

            plot_redundant <- ggplot() +
                geom_point(data = umap_emb %>% filter(colour == 0), aes(x = umap_1, y = umap_2), colour = "grey85", size = 0.5) +
                geom_point(data = umap_emb %>% filter(colour > 0) %>% arrange(colour), aes(x = umap_1, y = umap_2, colour = colour), size = 1) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                scale_colour_gradientn(colours = gradient, limits = c(0, 1)) +
                theme_void() +
                theme(legend.position = "right", plot.title = element_text(size = 20), legend.title = element_text(size = 18), legend.text = element_text(size = 16)) +
                ggtitle(paste0(method_name, ": ", module_type, " module ", selected_redundant))

        }
        
        umap_plot_list[[method_name]] <- patchwork::wrap_plots(plot_outlier, plot_redundant, ncol = 2) 
    }

    patchwork_plot <- wrap_plots(umap_plot_list, ncol = 1) + plot_annotation(title = paste0(dts_name, ": UMAP plots of outlier and redundant modules"), theme = theme(plot.title = element_text(size = 30)))
    pdf(file.path(output_dir, paste0("sup_umap_outliers_", dts_name, ".pdf")), width = 14, height = 20)
    print(patchwork_plot)
    dev.off()
}