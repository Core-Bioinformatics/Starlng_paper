library(Seurat)
library(dplyr)
library(ggplot2)
library(Starlng)
library(qualpalr)

psd_violin_plot <- function(psd_vec, module_cell_list, selected_modules = NULL) {
    if (is.null(selected_modules)) {
        selected_modules <- names(module_cell_list)
    }
    df <- NULL
    for (i in selected_modules) {
        temp_df <- data.frame(pseudotime = psd_vec[module_cell_list[[i]]])
        temp_df$module <- i
        if (is.null(df)) {
            df <- temp_df
        } else {
            df <- rbind(df, temp_df)
        }
    }
    cols <- qualpal(length(module_cell_list))$hex
    names(cols) <- names(module_cell_list)
    df$module <- factor(df$module, selected_modules)
    mds <- df %>% group_by(module) %>% summarise(med = median(pseudotime)) %>% arrange(med)
    df$module <- factor(df$module, levels = mds$module)

    return(
        ggplot(df, aes(x = module, y = pseudotime, fill = module)) +
            geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75), adjust = 2, size = 0.2) +
            theme_classic() +
            theme(
                legend.position = "none",
                axis.text = element_text(size = 18),
                axis.title = element_text(size = 18),
                axis.text.x = element_text(angle = 45, hjust = 1)
            ) +
            labs(x = "Module", y = "Pseudotime") +
            scale_fill_manual(values = cols[levels(df$module)])
    )
}

if (basename(getwd()) != "starlng_run") {
    setwd("starlng_run")
}

dts_combinations <- list(
    "masld_immune" = c("low_moran", "high_moran"),
    "cao_Liver_subset" = c(""),
    "cao_Pancreas" = c("")
)

enrichment_all <- qs2::qs_read(file.path("..", "comparison_with_other_methods", "comparison_files", "4_enrichment_analysis.qs2"), nthreads = 30)

for (dts_name in names(dts_combinations)) {
    pseudotime_vals <- read.csv(file.path("..", "data", paste0(dts_name, "_recommended_pseudotime.csv")))
    pseudotime_vals <- pseudotime_vals$recommended_pseudotime
    for (moran_type in dts_combinations[[dts_name]]) {
        if (moran_type == "") {
            prefix <- ""
        } else {
            prefix <- paste0(moran_type, "_")
        }
        app_dir <- file.path(paste0(dts_name, "_", prefix, "starlng_app"), "objects")
        output_dir <- file.path("plots_panels", paste0(dts_name, "_", prefix))
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        summ_path <- file.path(app_dir, "module_summaries.h5")

        trajectory_object <- qs2::qs_read(file.path(app_dir, "trajectory_object.qs2"), nthreads = 30)
        all_modules <- rhdf5::h5read(summ_path, "all_modules")
        n_modules <- as.character(max(as.numeric(all_modules)))
        module_names <- rhdf5::h5read(summ_path, paste0(n_modules, "/modules")) %>% as.character()
        module_list <- split(rhdf5::h5read(summ_path, "genes") %>% as.character(), rhdf5::h5read(summ_path, paste0(n_modules, "/clustering")))
        genes <- rhdf5::h5read(summ_path, "genes") %>% as.character()

        module_stats <- rhdf5::h5read(summ_path, paste0(n_modules, "/modules_stats_summary"))
        selected_modules <- module_stats %>%
            as.data.frame() %>%
            filter(is_outlier == "no") %>%
            arrange(median_pseudotime) %>%
            pull(module)
        avg_summ <- rhdf5::h5read(summ_path, paste0(n_modules, "/expression_summaries"))
        colnames(avg_summ) <- module_names
        scaled_avg_summ <- apply(avg_summ, 2, function(x) {
            x <- (x - min(x)) / (max(x) - min(x))
            x > 0.5
        })
        scaled_avg_summ <- split(scaled_avg_summ, col(scaled_avg_summ))

        # PSEUDOTIME
        pdf(file.path(output_dir, paste0("pseudotime_violin_plot.pdf")), width = 10, height = 6)
        print(psd_violin_plot(psd_vec = pseudotime_vals, module_cell_list = scaled_avg_summ[selected_modules]) +
            ggtitle("Per-module Pseudotime Distribution") +
            theme(
                plot.title = element_text(size = 25, hjust = 0.5)
            )
        )
        dev.off() 

        # UMAP Hubs
        umap_df <- qs2::qs_read(file.path(app_dir, "metadata.qs2"), nthreads = 30)
        umap_df <- umap_df[, c(ncol(umap_df) - 1, ncol(umap_df))]
        emb <- qs2::qs_read(file.path(app_dir, "gene_embedding.qs2"), nthreads = 30)
        gene_umap <- emb$umap
        gene_adj_matrix <- emb$adj_matrix
        module_medians <- get_module_centroid(
            module_expr = split(avg_summ, col(avg_summ)),
            cell_umap = umap_df,
            expression_threshold = 0.5,
            scale = TRUE
        )

        closest_nodes_to_module <- setNames(
            rhdf5::h5read(summ_path, paste0(n_modules, "/closest_nodes_to_module")),
            rhdf5::h5read(summ_path, paste0(n_modules, "/modules"))
        )
        gene_hub_stats <- rhdf5::h5read(summ_path, paste0(n_modules, "/gene_hub_stats"))
        gene_hub_stats$gene <- rhdf5::h5read(summ_path, "genes") %>% as.character()
        tfs <- rhdf5::h5read(summ_path, paste0(n_modules, "/tfs"))

        n_top_hubs <- 50 %/% length(selected_modules)

        gene_hub_stats <- gene_hub_stats %>%
            group_by(module) %>%
            slice_max(combined_score, n = n_top_hubs) %>%
            filter(module %in% selected_modules) %>%
            ungroup()

        module_adjacency <- get_module_transitions(
            trajectory_object,
            closest_nodes_to_module[selected_modules],
            start_node = NULL,
            similarity_values = module_medians[selected_modules, , drop = FALSE]
        )
        filtered_gene_adj <- get_filtered_gene_adjacency(
            gene_modules = module_list[selected_modules],
            module_adjacency = module_adjacency,
            gene_adjacency = gene_adj_matrix,
            hub_genes = gene_hub_stats,
            percentage_non_hub_nodes = 1,
            percentage_edges = ifelse(nrow(gene_adj_matrix) > 1500, 0.05, 0.25)
        )

        hub_nodes_df <- filtered_gene_adj$nodes_df %>%
            as.data.frame() %>%
            filter(is_hub)
        hub_label_df <- as.data.frame(gene_umap[hub_nodes_df$gene, , drop = FALSE])
        hub_label_df$gene <- hub_nodes_df$gene
        hub_label_df$module <- hub_nodes_df$module
        umap_cols <- colnames(hub_label_df)[1:2]

        # Draw hub labels with black text and module-coloured transparent boxes.
        hub_plot <- plot_gene_hub_umap(
            umap_df = gene_umap,
            filtered_gene_adj = filtered_gene_adj,
            legend_text_size = 10,
            point_size = 1,
            point_alpha = 0.02,
            edge_alpha = 0.2,
            hub_point_scale = 2.5,
            node_text_size = 0,
            label_box_alpha = 0,
            label_box_stroke = 0
        )
        hub_label_df$module <- factor(hub_label_df$module, levels = selected_modules)
        # if (nrow(hub_label_df) > 0) {
        #     hub_plot <- hub_plot +
        #         ggrepel::geom_label_repel(
        #             data = hub_label_df,
        #             aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]], label = gene, colour = module),
        #             size = 5.5,
        #             fill = ggplot2::alpha("white", 0),
        #             fontface = "bold",
        #             box.padding = 0.3,
        #             point.padding = 0.2,
        #             label.padding = unit(0.25, "lines"),
        #             max.overlaps = Inf,
        #             label.size = 1,
        #             min.segment.length = 0,
        #             seed = 42
        #         ) +
        #         ggrepel::geom_label_repel(
        #             data = hub_label_df,
        #             aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]], label = gene),
        #             size = 5.5,
        #             fill = ggplot2::alpha("white", 0),
        #             colour = "gray40",
        #             fontface = "bold",
        #             label.size = 0.01,
        #             max.overlaps = Inf,
        #             box.padding = 0.3,
        #             point.padding = 0.2,
        #             label.padding = unit(0.25, "lines"),
        #             min.segment.length = 0,
        #             segment.color = NA,
        #             seed = 42,
        #             show.legend = FALSE
        #         )
        # }
        if (nrow(hub_label_df) > 0) {
            hub_plot <- hub_plot +
                ggrepel::geom_text_repel(
                    data = hub_label_df,
                    aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]], label = gene, colour = module),
                    size = 5.5,
                    bg.color = "gray90",
                    bg.r = 0.1,
                    fontface = "bold",
                    box.padding = 0.3,
                    point.padding = 0.2,
                    max.overlaps = Inf,
                    min.segment.length = 0,
                    segment.color = NA,
                    show.legend = FALSE,
                    seed = 42
                )
        }

        pdf(file.path(output_dir, paste0("gene_hub_umap_plot.pdf")), width = 10, height = 8)
        print(hub_plot + ggtitle(paste("Top", n_top_hubs, "Hub Genes per Module")) +
            guides(linewidth = "none", colour = guide_legend(nrow = 2)) +
            theme(
                plot.title = element_text(size = 20, hjust = 0.5),
                legend.text = element_text(size = 10),
                legend.position = "bottom"
            )
        )
        dev.off()

        # STABILITY ASSESSMENT
        stab <- qs2::qs_read(file.path(app_dir, "full_stability_assessment.qs2"), nthreads = 30)
        if ("clusters_list" %in% names(stab)) {
            stab <- stab$clusters_list
        }
        df <- do.call(
            rbind,
            lapply(names(stab), function(x) {
                data.frame(
                    ecc = stab[[x]][[1]]$overall_ecc,
                    k = rep(x, length(stab[[x]][[1]]$overall_ecc))
                )
            })
        )
        df$k <- factor(df$k, levels = rev(names(stab)))

        pdf(file.path(output_dir, paste0("stability_ecc_n_neigh.pdf")), width = 5, height = 5)
        print(
        ggplot(df, aes(x = k, y = ecc)) +
            geom_boxplot() +
            theme_bw() +
            theme(
                axis.text = element_text(size = 18),
                axis.title = element_text(size = 18),
                plot.title = element_text(size = 20, hjust = 0.5)
            ) +
            labs(x = "number of nearest neighbours", y = "Overall ECC", title = "Stability Assessment of the\nnumber of nearest neighbours")
        )
        dev.off()

        best_config <- select_best_configuration(stab)
        best_nn <- best_config[1]
        best_qual <- best_config[2]

        df <- do.call(
            rbind,
            lapply(names(stab[[best_nn]][[best_qual]]$k), function(x) {
                data.frame(
                    mean_ecc = mean(stab[[best_nn]][[best_qual]]$k[[x]]$ecc),
                    iqr_ecc = IQR(stab[[best_nn]][[best_qual]]$k[[x]]$ecc),
                    freq = sum(sapply(stab[[best_nn]][[best_qual]]$k[[x]]$partitions, function(x) x$freq)),
                    k = x
                )
            })
        )
        rownames(df) <- df$k

        df$k <- factor(df$k, levels = names(stab[[best_nn]][[best_qual]]$k))
        cap_val <- 200
        df$freq[df$freq > cap_val] <- cap_val

        pdf(file.path(output_dir, paste0("stability_ecc_k.pdf")), width = 14)
        print(ggplot() +
            geom_point(data = df %>% filter(mean_ecc <= 0.9), aes(x = k, colour = iqr_ecc, y = freq), size = 1) +
            geom_point(data = df %>% filter(mean_ecc > 0.9), aes(x = k, colour = iqr_ecc, y = freq, size = mean_ecc)) +
            theme_bw() +
            scale_colour_viridis_c(direction = -1, name = "IQR\nECC") +
            scale_size_continuous(range = c(4, 7), name = "mean\nECC") +
            geom_hline(yintercept = 30, linetype = "dashed", colour = "red") +
            scale_x_discrete(
                limits = levels(df$k),
                labels = sapply(rownames(df), function(x) {
                    if (df[x, "mean_ecc"] < 0.9 || df[x, "freq"] < 30) {
                        return(paste0("<span style='color:red'>", x, "</span>"))
                    }
                    return(x)
                })
            ) +
            theme(
                axis.text = element_text(size = 20),
                axis.text.x = ggtext::element_markdown(),
                axis.title = element_text(size = 24),
                legend.text = element_text(size = 18),
                legend.title = element_text(size = 22),
                plot.title = element_text(size = 24, hjust = 0.5)
            ) +
            ggtitle("Stability Assessment of the number of modules")
        )
        dev.off()

        # ENRICHMENT
        suffix <- ifelse(moran_type == "", "", paste0("_", moran_type))
        current_enrichment <- enrichment_all[[dts_name]][[paste0("starlng", suffix)]]
        grouped_df <- do.call(rbind, lapply(selected_modules, function(mod_name) {
            x <- current_enrichment[[mod_name]]
            if (is.null(x)) {
                return(NULL)
            }
            x <- x$result
            if (is.null(x)) {
                return(NULL)
            }
            x <- x %>%
                filter(source == "GO:BP", term_size <= 500) 
            if (nrow(x) == 0) {
                return(NULL)
            }

            x$module <- mod_name
            return(x)
        }))

        gene_hub_stats <- rhdf5::h5read(summ_path, paste0(n_modules, "/gene_hub_stats"))
        gene_hub_stats$gene <- rhdf5::h5read(summ_path, "genes") %>% as.character()
        gene_hubs <- gene_hub_stats %>%
            group_by(module) %>%
            slice_max(combined_score, n = 15) %>%
            filter(module %in% selected_modules) %>%
            ungroup() %>%
            pull(gene)

        top_20_enrich_names <- grouped_df %>%
            group_by(module) %>%
            slice_min(p_value, n = 5) %>%
            ungroup() %>%
            arrange(p_value) %>%
            pull(term_name) %>%
            unique() %>%
            head(25)
        top_20_enrich_df <- grouped_df %>% filter(term_name %in% top_20_enrich_names) %>% arrange(p_value)
        top_20_enrich_df$intersection_hubs <- 0
        for (i in seq_len(nrow(top_20_enrich_df))) {
            term_genes <- strsplit(top_20_enrich_df$intersection[i], ",")[[1]]
            top_20_enrich_df$intersection_hubs[i] <- length(intersect(term_genes, gene_hubs)) / 15
        }
        top_20_enrich_df <- top_20_enrich_df[, c("term_name", "module", "intersection_hubs", "p_value")]
        top_20_enrich_df$p_value <- -log10(top_20_enrich_df$p_value)
        print(nrow(top_20_enrich_df))

        max_words_per_line <- 6
        top_20_enrich_df$term_name <- sapply(top_20_enrich_df$term_name, function(x) {
            words <- unlist(strsplit(x, " "))
            if (length(words) > max_words_per_line) {
                split_words <- split(words, ceiling(seq_along(words) / max_words_per_line))
                return(paste(sapply(split_words, paste, collapse = " "), collapse = "\n"))
            } else {
                return(x)
            }
        })

        unique_terms <- unique(top_20_enrich_df$term_name)
        term_line_counts <- sapply(unique_terms, function(lbl) {
            length(strsplit(lbl, "\n", fixed = TRUE)[[1]])
        })
        enrichment_dotplot_height <- max(10, min(35, 3 + 0.45 * sum(term_line_counts)))

        # sort the enriched terms based on the module ordering
        top_20_enrich_df$module <- factor(top_20_enrich_df$module, levels = selected_modules)
        top_20_enrich_df <- top_20_enrich_df %>% arrange(module, p_value)
        top_20_enrich_df$term_name <- factor(top_20_enrich_df$term_name, levels = unique(top_20_enrich_df$term_name))
        
        pdf(file.path(output_dir, paste0("enrichment_dotplot.pdf")), width = 12, height = enrichment_dotplot_height)
        print(
            ggplot(top_20_enrich_df, aes(x = module, y = term_name)) +
                geom_point(aes(size = intersection_hubs, colour = p_value)) +
                scale_colour_viridis_c(direction = 1, name = "-log10\n(p-value)") +
                scale_size_continuous(name = "%Hub Genes\nOverlap") +
                labs(x = "Module", y = "Enriched GO:BP Term") +
                scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
                theme_bw() +
                ggtitle(paste0("Enrichment Analysis\ntop 25 GO:BP terms")) +
                theme(
                    axis.text = element_text(size = 18),
                    axis.text.y = element_text(size = 16, lineheight = 0.95, margin = margin(r = 8)),
                    axis.title = element_text(size = 24),
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    plot.title = element_text(size = 30, hjust = 0.5),
                    legend.text = element_text(size = 18),
                    legend.title = element_text(size = 20)
                )
        )
        dev.off()




    }
}
