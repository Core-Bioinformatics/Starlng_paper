library(Seurat)
library(dplyr)
library(Starlng)
library(ggplot2)
library(ggpattern)

if (basename(getwd()) != "comparison_against_ground_truth") {
    setwd("comparison_against_ground_truth")
}
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
            "ground_truth_ecs",
            "ground_truth_discarded_genes"
        )))
}



n_tfs <- c(5, 10, 15, 20)
hdwgcna_result <- qs2::qs_read("hdwgcna_cluster_output.qs2")
hotspot_result <- jsonlite::fromJSON("hotspot_cluster_output.json")
scenic_result <- qs2::qs_read("scenic_cluster_output.qs2")
paga_cellrank_result <- jsonlite::fromJSON("paga_cellrank_cluster_output.json")

ecc_per_method_dataset <- NULL
stats_per_method_dataset <- NULL
plot_list <- list()
index <- 1

method_colours <- qualpalr::qualpal(7)$hex
names(method_colours) <- c("Starlng", "Hotspot", "hdWGCNA", "SCENIC", "PAGA+\nCellRank", "Starlng\n(low Moran's I)", "Starlng\n(high Moran's I)")
for (generation_type in c("", "dynamic_ds12_")) {
    # plot_list[[generation_type]] <- list()
    for (count_type in c("raw", "noisy")) {
        for (n_tf in n_tfs) {
            print(paste0("Processing ", n_tf, " TFs", generation_type, " ", count_type, " counts"))
            n_soft_clusters <- ceiling(n_tf * 0.6)
            
            gt <- readRDS(paste0("../data/sergio_ground_truth_", n_tf, "_tfs.rds"))
            gt_index <- gt$index
            gt_index_reverse <- setNames(names(gt_index), as.character(gt_index))
            gt_mapping <- gt$tf_mapping
            gt_mapping$tf_name <- gt_index_reverse[as.character(gt_mapping$tf_name)]
            gt_mapping$gene_name <- gt_index_reverse[as.character(gt_mapping$gene_name)]
            gt_mapping <- rbind(gt_mapping, data.frame(
                tf_name = unique(gt_mapping$tf_name),
                gene_name = unique(gt_mapping$tf_name),
                weight = 0,
                nregs = 0,
                hill = 0
            ))
            gt_cluster <- setNames(as.integer(factor(gt_mapping$tf_name)), gt_mapping$gene_name)


            # Starlng

            app_folder <- paste0("starlng/sergio_", generation_type, count_type, "_", n_tf, "_tfs_starlng_app")
            if (!file.exists(file.path(app_folder, "objects", "module_summaries.h5"))) {
                next
            }

            starlng_module_summs <- rhdf5::h5read(file.path(app_folder, "objects", "module_summaries.h5"), "/")
            used_genes <- as.character(starlng_module_summs$genes)
            available_n_modules <- as.character(starlng_module_summs$all_modules)
            starlng_modules <- do.call(cbind, lapply(available_n_modules, function(n_module) {
                starlng_module <- as.integer(rhdf5::h5read(file.path(app_folder, "objects", "module_summaries.h5"), paste0("/", n_module, "/clustering")))
                names(starlng_module) <- starlng_module_summs$genes
                return(starlng_module[used_genes])
            }))
            colnames(starlng_modules) <-  available_n_modules
            diff_genes <- setdiff(names(gt_cluster), used_genes)
            if (length(diff_genes) > 0) {
                starlng_modules <- rbind(starlng_modules, matrix(-1, nrow = length(diff_genes), ncol = ncol(starlng_modules), dimnames = list(diff_genes, colnames(starlng_modules))))
            }
            starlng_modules <- starlng_modules[names(gt_cluster), , drop = FALSE]

            ecs_scores_all <- apply(starlng_modules, 2, function(module) {
                ClustAssess::element_sim(module, gt_cluster)
            })
            which_cluster_max <- names(ecs_scores_all)[which.max(ecs_scores_all)]
            ecs_scores_all <- max(ecs_scores_all)
            starlng_ecs <- ClustAssess::element_sim_elscore(starlng_modules[, which_cluster_max], gt_cluster)
            starlng_modules <- starlng_modules[, which_cluster_max]
            starlng_n_clusters <- length(unique(starlng_modules))
            starlng_n_removed <- sum(starlng_modules == -1)

            temp_ecs_df <- data.frame(
                method = "starlng",
                n_tf = n_tf,
                count_type = count_type,
                generation_type = generation_type,
                ecs = starlng_ecs
            )
            temp_stats_df <- data.frame(
                method = "starlng",
                n_tf = n_tf,
                count_type = count_type,
                generation_type = generation_type,
                n_clusters = starlng_n_clusters,
                n_cells = starlng_n_removed,
                cells_type = "removed"
            )
            temp_stats_df <- rbind(temp_stats_df,
                data.frame(
                    method = "starlng",
                    n_tf = n_tf,
                    count_type = count_type,
                    generation_type = generation_type,
                    n_clusters = starlng_n_clusters,
                    n_cells = 0,
                    cells_type = "outliers"
                )
            )


            # hdWGCNA
            hdwgcna_modules <- hdwgcna_result[[as.character(n_tf)]][[count_type]]
            if (generation_type == "") {
                hdwgcna_modules <- hdwgcna_modules[[1]]
            } else {
                hdwgcna_modules <- hdwgcna_modules[[generation_type]]
            }
            hdwgcna_modules <- setNames(as.character(hdwgcna_modules$module), hdwgcna_modules$gene_name)
            hdwgcna_n_outliers <- sum(hdwgcna_modules == "grey")

            diff_genes <- setdiff(names(gt_cluster), names(hdwgcna_modules))
            if (length(diff_genes) > 0) {
                hdwgcna_modules <- c(hdwgcna_modules, setNames(rep("grey", length(diff_genes)), diff_genes))
            }
            hdwgcna_n_removed <- length(diff_genes)
            hdwgcna_modules <- hdwgcna_modules[names(gt_cluster)]

            hdwgcna_ecs <- ClustAssess::element_sim_elscore(as.integer(factor(hdwgcna_modules)), gt_cluster)
            hdwgcna_n_clusters <- length(unique(hdwgcna_modules))

            temp_ecs_df <- rbind(temp_ecs_df, data.frame(
                method = "hdwgcna",
                n_tf = n_tf,
                count_type = count_type,
                generation_type = generation_type,
                ecs = hdwgcna_ecs
            ))
            temp_stats_df <- rbind(temp_stats_df, data.frame(
                method = "hdwgcna",
                n_tf = n_tf,
                count_type = count_type,
                generation_type = generation_type,
                n_clusters = hdwgcna_n_clusters,
                n_cells = hdwgcna_n_removed,
                cells_type = "removed"
            ))
            temp_stats_df <- rbind(temp_stats_df, data.frame(
                method = "hdwgcna",
                n_tf = n_tf,
                count_type = count_type,
                generation_type = generation_type,
                n_clusters = hdwgcna_n_clusters,
                n_cells = hdwgcna_n_outliers,
                cells_type = "outliers"
            ))

            # hotspot
            hotspot_modules <- hotspot_result[[as.character(n_tf)]][[count_type]]
            if (generation_type == "") {
                hotspot_modules <- hotspot_modules[[1]]
            } else {
                hotspot_modules <- hotspot_modules[[generation_type]]
            }
            hotspot_modules <- setNames(hotspot_modules$clusters, hotspot_modules$genes)
            diff_genes <- setdiff(names(gt_cluster), names(hotspot_modules))
            hotspot_n_outliers <- sum(hotspot_modules == -1)
            if (length(diff_genes) > 0) {
                hotspot_modules <- c(hotspot_modules, setNames(rep(-1, length(diff_genes)), diff_genes))
            }
            hotspot_n_removed <- length(diff_genes)
            hotspot_modules <- hotspot_modules[names(gt_cluster)]
            hotspot_ecs <- ClustAssess::element_sim_elscore(hotspot_modules, gt_cluster)
            hotspot_n_clusters <- length(unique(hotspot_modules))

            temp_ecs_df <- rbind(temp_ecs_df, data.frame(
                method = "hotspot",
                n_tf = n_tf,
                count_type = count_type,
                generation_type = generation_type,
                ecs = hotspot_ecs
            ))
            temp_stats_df <- rbind(temp_stats_df, data.frame(
                method = "hotspot",
                n_tf = n_tf,
                count_type = count_type,
                generation_type = generation_type,
                n_clusters = hotspot_n_clusters,
                n_cells = hotspot_n_removed,
                cells_type = "removed"
            ))
            temp_stats_df <- rbind(temp_stats_df, data.frame(
                method = "hotspot",
                n_tf = n_tf,
                count_type = count_type,
                generation_type = generation_type,
                n_clusters = hotspot_n_clusters,
                n_cells = hotspot_n_outliers,
                cells_type = "outliers"
            ))

            # scenic
            scenic_modules <- scenic_result[[as.character(n_tf)]][[count_type]]
            if (generation_type == "") {
                scenic_modules <- scenic_modules[[1]]
            } else {
                scenic_modules <- scenic_modules[[generation_type]]
            }
            scenic_modules <- setNames(scenic_modules$tf_name, scenic_modules$gene_name)
            added_genes <- intersect(names(gt_cluster), names(scenic_modules))
            scenic_modules <- as.integer(factor(scenic_modules[added_genes]))
            names(scenic_modules) <- added_genes
            diff_genes <- setdiff(names(gt_cluster), names(scenic_modules))
            if (length(diff_genes) > 0) {
                scenic_modules <- c(scenic_modules, setNames(rep(-1, length(diff_genes)), diff_genes))
            }
            scenic_modules <- scenic_modules[names(gt_cluster)]
            scenic_ecs <- ClustAssess::element_sim_elscore(scenic_modules, gt_cluster)
            scenic_n_clusters <- length(unique(scenic_modules))
            scenic_n_removed <- sum(scenic_modules == -1)

            temp_ecs_df <- rbind(temp_ecs_df, data.frame(
                method = "scenic",
                n_tf = n_tf,
                count_type = count_type,
                generation_type = generation_type,
                ecs = scenic_ecs
            ))
            temp_stats_df <- rbind(temp_stats_df, data.frame(
                method = "scenic",
                n_tf = n_tf,
                count_type = count_type,
                generation_type = generation_type,
                n_clusters = scenic_n_clusters,
                n_cells = scenic_n_removed,
                cells_type = "removed"
            ))
            temp_stats_df <- rbind(temp_stats_df, data.frame(
                method = "scenic",
                n_tf = n_tf,
                count_type = count_type,
                generation_type = generation_type,
                n_clusters = scenic_n_clusters,
                n_cells = 0,
                cells_type = "outliers"
            ))

            # paga + cellrank
            paga_modules <- paga_cellrank_result[[as.character(n_tf)]][[count_type]]
            if (generation_type == "") {
                paga_modules <- paga_modules[[1]]
            } else {
                paga_modules <- paga_modules[[generation_type]]
            }
            paga_modules <- setNames(as.integer(paga_modules$clusters), paga_modules$genes)
            diff_genes <- setdiff(names(gt_cluster), names(paga_modules))
            if (length(diff_genes) > 0) {
                paga_modules <- c(paga_modules, setNames(rep(-1, length(diff_genes)), diff_genes))
            }
            paga_modules <- paga_modules[names(gt_cluster)]
            paga_ecs <- ClustAssess::element_sim_elscore(paga_modules, gt_cluster)
            paga_n_clusters <- length(unique(paga_modules))
            paga_n_outliers <- sum(paga_modules == -1)
            temp_ecs_df <- rbind(temp_ecs_df, data.frame(
                method = "paga_cellrank",
                n_tf = n_tf,
                count_type = count_type,
                generation_type = generation_type,
                ecs = paga_ecs
            ))
            temp_stats_df <- rbind(temp_stats_df, data.frame(
                method = "paga_cellrank",
                n_tf = n_tf,
                count_type = count_type,
                generation_type = generation_type,
                n_clusters = paga_n_clusters,
                n_cells = 0,
                cells_type = "outliers"
            ))
            temp_stats_df <- rbind(temp_stats_df, data.frame(
                method = "paga_cellrank",
                n_tf = n_tf,
                count_type = count_type,
                generation_type = generation_type,
                n_clusters = paga_n_clusters,
                n_cells = paga_n_outliers,
                cells_type = "removed"
            ))

            temp_stats_df$n_total_genes <- length(gt_cluster)
            temp_stats_df$percent_cells <- temp_stats_df$n_cells / temp_stats_df$n_total_genes * 100
            temp_stats_df$percent_cells[temp_stats_df$cells_type == "outliers"] <- temp_stats_df$percent_cells[temp_stats_df$cells_type == "removed"] + temp_stats_df$percent_cells[temp_stats_df$cells_type == "outliers"]

            if (is.null(ecc_per_method_dataset)) {
                ecc_per_method_dataset <- temp_ecs_df
                stats_per_method_dataset <- temp_stats_df
            } else {
                ecc_per_method_dataset <- rbind(ecc_per_method_dataset, temp_ecs_df)
                stats_per_method_dataset <- rbind(stats_per_method_dataset, temp_stats_df)
            }

            # plot
            so <- qs2::qs_read(paste0("../data/sergio_", generation_type, count_type, "_", n_tf, "_tfs_seurat.qs2"))
            cluster_colours <- setNames(qualpalr::qualpal(length(unique(so$gt_clusters)))$hex, as.character(unique(so$gt_clusters)))
            plot_list[[index]] <- DimPlot(so, reduction = "umap", group.by = "gt_clusters", cols = cluster_colours, label = TRUE) +
                labs(title = paste0("GT - ", n_tf, " TF\n ", ifelse(generation_type == "", "default", "trajectory-based"), " ", count_type, " counts")) +
                scale_fill_manual(values = cluster_colours) +
                theme(legend.position = "bottom")


            # starlng vs hotspot modules heatmap
            unique_starlng_clusters <- stringr::str_sort(unique(starlng_modules), numeric = TRUE)
            unique_hotspot_clusters <- stringr::str_sort(unique(hotspot_modules), numeric = TRUE)
            modules_jsi <- matrix(NA, nrow = length(unique_starlng_clusters), ncol = length(unique_hotspot_clusters))
            rownames(modules_jsi) <- unique_starlng_clusters
            colnames(modules_jsi) <- unique_hotspot_clusters
            for (i in unique_starlng_clusters) {
                for (j in unique_hotspot_clusters) {
                    common_genes <- intersect(names(starlng_modules)[starlng_modules == i], names(hotspot_modules)[hotspot_modules == j])
                    union_genes <- union(names(starlng_modules)[starlng_modules == i], names(hotspot_modules)[hotspot_modules == j])
                    jsi <- length(common_genes) / length(union_genes)
                    modules_jsi[as.character(i), as.character(j)] <- jsi
                }
            }
            jsi_df <- reshape2::melt(modules_jsi, varnames = c("starlng_cluster", "hotspot_cluster"), value.name = "jsi")
            jsi_df$starlng_cluster <- factor(jsi_df$starlng_cluster, levels = unique_starlng_clusters)
            jsi_df$hotspot_cluster <- factor(jsi_df$hotspot_cluster, levels = unique_hotspot_clusters)
            pdf(paste0("jsi_heatmap_starlng_hotspot_", generation_type, count_type, "_", n_tf, "_tfs.pdf"), width = 5 + length(unique_hotspot_clusters) * 0.5, height = 5 + length(unique_starlng_clusters) * 0.5)
            print(ggplot(jsi_df, aes(x = hotspot_cluster, y = starlng_cluster, fill = jsi)) +
                geom_tile() +
                scale_fill_viridis_c(limits = c(0, 1), option = "viridis") +
                geom_text(aes(label = round(jsi, 2)), color = "white", size = 5) +
                labs(title = paste0("Jaccard Similarity Index between\nStarlng and Hotspot modules\n", n_tf, " TFs ", ifelse(generation_type == "", "default", "trajectory-based"), " ", count_type, " counts"), x = "Hotspot Cluster", y = "Starlng Cluster", fill = "JSI") +
                theme_minimal() +
                theme(
                    axis.text = element_text(size = 12),
                    axis.title = element_text(size = 14),
                    plot.title = element_text(size = 16, hjust = 0.5)
                )
            )
            dev.off()
        
            index <- index + 1
        }
    }
}


ecc_per_method_dataset$generation_type <- ifelse(ecc_per_method_dataset$generation_type == "", "default", "trajectory-based")
ecc_per_method_dataset$generation_type <- factor(ecc_per_method_dataset$generation_type, levels = c("default", "trajectory-based"))
ecc_per_method_dataset$count_type <- factor(ecc_per_method_dataset$count_type, levels = c("raw", "noisy"))
ecc_per_method_dataset$n_tf <- factor(ecc_per_method_dataset$n_tf, levels = c(5, 10, 15, 20))

ecc_per_method_dataset$dataset_descr <- paste0(ecc_per_method_dataset$generation_type, " ", ecc_per_method_dataset$count_type, " ", ecc_per_method_dataset$n_tf, " TFs")
# order the dataset descr to be split by generation type, then count type, then n_tf
ecc_per_method_dataset$dataset_descr <- factor(ecc_per_method_dataset$dataset_descr, levels = paste0(rep(c("default", "trajectory-based"), each = 8), " ", rep(rep(c("raw", "noisy"), each = 4), 2), " ", rep(c(5, 10, 15, 20), 4), " TFs"))
ecc_per_method_dataset$method <- recode(ecc_per_method_dataset$method, "starlng" = "Starlng", "hotspot" = "Hotspot", "hdwgcna" = "hdWGCNA", "scenic" = "SCENIC", "paga_cellrank" = "PAGA+\nCellRank")


pdf("ecs_comparison.pdf", width = 12, height = 6)
ggplot(ecc_per_method_dataset, aes(x = dataset_descr, y = ecs, fill = method)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = method_colours) +
    geom_vline(xintercept = seq(4.5, 12.5, by = 4), linetype = "dashed", color = "grey") +
    labs(x = "Dataset", y = "ECS", fill = "Method", title = "Comparison of clustering methods against ground truth")
dev.off()


stats_per_method_dataset$generation_type <- ifelse(stats_per_method_dataset$generation_type == "", "default", "trajectory-based")
stats_per_method_dataset$generation_type <- factor(stats_per_method_dataset$generation_type, levels = c("default", "trajectory-based"))
stats_per_method_dataset$count_type <- factor(stats_per_method_dataset$count_type, levels = c("raw", "noisy"))
stats_per_method_dataset$n_tf <- factor(stats_per_method_dataset$n_tf, levels = c(5, 10, 15, 20))
stats_per_method_dataset$dataset_descr <- paste0(stats_per_method_dataset$generation_type, " ", stats_per_method_dataset$count_type, " ", stats_per_method_dataset$n_tf, " TFs")
stats_per_method_dataset$dataset_descr <- factor(stats_per_method_dataset$dataset_descr, levels = paste0(rep(c("default", "trajectory-based"), each = 8), " ", rep(rep(c("raw", "noisy"), each = 4), 2), " ", rep(c(5, 10, 15, 20), 4), " TFs"))
stats_per_method_dataset$method <- recode(stats_per_method_dataset$method, "starlng" = "Starlng", "hotspot" = "Hotspot", "hdwgcna" = "hdWGCNA", "scenic" = "SCENIC", "paga_cellrank" = "PAGA+\nCellRank")

pdf("percentage_outliers_comparison.pdf", width = 12, height = 6)
# plot percentage outliers and percentage removed stacekd, and the outlier should have the barplot striped
ggplot() +
    geom_bar_pattern(data = stats_per_method_dataset %>% filter(cells_type == "outliers"), aes(x = dataset_descr, y = percent_cells, fill = method, pattern = cells_type), stat = "identity", position = "dodge", pattern_color = "black", pattern_spacing = 0.04, pattern_density = 0.05, color = "black", pattern = "stripe") +
    geom_bar(data = stats_per_method_dataset %>% filter(cells_type == "removed"), aes(x = dataset_descr, y = percent_cells, fill = method), stat = "identity", position = "dodge", color = "black") +
    scale_fill_manual(values = method_colours) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_vline(xintercept = seq(4.5, 12.5, by = 4), linetype = "dashed", color = "grey", linewidth = 1) +
    labs(x = "Dataset", y = "Percentage of Cells", fill = "Method", pattern_fill = "Method", title = "Percentage of outlier and removed cells identified by each method across datasets") 
dev.off()

# ggplot(stats_per_method_dataset, aes(x = dataset_descr, y = percent_outliers, fill = method)) +
#     geom_bar(stat = "identity", position = "dodge") +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     scale_fill_manual(values = method_colours) +
#     geom_vline(xintercept = seq(4.5, 12.5, by = 4), linetype = "dashed", color = "grey") +
#     labs(x = "Dataset", y = "Percentage Outlier Cells", fill = "Method", title = "Percentage of outlier cells identified by each method across datasets")

n_genes <- stats_per_method_dataset %>%
    group_by(dataset_descr) %>%
    summarise(n_total_genes = n_total_genes[1]) %>%
    ungroup() %>%
    pull(n_total_genes) %>%
    sum()
method_ranking <- stats_per_method_dataset %>%
    group_by(method) %>%
    summarise(perc_removed = sum(n_cells) / n_genes * 100) %>%
    ungroup()
colnames(method_ranking) <- c("method_name", "value")
method_ranking$rank <- rank(method_ranking$value, ties.method = "min")
method_ranking$method_name <- gsub("\n", "", method_ranking$method_name)
method_ranking$dataset_name <- "sergio"
method_ranking$comparison_type <- "ground_truth_discarded_genes"

comparison_table <- rbind(comparison_table, method_ranking[, colnames(comparison_table)])


method_ranking <- ecc_per_method_dataset %>%
    group_by(method) %>%
    summarise(ecs = median(ecs)) %>%
    ungroup()
colnames(method_ranking) <- c("method_name", "value")
method_ranking$rank <- rank(-method_ranking$value, ties.method = "min")
method_ranking$method_name <- gsub("\n", "", method_ranking$method_name)
method_ranking$dataset_name <- "sergio"
method_ranking$comparison_type <- "ground_truth_ecs"

comparison_table <- rbind(comparison_table, method_ranking[, colnames(comparison_table)])

write.csv(comparison_table, comparison_table_file, row.names = FALSE)
