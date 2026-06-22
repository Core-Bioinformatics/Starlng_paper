# single-cell analysis package
library(Seurat)
library(ClustAssess)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)


setwd("comparison_against_ground_truth")

dts_path <- "../data/immuneCellsSCTransformed.rds"
dts_prefix <- "masld_immune"

library(Seurat)
library(dplyr)
library(Starlng)

library(WGCNA)
library(hdWGCNA)


enableWGCNAThreads(4)

if (basename(getwd()) != "comparison_against_ground_truth") {
    setwd("comparison_against_ground_truth")
}
n_tfs <- c(5, 10, 15, 20)

cluster_output <- list()
for (n_tf in n_tfs) {
    cluster_output[[as.character(n_tf)]] <- list()
    for (count_type in c("raw", "noisy")) {
        cluster_output[[as.character(n_tf)]][[count_type]] <- list()
        for (generation_type in c("", "dynamic_ds12_")) {
            print(paste0("Processing ", n_tf, " TFs", generation_type, " ", count_type, " counts"))
            if (generation_type %in% names(cluster_output[[as.character(n_tf)]][[count_type]])) {
                next
            }
            n_soft_clusters <- ceiling(n_tf * 0.6)


            input_file <- file.path("..", "data", paste0("sergio_", generation_type, count_type, "_", n_tf, "_tfs_seurat.qs2"))

            so <- qs2::qs_read(input_file)
            expr_matrix <- GetAssayData(so, assay = "RNA", layer = "data")
            colnames(expr_matrix) <- as.character(colnames(expr_matrix))
            rownames(expr_matrix) <- as.character(rownames(expr_matrix))
            cell_mask <- apply(expr_matrix, 2, function(x) sum(as.numeric(x)) > 0)
            so <- so[, cell_mask]

            modules <- list()
            so <- SetupForWGCNA(
                so,
                gene_select = "fraction",
                fraction = 0.05,
                wgcna_name = dts_prefix
            )

            selected_meta <- "gt_clusters"
            so <- MetacellsByGroups(
                seurat_obj = so,
                group.by = c(selected_meta),
                reduction = "pca",
                layer = "data",
                slot = "data",
                k = 25,
                max_shared = 10,
                min_cells = 40,
                ident.group = selected_meta,
                wgcna_name = dts_prefix

            )
    
            so <- NormalizeMetacells(so)
            so <- SetDatExpr(
                seurat_obj = so,
                group.by = selected_meta,
                group_name = so@misc[[dts_prefix]]$wgcna_params$metacell_stats$gt_clusters,
                assay = "RNA",
                slot = "data",
                layer = "data"
            )

            so <- TestSoftPowers(
                seurat_obj = so,
                networkType = "signed"
            )

            # plot_list <- PlotSoftPowers(so)
            # wrap_plots(plot_list, ncol = 2)

            power_table <- GetPowerTable(so)
            head(power_table)

            is_power_inf <- TRUE
            try({
                so <- ConstructNetwork(
                    seurat_obj = so,
                    tom_name = dts_prefix,
                    overwrite_tom = TRUE,
                    wgcna_name = dts_prefix
                )
                is_power_inf <- FALSE
            }, silent = TRUE)
            if (is_power_inf) {
                if (!is.null(power_table) && nrow(power_table) > 0) {
                    ma <- power_table$Power[which.min(power_table$median.k.)]
                } else {
                    ma <- 3
                }
                so <- ConstructNetwork(
                    seurat_obj = so,
                    soft_power = ma,
                    tom_name = dts_prefix,
                    overwrite_tom = TRUE,
                    wgcna_name = dts_prefix
                )
            }

            PlotDendrogram(so)
            cluster_output[[as.character(n_tf)]][[count_type]][[generation_type]] <- GetModules(so)
        }
    }
}
qs2::qs_save(cluster_output, "hdwgcna_cluster_output.qs2")
