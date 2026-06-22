library(ggplot2)
library(dplyr)

if (basename(getwd()) != "moran_i_impact") {
    setwd("moran_i_impact")
}

apps_folder <- file.path("..", "starlng_run")
dataset_name <- "cao_lung"

gene_info <- read.csv(file.path(apps_folder, paste0(dataset_name, "_starlng_app"), "objects/genes_info.csv"), row.names = 1)
expr_matrix <- rhdf5::h5read(file.path(apps_folder, paste0(dataset_name, "_starlng_app"), "objects/expression.h5"), "expression_matrix")
colnames(expr_matrix) <- as.character(rhdf5::h5read(file.path(apps_folder, paste0(dataset_name, "_starlng_app"), "objects/expression.h5"), "cells"))
rownames(expr_matrix) <- as.character(rhdf5::h5read(file.path(apps_folder, paste0(dataset_name, "_starlng_app"), "objects/expression.h5"), "genes"))
pseudotime_vals <- qs2::qs_read(file.path(apps_folder, paste0(dataset_name, "_starlng_app"), "objects/recommended_pseudotime.qs2"), nthreads = 30)$recommended_pseudotime
cell_umap <- qs2::qs_read(file.path(apps_folder, paste0(dataset_name, "_starlng_app"), "objects/metadata.qs2"), nthreads = 30)
cell_umap <- cell_umap[, c(ncol(cell_umap)-1, ncol(cell_umap)), drop = FALSE]
cell_umap <- cell_umap[names(pseudotime_vals), , drop = FALSE]

# pseudotime iqr per gene - all expressed cells, top 50% expressed cells, top 25% expressed cells
expr_matrix <- expr_matrix[rownames(gene_info), names(pseudotime_vals), drop = FALSE]
gene_info$pseudotime_iqr_all_expressed <- apply(expr_matrix, 1, function(x) {
    expressed_cells <- which(x > 0)
    if (length(expressed_cells) < 2) {
        return(NA)
    }
    return(IQR(pseudotime_vals[expressed_cells]))
})
gene_info$pseudotime_iqr_top50_expressed <- apply(expr_matrix, 1, function(x) {
    expressed_cells <- which(x > 0)
    if (length(expressed_cells) < 2) {
        return(NA)
    }
    top50_threshold <- quantile(x[expressed_cells], probs = 0.5)
    top50_cells <- expressed_cells[which(x[expressed_cells] >= top50_threshold)]
    if (length(top50_cells) < 2) {
        return(NA)
    }
    return(IQR(pseudotime_vals[top50_cells]))
})
gene_info$pseudotime_iqr_top25_expressed <- apply(expr_matrix, 1, function(x) {
    expressed_cells <- which(x > 0)
    if (length(expressed_cells) < 2) {
        return(NA)
    }
    top25_threshold <- quantile(x[expressed_cells], probs = 0.75)
    top25_cells <- expressed_cells[which(x[expressed_cells] >= top25_threshold)]
    if (length(top25_cells) < 2) {
        return(NA)
    }
    return(IQR(pseudotime_vals[top25_cells]))
})

# gene_info$umap_distance_all_expressed <- apply(expr_matrix, 1, function(x) {
#     expressed_cells <- which(x > 0)
#     if (length(expressed_cells) < 2) {
#         return(NA_real_)
#     }
#     avg_dist <- dist(cell_umap[expressed_cells, , drop = FALSE])
#     return(IQR(avg_dist))
# })
gene_info$umap_distance_top50_expressed <- apply(expr_matrix, 1, function(x) {
    expressed_cells <- which(x > 0)
    if (length(expressed_cells) < 2) {
        return(NA_real_)
    }
    top50_threshold <- quantile(x[expressed_cells], probs = 0.5)
    top50_cells <- expressed_cells[which(x[expressed_cells] >= top50_threshold)]
    if (length(top50_cells) < 2) {
        return(NA_real_)
    }
    
    avg_dist <- dist(cell_umap[top50_cells, , drop = FALSE])
    return(IQR(avg_dist))
})
# gene_info$umap_distance_top25_expressed <- apply(expr_matrix, 1, function(x) {
#     expressed_cells <- which(x > 0)
#     if (length(expressed_cells) < 2) {
#         return(NA_real_)
#     }
#     top25_threshold <- quantile(x[expressed_cells], probs = 0.75)
#     top25_cells <- expressed_cells[which(x[expressed_cells] >= top25_threshold)]
#     if (length(top25_cells) < 2) {
#         return(NA_real_)
#     }
#     avg_dist <- dist(cell_umap[top25_cells, , drop = FALSE])
#     return(IQR(avg_dist))
# })

# bins
# split the genes into more even bins based on the distribution of Moran's I values
bins <- seq(-1, 1, by = 0.05)
# bins <- quantile(gene_info$morans_I, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
# bins[1] <- bins[1] - 1e-5
# bins[length(bins)] <- bins[length(bins)] + 1e-5
gene_info$morans_I_bin <- cut(gene_info$morans_I, breaks = bins, include.lowest = TRUE)
gene_info_summary <- gene_info %>%
    group_by(morans_I_bin) %>%
    summarise(count = n()) %>%
    ungroup()


plotlist <- list(
    ggplot(gene_info_summary, aes(x = morans_I_bin, y = count)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_y_log10() +
        xlab("Moran's I bins") +
        ylab("Number of genes") +
        ggtitle(paste0("Moran's I distribution - ", dataset_name)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(gene_info, aes(x = morans_I_bin, y = umap_distance_top50_expressed)) +
        geom_boxplot() +
        theme_bw() +
        xlab("Moran's I bins") +
        ylab("UMAP Distance (top 50% expressed cells)") +
        ggtitle(paste0("Moran's I vs UMAP Distance IQR - ", dataset_name )) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(gene_info, aes(x = morans_I_bin, y = pseudotime_iqr_top50_expressed)) +
        geom_boxplot() +
        theme_bw() +
        xlab("Moran's I bins") +
        ylab("Pseudotime IQR (top 50% expressed cells)") +
        ggtitle(paste0("Moran's I vs Pseudotime IQR - ", dataset_name)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

combined_plot <- cowplot::plot_grid(plotlist = plotlist, nrow = 1)

pdf(paste0("moran_i_score_distribution_", dataset_name, ".pdf"), width = 15, height = 5)
print(combined_plot)
dev.off()
