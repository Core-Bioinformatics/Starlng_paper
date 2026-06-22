library(dplyr)
library(gprofiler2)
library(Starlng)

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
            input_file <- file.path("scenic", paste0("sergio_", generation_type, count_type, "_", n_tf, "_tfs_reg.csv"))
            reg_df <- read.csv(input_file)
            reg_df <- reg_df[-(1:2), ]
            tfs <- reg_df[,1]
            associated_genes <- reg_df$Enrichment.6

            # associated_genes <- lapply(associated_genes, function(x) {
            tf_gene_mapping <- do.call(rbind, lapply(seq_along(tfs), function(i) {
                tf <- tfs[i]
                x <- associated_genes[i]
                atomic_gene <- do.call(rbind, lapply(strsplit(x, "),")[[1]], function(y) {
                    gene_name <- strsplit(y, "'")[[1]][2]
                    score_name <- strsplit(y, "',")[[1]][2]
                    data.frame(tf_name = tf, gene_name = gene_name, weight = score_name)
                }))

                return(atomic_gene)
            }))

            # for each gene, keep the interaction with highest wieight
            tf_gene_mapping_no_overlap <- tf_gene_mapping %>%
                group_by(gene_name) %>%
                slice_max(weight, n = 1, with_ties = FALSE) %>%
                ungroup()
                
            unique_tf_names <- setdiff(unique(tf_gene_mapping_no_overlap$tf_name), unique(tf_gene_mapping_no_overlap$gene_name))
            if (length(unique_tf_names) > 0) {
                tf_gene_mapping_no_overlap <- rbind(
                    tf_gene_mapping_no_overlap,
                    data.frame(tf_name = unique_tf_names, gene_name = unique_tf_names, weight = NA)
                )
            }

            cluster_output[[as.character(n_tf)]][[count_type]][[generation_type]] <- tf_gene_mapping_no_overlap
        }
    }
}

qs2::qs_save(cluster_output, "scenic_cluster_output.qs2")
