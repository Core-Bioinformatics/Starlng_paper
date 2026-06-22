library(dplyr)
library(gprofiler2)
library(Starlng)

if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}


reg_df <- read.csv("scenic_results/reg.csv")
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

write.csv(tf_gene_mapping, "scenic_results/tf_gene_mapping.csv", row.names = FALSE)

# for each gene, keep the interaction with highest wieight
tf_gene_mapping_no_overlap <- tf_gene_mapping %>%
    group_by(gene_name) %>%
    slice_max(weight, n = 1, with_ties = FALSE) %>%
    ungroup()
write.csv(tf_gene_mapping_no_overlap, "scenic_results/tf_gene_mapping_no_overlap.csv", row.names = FALSE)

# create list of associated genes for each TF, with no overlap between TFs
tf_gene_mapping_no_overlap <- tf_gene_mapping_no_overlap %>%
    group_by(tf_name) %>%
    summarise(genes = list(gene_name)) %>%
    ungroup()
tf_order_names <- tf_gene_mapping_no_overlap$tf_name
tf_gene_mapping_no_overlap <- as.list(tf_gene_mapping_no_overlap$genes)
names(tf_gene_mapping_no_overlap) <- tf_order_names

expr_matrix_path <- file.path("../starlng_run/masld_immune_high_moran_starlng_app", "objects", "expression.h5")
expr_matrix <- rhdf5::h5read(expr_matrix_path, "expression_matrix")
rownames(expr_matrix) <- as.character(rhdf5::h5read(expr_matrix_path, "genes"))
colnames(expr_matrix) <- as.character(rhdf5::h5read(expr_matrix_path, "cells"))

for (i in seq_along(tf_gene_mapping_no_overlap)) {
    tf_gene_mapping_no_overlap[[i]] <- list(
        genes = tf_gene_mapping_no_overlap[[i]]
    )
    temp_enrich <- gprofiler2::gost(
        query = tf_gene_mapping_no_overlap[[i]]$genes,
        organism = "hsapiens",
        sources = c("GO", "TF", "KEGG", "REAC"),
        significant = TRUE,
        domain_scope = "custom",
        custom_bg = as.vector(rownames(expr_matrix)),
        evcodes = TRUE
    )
    if (!is.null(temp_enrich)) {
        temp_enrich <- temp_enrich$result
    }


    tf_gene_mapping_no_overlap[[i]]$enrichment <- temp_enrich
    tf_gene_mapping_no_overlap[[i]]$avg_expr <- voting_scheme(
        expression_matrix = expr_matrix[tf_gene_mapping_no_overlap[[i]]$genes, , drop = FALSE],
        genes = tf_gene_mapping_no_overlap[[i]]$genes,
        thresh_percentile = 0,
        thresh_value = 0,
        n_coexpressed_thresh = 1,
        summary_function = mean
    )
}
str(tf_gene_mapping_no_overlap)

qs2::qs_save(tf_gene_mapping_no_overlap, "scenic_results/scenic_results.qs2", nthreads = 5)
