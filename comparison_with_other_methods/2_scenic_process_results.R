library(dplyr)

# dts_names <- c(paste0("cao_", c("Liver", "Pancreas", "Lung")), "masld_immune")
dts_names <- "cao_Liver_subset"
if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}
results_folder <- "scenic_results"

for (dts_name in dts_names) {
    reg_df <- read.csv(file.path(results_folder, paste0(dts_name, "_reg.csv")))
    reg_df <- reg_df[-(1:2), ]
    tfs <- reg_df[,1]
    associated_genes <- reg_df$Enrichment.6
    tf_gene_mapping <- do.call(rbind, lapply(seq_along(tfs), function(i) {
        tf <- tfs[i]
        x <- associated_genes[i]
        atomic_gene <- do.call(rbind, lapply(strsplit(x, "),")[[1]], function(y) {
            gene_name <- strsplit(y, "'")[[1]][2]
            score_name <- strsplit(y, "',")[[1]][2]
            score_name <- strsplit(score_name, "]")[[1]][1]
            score_name <- strsplit(score_name, ")")[[1]][1]
            score_name <- gsub("\\[", "", score_name)
            score_name <- gsub("\\)", "", score_name)
            score_name <- as.numeric(score_name)
            data.frame(tf_name = tf, gene_name = gene_name, weight = score_name)
        }))

        return(atomic_gene)
    }))

    tf_gene_mapping_no_overlap <- tf_gene_mapping %>%
        group_by(gene_name) %>%
        slice_max(weight, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        as.data.frame() %>%
        arrange(tf_name, gene_name)
    tf_gene_mapping_no_overlap$tf_name <- as.integer(factor(tf_gene_mapping_no_overlap$tf_name))
    colnames(tf_gene_mapping_no_overlap) <- c("module", "gene", "score")
    

    write.csv(tf_gene_mapping_no_overlap[, c("gene", "module", "score")], file.path(results_folder, paste0(dts_name, "_processed.csv")), row.names = FALSE)
}
