library(dplyr)

# dts_names <- c(paste0("cao_", c("Liver", "Pancreas", "Lung")), "masld_immune_high_moran", "masld_immune_low_moran")
dts_names <- c("cao_Liver_subset")
if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}
results_folder <- "starlng_results"
if (!dir.exists(results_folder)) {
    dir.create(results_folder, recursive = TRUE)
}

for (dts_name in dts_names) {
    app_path <- file.path("..", "starlng_run", paste0(tolower(dts_name), "_starlng_app"), "objects", "module_summaries.h5")
    genes <- as.character(rhdf5::h5read(app_path, "genes"))
    nmodules <- rhdf5::h5read(app_path, "all_modules")
    nmodules <- as.character(max(as.numeric(nmodules)))

    clustering <- rhdf5::h5read(app_path, paste0(nmodules, "/gene_hub_stats"))
    clustering$gene <- genes
    clustering <- clustering[, c("gene", "module", "combined_score")]
    colnames(clustering) <- c("gene", "module", "score")
    write.csv(clustering, file.path(results_folder, paste0(dts_name, "_processed.csv")), row.names = FALSE)
}



