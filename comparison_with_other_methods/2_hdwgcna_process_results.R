library(dplyr)

# dts_names <- c(paste0("cao_", c("Liver", "Pancreas", "Lung")), "masld_immune")
dts_names <- "cao_Liver_subset"
if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}
results_folder <- "hdwgcna_results"

for (dts_name in dts_names) {
    modules <- read.csv(file.path(results_folder, paste0(dts_name, "_modules.csv")))
    colnames(modules) <- c("gene", "module", "score")
    modules$module <- as.integer(factor(modules$module))
    write.csv(modules, file.path(results_folder, paste0(dts_name, "_processed.csv")), row.names = FALSE)
}