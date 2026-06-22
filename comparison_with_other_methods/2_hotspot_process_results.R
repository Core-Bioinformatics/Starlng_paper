library(dplyr)

# dts_names <- c(paste0("cao_", c("Liver", "Pancreas", "Lung")), "masld_immune")
dts_names <- c("cao_Liver_subset")
if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}
results_folder <- "hotspot_results"

for (dts_name in dts_names) {
    modules <- read.csv(file.path(results_folder, paste0(dts_name, "_hotspot_modules.csv")))
    modules <- modules[, c("gene", "module", "C")]
    colnames(modules) <- c("gene", "module", "score")
    modules <- modules %>% dplyr::filter(module != -1)
    modules$module <- as.integer(factor(modules$module))
    write.csv(modules, file.path(results_folder, paste0(dts_name, "_processed.csv")), row.names = FALSE)
}