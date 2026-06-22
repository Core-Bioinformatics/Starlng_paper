library(dplyr)

if (basename(getwd()) != "data") {
    setwd("data")
}

dts_name <- "cao_Liver_subset"
starlng_run_path <- file.path("..", "starlng_run", paste0(dts_name, "_starlng_app"), "objects")
recommended_pseudotime <- qs2::qs_read(file.path(starlng_run_path, "recommended_pseudotime.qs2"), nthreads = 30)

write.csv(recommended_pseudotime, file = paste0(dts_name, "_recommended_pseudotime.csv"), row.names = TRUE)
