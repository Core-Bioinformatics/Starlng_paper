library(ggplot2)
library(dplyr)
library(Seurat)
library(Starlng)
devtools::load_all("/mnt/d/Starlng")

if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}

# dts_names <- c(paste0("cao_", c("Liver", "Pancreas", "Lung")), "masld_immune")
dts_names <- c("cao_Liver_subset")
# method_names <- c("starlng", "hotspot", "hdwgcna", "paga_cellrank")
method_names <- c("scenic")
output_dir <- "comparison_files"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

output_file_analysis <- file.path(output_dir, "4_enrichment_analysis.qs2")
output_file_top_10 <- file.path(output_dir, "4_enrichment_top_10.qs2")

if (file.exists(output_file_analysis)) {
    enrichment_analysis <- qs2::qs_read(output_file_analysis, nthreads = 30)
} else {
    enrichment_analysis <- list()
}

if (file.exists(output_file_top_10)) {
    enrichment_top_10 <- qs2::qs_read(output_file_top_10, nthreads = 30)
} else {
    enrichment_top_10 <- list()
}
# target_sources <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "TF")
target_sources <- c("GO:BP", "TF")

for (dts_name in dts_names) {
    if (!dts_name %in% names(enrichment_analysis)) {
        enrichment_analysis[[dts_name]] <- list()
        enrichment_top_10[[dts_name]] <- list()
    }

    so <- qs2::qs_read(file.path("..", "data", paste0(dts_name, "_filtered_normalized.qs2")), nthreads = 30)
    expr_matrix <- GetAssayData(so, assay = "RNA", layer = "data")
    background_genes <- rownames(expr_matrix)
    rm(expr_matrix)
    rm(so)
    gc()

    if (dts_name == "masld_immune") {
        actual_methods <- c(method_names[-1], "starlng_high_moran", "starlng_low_moran")
    } else {
        actual_methods <- method_names
    }
    for (method_name in actual_methods) {
        if (startsWith(method_name, "starlng") && dts_name == "masld_immune") {
            prefix <- gsub("starlng_", "", method_name)
            input_file <- file.path(paste0("starlng", "_results"), paste0(dts_name, "_", prefix, "_processed.csv"))
        } else {
            input_file <- file.path(paste0(method_name, "_results"), paste0(dts_name, "_processed.csv"))
        }
        
        module_df <- read.csv(input_file)
        module_names <- stringr::str_sort(unique(module_df$module), numeric = TRUE)
        if (!method_name %in% names(enrichment_analysis[[dts_name]])) {
            enrichment_analysis[[dts_name]][[method_name]] <- list()
            enrichment_top_10[[dts_name]][[method_name]] <- list()
        }

        highest_existing_module <- as.character(max(as.numeric(names(enrichment_analysis[[dts_name]][[method_name]])), na.rm = TRUE))
        module_names <- module_names[as.numeric(module_names) > as.numeric(highest_existing_module)]
        if (length(module_names) == 0) {
            next
        }
        for (module in module_names) {
            print(paste0("Processing ", dts_name, " with method ", method_name, " for module ", module))
            if (module %in% names(enrichment_analysis[[dts_name]][[method_name]]) || module %in% names(enrichment_top_10[[dts_name]][[method_name]])) {
                next
            }
            module_genes <- module_df$gene[module_df$module == module]
            enrichment_result <- gprofiler2::gost(
                query = module_genes,
                organism = "hsapiens",
                sources = target_sources,
                significant = TRUE,
                domain_scope = "custom_annotated",
                custom_bg = background_genes,
                evcodes = TRUE,
                correction_method = "fdr"
            )
            enrichment_analysis[[dts_name]][[method_name]][[module]] <- enrichment_result

            top_10_genes <- module_df[module_df$module == module, ] %>%
                dplyr::arrange(desc(score)) %>%
                dplyr::slice_head(n = 10) %>%
                dplyr::pull(gene)
            enrichment_top_10[[dts_name]][[method_name]][[module]] <- gprofiler2::gost(
                query = top_10_genes,
                organism = "hsapiens",
                sources = target_sources,
                significant = TRUE,
                domain_scope = "custom_annotated",
                custom_bg = background_genes,
                evcodes = TRUE,
                correction_method = "fdr"
            )
            qs2::qs_save(enrichment_analysis, output_file_analysis, nthreads = 30)
            qs2::qs_save(enrichment_top_10, output_file_top_10, nthreads = 30)
        }
    }
}

qs2::qs_save(enrichment_analysis, output_file_analysis, nthreads = 30)
qs2::qs_save(enrichment_top_10, output_file_top_10, nthreads = 30)


