library(dplyr)
library(ggplot2)
library(tidyr)
library(Starlng)
library(patchwork)

if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}

input_dir <- "comparison_files"
output_dir <- "panels"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

enrichment_all <- qs2::qs_read(file.path(input_dir, "4_enrichment_analysis.qs2"), nthreads = 30)
enrichment_top <- qs2::qs_read(file.path(input_dir, "4_enrichment_top_10.qs2"), nthreads = 30)
outlier_assessment <- qs2::qs_read(file.path(input_dir, "3_module_stats_summary.qs2"), nthreads = 30)


method_colours <- qualpalr::qualpal(7)$hex
names(method_colours) <- c("Starlng", "Hotspot", "hdWGCNA", "SCENIC", "PAGA+\nCellRank", "Starlng\n(low Moran's I)", "Starlng\n(high Moran's I)")
max_n_modules <- 10

for (dts_name in names(enrichment_all)) {
    if (dts_name %in% c("cao_Liver", "cao_Lung")) {
        next
    }
    all_terms_df <- NULL
    top_terms_df <- NULL
    enrichment_plots <- list()
    for (method_name in names(enrichment_all[[dts_name]])) {
        considered_modules <- outlier_assessment[[dts_name]][[method_name]] %>%
            filter(is_outlier == "no") %>%
            arrange(median_pseudotime) %>%
            pull(module)
        # select max_n_modules and make sure you cover the entire range


        enrich_method <- enrichment_all[[dts_name]][[method_name]]
        searched_modules <- intersect(considered_modules, names(enrich_method))
        if (length(searched_modules) == 0) {
            next
        }
        if (length(searched_modules) > max_n_modules) {
            cuts <- seq(1, length(searched_modules), length.out = max_n_modules + 1)
            searched_modules <- searched_modules[round(cuts)]
        }

        method_enrich_df <- filter_enrichment_results(enrich_method[searched_modules])
        enrichment_plots[[method_name]] <- plot_enrichment_top_terms(
            method_enrich_df %>% filter(source == "GO:BP", term_size <= 500),
            colour_column = "intersection_size",
            top_n = 2,
            point_size_range = c(3, 8),
            words_per_line = 7


        ) + ggtitle(paste0(dts_name, " - ", method_name, " top 2 enriched GO:BP terms")) +
        labs(
            y = "",
        ) +
        coord_flip() +
        theme(
            axis.text.x = element_text(angle = 15, hjust = 0)
        ) +
        scale_y_discrete(position = "right") +
        theme(
            legend.box = "horizontal"
        ) +
        scale_colour_gradient(name = "-log10\n(p-value)") +
        scale_size_continuous(name = "Intersect\nSize")
        
    }

    pdf(file.path(output_dir, paste0("sup_enrichment_dotplots_", dts_name, ".pdf")), width = 14, height = 20)
    print(wrap_plots(enrichment_plots, ncol = 1))
    dev.off()
}

# pdf(file.path(output_dir, "enrichment_dotplots.pdf"), width = 16, height = 12)
# dev.off()