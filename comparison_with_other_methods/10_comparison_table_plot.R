library(ggplot2)
library(dplyr)
library(tidyr)

if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}

input_file <- file.path("..", "comparison_table.csv")
comparison_table <- read.csv(input_file, stringsAsFactors = FALSE)

real_datasets <- setdiff(unique(comparison_table$dataset_name), "sergio")
target_categories <- c(
    "selected_module_proportion",
    "selected_module_total_coverage",
    "selected_module_max_iqr",
    "enrichment_top10hubs_modules_GO:BP",
    "pseudotime_jumps",
    "ground_truth_discarded_genes",
    "ground_truth_ecs"
)

htmp_dts <- list()
for (dts in real_datasets) {
    dts_list <- c(dts)

    sub_df <- comparison_table %>%
        filter(dataset_name == dts, comparison_type %in% target_categories)

    if (dts == "masld_immune") {
        added_sergio <- comparison_table %>%
            filter(dataset_name == "sergio", comparison_type %in% target_categories)
        added_sergio$method_name[added_sergio$method_name == "Starlng"] <- "Starlng(high Moran's I)"
        added_sergio <- rbind(added_sergio, added_sergio %>% filter(method_name == "Starlng(high Moran's I)") %>% mutate(method_name = "Starlng(low Moran's I)"))
        sub_df <- rbind(sub_df, added_sergio) %>% filter(method_name != "Starlng")
    }
    sub_df$comparison_type <- factor(sub_df$comparison_type, levels = target_categories)
    sub_df$comparison_type <- recode(sub_df$comparison_type,
        "selected_module_proportion" = "Percentage\nselected modules",
        "selected_module_total_coverage" = "Modules\nTotal Coverage",
        "selected_module_max_iqr" = "Modules\nMax IQR",
        "enrichment_top10hubs_modules_TF" = "N TFs\nfrom top 10 hubs",
        "enrichment_top10hubs_modules_GO:BP" = "N GO:BP terms\nfrom top 10 hubs",
        "pseudotime_jumps" = "Pseudotime\nJumps",
        "ground_truth_discarded_genes" = "Synthetic\n% Discarded Genes",
        "ground_truth_ecs" = "Synthetic\nECS to Ground Truth"
    )

    # create a heatmap based on the ranks, show the rank in each cell and in paranthesis the value
    # columns - criteria; rows - methods
    # add a column on the right with the median rank and sum ranks

    summ_rank_df <- sub_df %>%
        group_by(method_name) %>%
        summarise(median_rank = median(rank), sum_rank = sum(rank), .groups = "drop") %>%
        ungroup() 
    method_order <- summ_rank_df %>%
        arrange(median_rank, sum_rank) %>%
        pull(method_name)
    sub_df$method_name <- factor(sub_df$method_name, levels = rev(method_order))
    summ_rank_df$method_name <- factor(summ_rank_df$method_name, levels = rev(method_order))

    main_htmp <- ggplot(sub_df, aes(x = comparison_type, y = method_name, fill = rank)) +
        geom_tile() +
        geom_text(aes(label = paste0(rank, "\n(", round(value, 2), ")")), size = 6) +
        scale_fill_gradient(high = "white", low = "steelblue") +
        theme_minimal() +
        scale_x_discrete(position = "top") +
        theme(legend.position = "bottom") +
        labs(title = paste0("Overall Method Comparison for dataset: ", dts), fill = "Rank", x = "", y = "") +
        theme(
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 15),
            plot.title = element_text(size = 16),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12)
        )
    rnk_median_htmp <- ggplot(summ_rank_df, aes(x = 1, y = method_name, fill = median_rank)) +
        geom_tile() +
        geom_text(aes(label = round(median_rank, 2)), size = 7) +
        theme_void() +
        scale_fill_gradient(high = "#a8d8b8", low = "#045204") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
        labs(title = "Median\nRank", fill = "", x = "", y = "") 
    rnk_sum_htmp <- ggplot(summ_rank_df, aes(x = 1, y = method_name, fill = sum_rank)) +
        geom_tile() +
        geom_text(aes(label = round(sum_rank, 2)), size = 7) +
        theme_void() +
        scale_fill_gradient(high = "#ffe6e6", low = "#af4d4d") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
        labs(title = "Sum\nRank", fill = "", x = "", y = "")

    united_htmp <- patchwork::wrap_plots(main_htmp, rnk_median_htmp, rnk_sum_htmp, ncol = 3, widths = c(7, 1, 1))

    ncols <- length(unique(sub_df$comparison_type))

    pdf(file.path("panels", paste0("main_comparison_table_plot_", dts, ".pdf")), width = ncols * 0.7 + 8, height = 7)
    print(united_htmp)
    dev.off()

    htmp_dts[[dts]] <- united_htmp
}

htmp_all <- patchwork::wrap_plots(htmp_dts, ncol = 1)

# pdf(file.path("panels", "main_comparison_table_plot_all.pdf"), width = 14, height = 20)
# print(htmp_all)
# dev.off()
