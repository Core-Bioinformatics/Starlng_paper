library(dplyr)
library(ggplot2)
library(Starlng)
library(patchwork)
library(Seurat)
library(grid)

create_intervals <- function(stats_df, min_psd, max_psd) {
    stats_df <- stats_df %>% arrange(.data$median_psd)

    intervals <- list()
    intervals[[as.character(sprintf("%.4f", min_psd))]] <- c(as.character(sprintf("%.4f", max_psd)), "0")

    n_intervals <- nrow(stats_df)
    to_fit_intervals <- lapply(seq_len(n_intervals), function(i) {
        c(stats_df$q1[i], stats_df$q3[i])
    })

    while (n_intervals > 0) {
        q1 <- as.numeric(sprintf("%.4f", to_fit_intervals[[1]][1]))
        q3 <- as.numeric(sprintf("%.4f", to_fit_intervals[[1]][2]))
        to_fit_intervals[[1]] <- NULL
        n_intervals <- n_intervals - 1

        for (curr_q1 in names(intervals)) {
            curr_q3 <- as.numeric(intervals[[curr_q1]][1])
            index <- as.integer(intervals[[curr_q1]][2])
            curr_q1 <- as.numeric(curr_q1)

            # check if q1 - q3 intersects with curr_q1 - curr_q3
            if (q3 <= curr_q1 || q1 >= curr_q3) {
                next
            }

            # print(glue::glue("We have intersection between {q1}-{q3} and {curr_q1}-{curr_q3}"))

            intervals[[sprintf("%.4f", curr_q1)]] <- NULL

            if (q1 == curr_q1 && q3 == curr_q3) {
                intervals[[sprintf("%.4f", q1)]] <- c(sprintf("%.4f", q3), as.character(index + 1))
                next
            }

            if (q1 < curr_q1) {
                if (q3 < curr_q3) {
                    # intervals[[sprintf("%.4f", q1)]] <- c(sprintf("%.4f", curr_q1), as.character(index))
                    intervals[[sprintf("%.4f", curr_q1)]] <- c(sprintf("%.4f", q3), as.character(index + 1))
                    intervals[[sprintf("%.4f", q3)]] <- c(sprintf("%.4f", curr_q3), as.character(index))

                    n_intervals <- n_intervals + 1
                    to_fit_intervals[[n_intervals]] <- c(q1, curr_q1)
                } else {
                    # intervals[[sprintf("%.4f", q1)]] <- c(sprintf("%.4f", curr_q1), as.character(index))
                    intervals[[sprintf("%.4f", curr_q1)]] <- c(sprintf("%.4f", curr_q3), as.character(index + 1))

                    n_intervals <- n_intervals + 1
                    to_fit_intervals[[n_intervals]] <- c(q1, curr_q1)

                    if (q3 != curr_q3) {
                        # intervals[[sprintf("%.4f", curr_q3)]] <- c(sprintf("%.4f", q3), as.character(index))

                        n_intervals <- n_intervals + 1
                        to_fit_intervals[[n_intervals]] <- c(curr_q3, q3)
                    }
                }
                break
            }

            if (q1 > curr_q1) {
                if (q3 < curr_q3) {
                    intervals[[sprintf("%.4f", curr_q1)]] <- c(sprintf("%.4f", q1), as.character(index))
                    intervals[[sprintf("%.4f", q1)]] <- c(sprintf("%.4f", q3), as.character(index + 1))
                    intervals[[sprintf("%.4f", q3)]] <- c(sprintf("%.4f", curr_q3), as.character(index))

                    # n_intervals <- n_intervals + 1
                    # to_fit_intervals[[n_intervals]] <- c(curr_q1, q1)
                } else {
                    intervals[[sprintf("%.4f", curr_q1)]] <- c(sprintf("%.4f", q1), as.character(index))
                    intervals[[sprintf("%.4f", q1)]] <- c(sprintf("%.4f", curr_q3), as.character(index + 1))

                    # n_intervals <- n_intervals + 1
                    # to_fit_intervals[[n_intervals]] <- c(curr_q1, q1)

                    if (q3 != curr_q3) {
                        # intervals[[sprintf("%.4f", curr_q3)]] <- c(sprintf("%.4f", q3), as.character(index))

                        n_intervals <- n_intervals + 1
                        to_fit_intervals[[n_intervals]] <- c(curr_q3, q3)
                    }
                }
                break
            }

            if (q1 == curr_q1) {
                if (q3 < curr_q3) {
                    intervals[[sprintf("%.4f", q1)]] <- c(sprintf("%.4f", q3), as.character(index + 1))
                    intervals[[sprintf("%.4f", q3)]] <- c(sprintf("%.4f", curr_q3), as.character(index))
                } else {
                    intervals[[sprintf("%.4f", q1)]] <- c(sprintf("%.4f", curr_q3), as.character(index + 1))
                    # intervals[[sprintf("%.4f", curr_q3)]] <- c(sprintf("%.4f", q3), as.character(index))

                    n_intervals <- n_intervals + 1
                    to_fit_intervals[[n_intervals]] <- c(curr_q3, q3)
                }
            }
            break
        }
    }

    order_intervals <- order(as.numeric(names(intervals)))
    return(intervals[order_intervals])
}

compute_jumps <- function(stat_df, min_psd) {
    stat_df <- stat_df %>% arrange(.data$median_psd)
    k <- nrow(stat_df)
    jumps <- rep(NA_real_, k)

    for (i in seq_len(k)) {
        current_q1 <- stat_df$q1[i]
        current_q3 <- stat_df$q3[i]
        prev_q3 <- ifelse(i == 1, min_psd, stat_df$q3[i - 1])
        prev_q1 <- ifelse(i == 1, min_psd, stat_df$q1[i - 1])

        if (prev_q3 < current_q1) {
            jumps[i] <- current_q1 - prev_q3
        } else {
            jumps[i] <- (current_q1 + current_q3 - prev_q3 - prev_q1) / 2
        }
        # if (i == 1) {
        #     jumps[i] <- stat_df$q1[i] - min_psd
        # } else {
        #     jumps[i] <- stat_df$q1[i] - stat_df$q3[i - 1]
        # }
    }

    jumps
}

plot_jumps <- function(stat_df, min_psd) {
    stat_df <- stat_df %>% arrange(.data$median_psd)
    k <- nrow(stat_df)
    jumps <- compute_jumps(stat_df, min_psd)

    sum_jumps <- sum(abs(jumps))
    ggplot(stat_df, aes(x = seq_len(nrow(stat_df)), y = jumps)) +
        geom_histogram(binwidth = 0.1, stat = "identity") +
        annotate("text", label = paste0("Sum: ", format(round(sum_jumps, 2), nsmall = 2)), x = nrow(stat_df) / 2.5, y = max(jumps) * 1.1, hjust = 0, vjust = 1, size = 3) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 10),
            legend.position = "none"
        ) + xlab("") 
}

compute_n_hits <- function(stats_df, intervals, pseudotime_vals) {
    stats_df <- stats_df %>% arrange(.data$median_psd)
    n_hits <- list()
    for (interv_start in names(intervals)) {
        interv_end <- intervals[[interv_start]][1]
        index <- intervals[[interv_start]][2]

        n_cells <- sum(pseudotime_vals > as.numeric(interv_start) & pseudotime_vals <= as.numeric(interv_end))

        if (index %in% names(n_hits)) {
            n_hits[[index]] <- n_hits[[index]] + n_cells
        } else {
            n_hits[[index]] <- n_cells
        }
    }

    total_hits <- sum(unlist(n_hits))

    for (index in names(n_hits)) {
        n_hits[[index]] <- n_hits[[index]] / total_hits
    }

    unlist(n_hits)
}

expand_ylim <- function(ylim, pad_fraction = 0.05, clamp_lower = NULL) {
    if (any(!is.finite(ylim))) {
        return(NULL)
    }

    if (ylim[1] == ylim[2]) {
        pad <- ifelse(ylim[1] == 0, 0.1, abs(ylim[1]) * pad_fraction)
        ylim <- c(ylim[1] - pad, ylim[2] + pad)
    } else {
        pad <- diff(ylim) * pad_fraction
        ylim <- c(ylim[1] - pad, ylim[2] + pad)
    }

    if (!is.null(clamp_lower)) {
        ylim[1] <- max(clamp_lower, ylim[1])
    }

    ylim
}


plot_n_hits <- function(stats_df, intervals, pseudotime_vals) {
    stats_df <- stats_df %>% arrange(.data$median_psd)
    n_hits <- list()
    for (interv_start in names(intervals)) {
        interv_end <- intervals[[interv_start]][1]
        index <- intervals[[interv_start]][2]

        n_cells <- sum(pseudotime_vals > as.numeric(interv_start) & pseudotime_vals <= as.numeric(interv_end))

        if (index %in% names(n_hits)) {
            n_hits[[index]] <- n_hits[[index]] + n_cells
        } else {
            n_hits[[index]] <- n_cells
        }
    }

    total_hits <- sum(unlist(n_hits))

    for (index in names(n_hits)) {
        n_hits[[index]] <- n_hits[[index]] / total_hits
    }

    # print(n_hits)
    # print(sum(unlist(n_hits)))
    # print(nrow(temp_df))

    ggplot(data.frame(x = as.integer(names(n_hits)), y = unlist(n_hits)), aes(x = .data$x, y = .data$y)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 10),
        legend.position = "none"
    ) +
    xlab("#module overlaps") +
    ylab("% cells")
}



if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}

input_dir <- "comparison_files"
output_dir <- "panels"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}
comparison_table_file <- file.path("..", "comparison_table.csv")
if (!file.exists(comparison_table_file)) {
    comparison_table <- data.frame(
        method_name = character(),
        dataset_name = character(),
        comparison_type = character(),
        value = numeric(),
        rank = numeric()
    )
} else {
    comparison_table <- read.csv(comparison_table_file, stringsAsFactors = FALSE)
    # remove rows with comparison_type starting with "enrichment_all_modules_" or "enrichment_top10hubs_modules_"
    comparison_table <- comparison_table %>%
        filter(!(comparison_type %in% c(
            "pseudotime_jumps"
        )))
}

module_stats <- qs2::qs_read(file.path(input_dir, "3_module_stats_summary.qs2"), nthreads = 30)
avg_summs <- qs2::qs_read(file.path(input_dir, "3_avg_summaries.qs2"), nthreads = 30)

trend_plot_list <- list()
heatmap_plot_list <- list()

for (dts_name in names(module_stats)) {
    if (dts_name %in% c("cao_Liver", "cao_Lung")) {
        next
    }
    jumps_plot_list <- list()
    n_hits_plot_list <- list()
    vln_plot_list <- list()

    pseudotime_vals <- read.csv(file.path("..", "data", paste0(dts_name, "_recommended_pseudotime.csv")))
    pseudotime_vals <- pseudotime_vals$recommended_pseudotime
    pseudotime_ylim <- range(pseudotime_vals, na.rm = TRUE)
    jumps_ylim <- c(Inf, -Inf)
    n_hits_ylim <- c(Inf, -Inf)

    method_ranking <- NULL
    for (method_name in names(module_stats[[dts_name]])) {
        print(paste0("Processing ", dts_name, " - ", method_name))
        if (method_name %in% trend_plot_list[[dts_name]]) {
            next
        }
        if (startsWith(method_name, "starlng") && dts_name == "masld_immune") {
            prefix <- gsub("starlng_", "", method_name)
            input_file <- file.path(paste0("starlng", "_results"), paste0(dts_name, "_", prefix, "_processed.csv"))
        } else {
            input_file <- file.path(paste0(method_name, "_results"), paste0(dts_name, "_processed.csv"))
        }
        method_stats <- module_stats[[dts_name]][[method_name]]
        selected_modules <- method_stats %>%
            filter(is_outlier == "no") %>%
            arrange(median_pseudotime) %>%
            pull(module)
        method_avg <- avg_summs[[dts_name]][[method_name]][selected_modules]
        method_clust_df <- read.csv(input_file)
        method_clust_df <- method_clust_df %>% filter(module %in% selected_modules)

        module_psd_df_stat <- NULL
        module_psd_df <- NULL
        for (module_name in selected_modules) {
            avg_expr <- avg_summs[[dts_name]][[method_name]][[module_name]]
            module_mask <- avg_expr > 0.5
            module_psd <- pseudotime_vals[module_mask]
            module_psd_df <- rbind(module_psd_df, data.frame(cluster = module_name, pseudotime = module_psd))
            module_psd_df_stat <- rbind(module_psd_df_stat, data.frame(
                cluster = module_name,
                median_psd = median(module_psd),
                q1 = quantile(module_psd, 0.25),
                q3 = quantile(module_psd, 0.75)
            ))
        }
        module_psd_df_stat <- module_psd_df_stat %>% arrange(median_psd)
        module_order_psd <- module_psd_df_stat$cluster
        module_psd_df$cluster <- factor(module_psd_df$cluster, levels = module_order_psd)

        jumps_vals <- compute_jumps(module_psd_df_stat, min(pseudotime_vals))
        n_hits_vals <- compute_n_hits(module_psd_df_stat, create_intervals(module_psd_df_stat, min(pseudotime_vals), max(pseudotime_vals)), pseudotime_vals)
        jumps_ylim[1] <- min(jumps_ylim[1], min(jumps_vals, na.rm = TRUE))
        jumps_ylim[2] <- max(jumps_ylim[2], max(jumps_vals, na.rm = TRUE))
        n_hits_ylim[1] <- min(n_hits_ylim[1], min(n_hits_vals, na.rm = TRUE))
        n_hits_ylim[2] <- max(n_hits_ylim[2], max(n_hits_vals, na.rm = TRUE))



        # histogram of jumps

        jumps_plot_list[[method_name]] <- plot_jumps(module_psd_df_stat, min(pseudotime_vals)) +
            ggtitle(paste0("Jumps histogram for\n", method_name, " on ", dts_name))
        n_hits_plot_list[[method_name]] <- plot_n_hits(module_psd_df_stat, create_intervals(module_psd_df_stat, min(pseudotime_vals), max(pseudotime_vals)), pseudotime_vals) +
            ggtitle(paste0("Number of hits for\n", method_name, " on ", dts_name))

        vln_plot_list[[method_name]] <- ggplot(module_psd_df, aes(x = .data$cluster, y = .data$pseudotime)) +
            geom_boxplot(outlier.size = 0.1) +
            theme_classic() +
            theme(
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                axis.text.y = element_text(size = 10),
                legend.position = "none"
            ) +
            xlab("Module") +
            ylab("Pseudotime") +
            ggtitle(paste0("Pseudotime distribution for", method_name, " on ", dts_name))

        jumps_vals <- compute_jumps(module_psd_df_stat, min(pseudotime_vals))
        method_ranking <- rbind(method_ranking, data.frame(
            method_name = method_name,
            dataset_name = dts_name,
            comparison_type = "pseudotime_jumps",
            value = sum(abs(jumps_vals)),
            rank = NA
        ))
    }

    row1 <- wrap_plots(vln_plot_list) + plot_layout(nrow = length(vln_plot_list)) &
        coord_cartesian(ylim = expand_ylim(pseudotime_ylim, pad_fraction = 0.03))
    row2 <- wrap_plots(jumps_plot_list) + plot_layout(ncol = length(jumps_plot_list)) &
        coord_cartesian(ylim = expand_ylim(jumps_ylim, pad_fraction = 0.05))
    row3 <- wrap_plots(n_hits_plot_list) + plot_layout(ncol = length(n_hits_plot_list)) &
        coord_cartesian(ylim = expand_ylim(n_hits_ylim, pad_fraction = 0.05, clamp_lower = 0))

    fig <- wrap_plots(row1, row2, row3, nrow = 3, heights = c(4.5, 1, 1)) 
    pdf(paste0("panels/sup_jumps_n_hits_", dts_name, ".pdf"), width = 13, height = 20)
    print(fig)
    dev.off()

    # add the jumps values to the comparison_table
    method_ranking$rank <- rank(method_ranking$value, ties.method = "min")
    method_ranking$method_name <- recode(method_ranking$method_name,
        "starlng" = "Starlng",
        "hotspot" = "Hotspot",
        "hdwgcna" = "hdWGCNA",
        "scenic" = "SCENIC",
        "paga_cellrank" = "PAGA+CellRank",
        "starlng_low_moran" = "Starlng(low Moran's I)",
        "starlng_high_moran" = "Starlng(high Moran's I)"
    )
    comparison_table <- rbind(comparison_table, method_ranking)
}

write.csv(comparison_table, comparison_table_file, row.names = FALSE)
