library(dplyr)
library(igraph)

if (basename(getwd()) != "data") {
    setwd("data")
}

tf_gene_mapping <- read.csv("../comparison_with_other_methods/scenic_results/tf_gene_mapping_no_overlap.csv")
tf_gene_mapping$weight <- sapply(tf_gene_mapping$weight, function(x) {
    round(as.numeric(strsplit(strsplit(x, ")")[[1]][1], "]")[[1]][1]), 4)
})
tf_gene_mapping <- tf_gene_mapping %>%
    filter(weight > 0)
ngenes_per_tf <- tf_gene_mapping %>%
    group_by(tf_name) %>%
    summarise(n_genes = length(gene_name)) %>%
    ungroup() %>%
    arrange(desc(n_genes)) %>%
    filter(n_genes >= 10)
ngenes_per_tf

for (n_tfs in c(5, 10, 15, 20)) {
    tf_select_prob <- 1 / ngenes_per_tf$n_genes^0.5
    tf_select_prob <- tf_select_prob / sum(tf_select_prob)

    set.seed(42)
    selected_tfs <- sample(ngenes_per_tf$tf_name, size = n_tfs, replace = FALSE, prob = tf_select_prob)


    tf_gene_mapping_selected <- tf_gene_mapping %>%
        filter(tf_name %in% selected_tfs)
    gene_names <- setdiff(tf_gene_mapping_selected$gene_name, tf_gene_mapping_selected$tf_name)
    tf_gene_mapping_selected <- tf_gene_mapping_selected %>%
        filter(gene_name %in% gene_names)
    all_genes <- c(gene_names, selected_tfs)
    all_genes <- setNames(seq_along(all_genes) - 1, all_genes)

    tf_gene_mapping_selected$gene_name <- all_genes[tf_gene_mapping_selected$gene_name]
    tf_gene_mapping_selected$tf_name <- all_genes[tf_gene_mapping_selected$tf_name]
    tf_gene_mapping_selected$nregs <- 1
    tf_gene_mapping_selected$hill <- 2


    write.table(tf_gene_mapping_selected[, c("gene_name", "nregs", "tf_name", "weight", "hill")], paste0("input_sergio_file_targets_", n_tfs, "_tfs.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
    saveRDS(
        list(
            index = all_genes,
            tf_mapping = tf_gene_mapping_selected
        ),
        paste0("sergio_ground_truth_", n_tfs, "_tfs.rds")
    )

    # first approach - one cluster per tf, very clear separation
    reg_cluster_df <- matrix(0, nrow = n_tfs, ncol = n_tfs + 1)
    set.seed(42)
    reg_cluster_df[, 1] <- sample(all_genes[selected_tfs], n_tfs, replace = FALSE)
    for (i in seq_along(selected_tfs)) {
        reg_cluster_df[i, i+1] <- 1
    }
    write.table(reg_cluster_df, paste0("input_sergio_file_hard_clusters_", n_tfs, "_tfs.csv"), sep = ",", row.names = FALSE, col.names = FALSE)



    # second approach - some tfs will be transitory
    nclusters <- ceiling(0.6 * n_tfs) 
    print(paste0("nclusters: ", nclusters))
    reg_cluster_df <- matrix(0, nrow = n_tfs, ncol = nclusters + 1)
    set.seed(42)
    reg_cluster_df[, 1] <- sample(all_genes[selected_tfs], n_tfs, replace = FALSE)
    for (i in seq_len(nclusters)) {
        for (j in seq_len(nclusters)) {
            if (i == j) {
                reg_cluster_df[i, j+1] <- 3
                next
            } 
            reg_cluster_df[i, j+1] <- 0.5
        }
    }
    for(i in seq(nclusters + 1, n_tfs)) {
        two_indices <- sample(seq_len(nclusters), 2, replace = FALSE)
        reg_cluster_df[i, two_indices + 1] <- 2
        reg_cluster_df[i, -c(1, two_indices + 1)] <- 0.5
    }
    write.table(reg_cluster_df, paste0("input_sergio_file_soft_", nclusters, "_clusters_", n_tfs, "_tfs.csv"), sep = ",", row.names = FALSE, col.names = FALSE)

    nlinks <- colSums(reg_cluster_df[, -1] == 2)
    g <- igraph::make_empty_graph(n = nclusters, directed = FALSE)
    while (sum(nlinks) > 0) {
        cluster_idx <- which.max(nlinks)
        
        link_prob <- nlinks / sum(nlinks)
        link_prob[cluster_idx] <- 0
        linked_clusters <- E(g)[.from(cluster_idx)]$to
        link_prob[linked_clusters] <- 0
        if (sum(link_prob) == 0) {
            break
        }
        link_prob <- link_prob / sum(link_prob)
        while (TRUE) {
            link_idx <- sample(seq_len(nclusters), size = 1, prob = link_prob)
            test_g <- igraph::add_edges(g, c(cluster_idx, link_idx))
            if (igraph::girth(test_g)$girth == Inf) {
                break
            }
            link_prob[link_idx] <- 0
            if (sum(link_prob) == 0) {
                break
            }
            link_prob <- link_prob / sum(link_prob)
        }

        if (sum(link_prob) != 0) {
            g <- igraph::add_edges(g, c(cluster_idx, link_idx))
        }

        nlinks[cluster_idx] <- nlinks[cluster_idx] - 1
    }
    # remove duplicate edges
    g <- igraph::simplify(g)
    while (igraph::count_components(g) > 1) {
        components <- igraph::components(g)$membership
        comp_table <- table(components)
        largest_comp <- as.numeric(names(comp_table)[which.max(comp_table)])
        other_comps <- setdiff(unique(components), largest_comp)
        for (comp in other_comps) {
            comp_nodes <- which(components == comp)
            link_idx <- sample(which(components == largest_comp), size = 1)
            g <- igraph::add_edges(g, c(comp_nodes[1], link_idx))
        }
    }
    E(g)$weight <- runif(igraph::ecount(g), min = 0.2, max = 1.1)
    plot(g)

    g_adj <- igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE, type = "upper")
    write.table(g_adj, paste0("input_sergio_file_grn_", nclusters, "_clusters_", n_tfs, "_tfs.tab"), sep = "\t", row.names = FALSE, col.names = FALSE)
}
