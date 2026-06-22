# single-cell analysis package
library(Seurat)
library(ClustAssess)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

enableWGCNAThreads(5)

if (basename(getwd()) != "comparison_with_other_methods") {
    setwd("comparison_with_other_methods")
}
output_dir <- "hdwgcna_results"
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

dts_path <- "../data"
dts_prefix <- "cao_Liver_subset"
var_group <- "seurat_clusters"
so <- qs2::qs_read(file.path(dts_path, paste0(dts_prefix, "_filtered_normalized.qs2")))

set.seed(42)
so <- SetupForWGCNA(
    so,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = dts_prefix
)

selected_meta <- var_group
so <- MetacellsByGroups(
    seurat_obj = so,
    group.by = c(selected_meta),
    reduction = "pca",
    layer = "data",
    slot = "data",
    k = 25,
    max_shared = 10,
    ident.group = selected_meta,
    wgcna_name = dts_prefix

)

so <- NormalizeMetacells(so)
so <- SetDatExpr(
    seurat_obj = so,
    group.by = selected_meta,
    group_name = so@misc[[dts_prefix]]$wgcna_params$metacell_stats[[selected_meta]],
    assay = "RNA",
    slot = "data",
    layer = "data"
)

so <- TestSoftPowers(
    seurat_obj = so,
    networkType = "signed"
)

# plot_list <- PlotSoftPowers(so)
# wrap_plots(plot_list, ncol = 2)

power_table <- GetPowerTable(so)
head(power_table)

is_power_inf <- TRUE
try({
    so <- ConstructNetwork(
        seurat_obj = so,
        tom_name = dts_prefix,
        overwrite_tom = TRUE,
        wgcna_name = dts_prefix
    )
    is_power_inf <- FALSE
}, silent = TRUE)
if (is_power_inf) {
    if (!is.null(power_table) && nrow(power_table) > 0) {
        ma <- power_table$Power[which.min(power_table$median.k.)]
    } else {
        ma <- 3
    }
    so <- ConstructNetwork(
        seurat_obj = so,
        soft_power = ma,
        tom_name = dts_prefix,
        overwrite_tom = TRUE,
        wgcna_name = dts_prefix
    )
}

PlotDendrogram(so)
 
so <- ModuleEigengenes(so)
so <- ModuleConnectivity(so)

PlotKMEs(so, ncol = 3)
GetModules(so) %>% filter(module != "grey") %>% nrow
write.csv(GetHubGenes(so, n_hubs = -1), file = file.path(output_dir, paste0(dts_prefix, "_modules.csv")), row.names = FALSE)
so <- ModuleExprScore(so)

plot_list <- ModuleFeaturePlot(
    so,
    feature = 'hMEs',
    order = TRUE
)

patchwork::wrap_plots(plot_list, ncol = 5)

