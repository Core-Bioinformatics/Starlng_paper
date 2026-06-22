import hotspot
import scanpy as sc
import numpy as np
import h5py
import os
import pandas as pd
from scipy.sparse import csc_matrix
import json

data_folder = "/mnt/d/Starlng_paper/data"

n_tfs = [5, 10, 15, 20]

results_dictionary = {}
for n_tf in n_tfs:
    results_dictionary[n_tf] = {}
    for count_type in ["raw", "noisy"]:
        results_dictionary[n_tf][count_type] = {}
        for generation_type in ["", "dynamic_ds12_"]:
            input_file = os.path.join(data_folder, f"sergio_{generation_type}{count_type}_{n_tf}_tfs_seurat.h5ad")
            adata = sc.read_h5ad(input_file)
            adata.obs["total_counts"] = adata.X.sum(axis=1)
            adata.layers["counts"] = adata.X.copy()
            sc.pp.filter_genes(adata, min_cells = 10)
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            adata.layers["log_normalized"] = adata.layers["counts"]
            sc.pp.scale(adata)
            sc.pp.pca(adata, mask_var = None)

            adata.layers["counts"] = adata.layers["counts"].todense()
            hs = hotspot.Hotspot(
                adata,
                layer_key = "counts",
                model = "danb",
                latent_obsm_key = "X_pca",
                umi_counts_obs_key = "total_counts"
            )

            np.random.seed(42)
            hs.create_knn_graph(weighted_graph = False, n_neighbors = 30)
            hs_results = hs.compute_autocorrelations(jobs = 4)

            hs_genes = hs_results.loc[hs_results.FDR < 0.05].sort_values('Z', ascending=False)
            print(f"Number of genes with FDR < 0.05: {hs_genes.shape[0]}")
            hs_genes = hs_genes.index
            lcz = hs.compute_local_correlations(hs_genes, jobs = 30)

            modules = hs.create_modules(
                min_gene_threshold = 20,
                fdr_threshold = 0.05,
                core_only = True
            )
            
            results_dictionary[n_tf][count_type][generation_type] = {
                "genes": modules.index.tolist(),
                "clusters": modules.tolist()
            }
# write to json
output_path = "/mnt/d/Starlng_paper/comparison_against_ground_truth/hotspot_cluster_output.json"
with open(output_path, "w") as f:
    json.dump(results_dictionary, f)
