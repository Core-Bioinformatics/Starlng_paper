import hotspot
import scanpy as sc
import numpy as np
import h5py
import os
import sys
import pandas as pd
from scipy.sparse import csc_matrix
import anndata as ad
import pickle

data_folder="/mnt/d/Starlng_paper/data"
output_folder="/mnt/d/Starlng_paper/comparison_with_other_methods/hotspot_results"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

dts_name = "cao_Liver_subset"
input_path = os.path.join(data_folder, f"{dts_name}.h5ad")
output_path_pickle = os.path.join(output_folder, f"{dts_name}_hotspot_results.pkl")
output_path_module_csv = os.path.join(output_folder, f"{dts_name}_hotspot_modules.csv")

adata = sc.read_h5ad(input_path)
adata.X = adata.layers["counts"]
adata.obs["total_counts"] = adata.X.sum(axis=1)
adata.layers["counts"] = adata.X.copy()
sc.pp.filter_cells(adata, min_genes = 5)
sc.pp.filter_genes(adata, min_cells = 5)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["log_normalized"] = adata.layers["counts"]
sc.pp.scale(adata)

adata.layers["counts"] = adata.layers["counts"].todense()

hs = hotspot.Hotspot(
    adata,
    layer_key = "counts",
    model = "danb",
    latent_obsm_key = "pca",
    umi_counts_obs_key = "total_counts"
)

np.random.seed(42)
hs.create_knn_graph(weighted_graph = False, n_neighbors = 30)
hs_results = hs.compute_autocorrelations(jobs = 60)

hs_genes = hs_results.loc[hs_results.FDR < 0.05].sort_values('Z', ascending=False)
print(f"Number of genes with FDR < 0.05: {hs_genes.shape[0]}")
hs_genes = hs_genes.index

required_recursion_limit = max(10000, (2 * len(hs_genes)) + 1000)
sys.setrecursionlimit(required_recursion_limit)

lcz = hs.compute_local_correlations(hs_genes, jobs = 60)

modules = hs.create_modules(
    min_gene_threshold = 20,
    fdr_threshold = 0.05,
    core_only = True
)

score_results = hs.results
modules_result = hs.modules
score_results = score_results.loc[modules_result.index]
score_results["module"] = modules_result
score_results["gene"] = score_results.index

score_results.to_csv(output_path_module_csv, index=False)
# with open(output_path_pickle, "wb") as f:
    # pickle.dump(hs, f)
