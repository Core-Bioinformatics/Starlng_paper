import cellrank as cr
import scanpy as sc
import numpy as np
import h5py
import os
import pandas as pd
from scipy.sparse import csc_matrix, csr_matrix

# Compatibility patch for scanpy<=1.9.3 with scipy>=1.13 in PAGA internals.
import scanpy._utils as sc_utils


def _get_sparse_from_igraph_compat(graph, weight_attr=None):
    edges = graph.get_edgelist()
    weights = (
        np.array(graph.es[weight_attr]).astype(np.float64)
        if weight_attr is not None
        else np.ones(len(edges), dtype=np.float64)
    )
    shape = graph.vcount()
    if not graph.is_directed():
        edges = edges + [(v, u) for u, v in edges]
        weights = np.concatenate([weights, weights])
    if len(edges) > 0:
        rows, cols = zip(*edges)
        return csr_matrix((weights, (rows, cols)), shape=(shape, shape))
    return csr_matrix((shape, shape))


sc_utils.get_sparse_from_igraph = _get_sparse_from_igraph_compat

dts_name = "cao_Liver_subset"
data_folder = "/mnt/d/Starlng_paper/data"
output_folder = "/mnt/d/Starlng_paper/comparison_with_other_methods/paga_results"
path = f"{data_folder}/{dts_name}.h5ad"
psd_recommendation = pd.read_csv(f"{data_folder}/{dts_name.lower()}_recommended_pseudotime.csv", index_col=0)
starting_metadata = psd_recommendation["recommended_mtd_name"][0]
starting_group = psd_recommendation["recommended_mtd_group"][0]
grouping_metadata_name = "seurat_clusters" if dts_name != "masld_immune" else "stable_24_clusters"

adata = sc.read_h5ad(path)
adata.X = adata.layers["counts"].copy()
adata.obs["total_counts"] = adata.X.sum(axis=1)
adata.layers["counts"] = adata.X.copy()
sc.pp.filter_cells(adata, min_genes = 5)
sc.pp.filter_genes(adata, min_cells = 5)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["log_normalized"] = adata.X.copy()
sc.pp.scale(adata)

sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30, random_state = 42, use_rep = "pca")
sc.tl.draw_graph(adata, random_state = 42, init_pos = "umap")

sc.tl.paga(adata, groups=grouping_metadata_name)
sc.pl.paga(adata)

adata.uns["iroot"] = np.argmin(np.linalg.norm(adata.obsm["umap"] - adata[adata.obs[starting_metadata] == starting_group].obsm["umap"].mean(axis=0), axis=1))
sc.tl.dpt(adata)

pt_kernel = cr.kernels.PseudotimeKernel(adata, time_key="dpt_pseudotime")
pt_kernel.compute_transition_matrix()

pt_es = cr.estimators.GPCCA(pt_kernel)
pt_es.fit(cluster_key = grouping_metadata_name)

pt_es.set_initial_states({"0": adata.obs_names[adata.obs[starting_metadata] == str(starting_group)]}, cluster_key = starting_metadata, n_cells = 15)
pt_es.predict_terminal_states(allow_overlap=True)
pt_es.compute_fate_probabilities()

adata.X = adata.layers["log_normalized"].copy()
model = cr.models.GAM(adata, n_knots = 6)

if "var.features" in adata.var.keys():
    var_features = adata.var["var.features"].copy()
    mask = var_features.values == "NA"
    var_features = var_features[~mask]
    var_features = var_features.values.tolist()
else:
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", layer = "log_normalized")
    var_features = adata.var[adata.var["highly_variable"]].index.tolist()

lineage_name_term = list(map(str, adata.obs["term_states_fwd"].cat.categories.to_list()))
lineage_name_initial = list(map(str, adata.obs["init_states_fwd"].cat.categories.to_list()))
print(lineage_name_term)
print(lineage_name_initial)
if not lineage_name_term:
    raise ValueError("No terminal lineage names were found in `term_states_fwd`.")

selected_lineage = starting_group if starting_group in lineage_name_term else lineage_name_term[0]

cr.pl.cluster_trends(
    adata,
    model=model,
    lineage=selected_lineage,
    genes=var_features, # the tutorial recommendation
    time_key="dpt_pseudotime",
    n_jobs=30,
    random_state=0,
    clustering_kwargs={"resolution": 1, "random_state": 42},
    # neighbors_kwargs={"n_neighbors": 2, "random_state": 42},
    recompute=True
)

uns_keys = list(adata.uns.keys())
uns_keys = [k for k in uns_keys if k.startswith("lineage_")]
uns_keys = uns_keys[0]

cluster_mean_combined = adata.uns[uns_keys].obs
mean_info = adata.var.loc[cluster_mean_combined.index]
cluster_mean_combined["mean"] = mean_info["mean"].values

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
cluster_mean_combined.to_csv(f"{output_folder}/{dts_name.lower()}_modules.csv", index=True)