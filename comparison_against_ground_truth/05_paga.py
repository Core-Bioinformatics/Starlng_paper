import cellrank as cr
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
            sc.pp.filter_cells(adata, min_genes = 5)
            sc.pp.filter_genes(adata, min_cells = 5)
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            adata.layers["log_normalized"] = adata.X.copy()
            sc.pp.scale(adata)

            sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30, random_state = 42)
            sc.tl.umap(adata, random_state = 42)

            sc.tl.draw_graph(adata, random_state = 42, init_pos = "X_umap")
            sc.tl.paga(adata, groups="gt_clusters")

            adata.uns["iroot"] = np.argmin(np.linalg.norm(adata.obsm["X_umap"] - adata[adata.obs["gt_clusters"] == "0"].obsm["X_umap"].mean(axis=0), axis=1))
            sc.tl.dpt(adata)

            pt_kernel = cr.kernels.PseudotimeKernel(adata, time_key="dpt_pseudotime")
            pt_kernel.compute_transition_matrix()

            pt_es = cr.estimators.GPCCA(pt_kernel)
            pt_es.fit(cluster_key = "gt_clusters")

            pt_es.set_initial_states({"0": adata.obs_names[adata.obs["gt_clusters"] == "cluster1"]}, cluster_key = "gt_clusters", n_cells = 15)
            pt_es.predict_terminal_states(allow_overlap=True)
            pt_es.compute_fate_probabilities()

            adata.X = adata.layers["log_normalized"].copy()
            model = cr.models.GAM(adata, n_knots = 6)

            trend_genes = adata.var_names
            # get the lowest lineage
            lineage_name_term = list(map(str, adata.obs["term_states_fwd"].cat.categories.to_list()))
            lineage_name_initial = list(map(str, adata.obs["init_states_fwd"].cat.categories.to_list()))
            print(lineage_name_term)
            print(lineage_name_initial)
            if not lineage_name_term:
                raise ValueError("No terminal lineage names were found in `term_states_fwd`.")

            selected_lineage = "cluster1" if "cluster1" in lineage_name_term else lineage_name_term[0]
            cr.pl.cluster_trends(
                adata,
                model=model,
                lineage=selected_lineage,
                genes=trend_genes,
                time_key="dpt_pseudotime",
                n_jobs=8,
                random_state=0,
                # pca_kwargs={"n_comps": 1, "svd_solver": "auto", "random_state": 0},
                clustering_kwargs={"resolution": 0.4, "random_state": 0},
                # neighbors_kwargs={"n_neighbors": 2, "random_state": 42},
                recompute=True,
            )

            trend_key = f"lineage_{selected_lineage}_trend"
            results_dictionary[n_tf][count_type][generation_type] = {
                "lineage": selected_lineage,
                "genes": list(map(str, adata.uns[trend_key].obs["clusters"].index.tolist())),
                "clusters": adata.uns[trend_key].obs["clusters"].tolist()
            }


# write to json
output_path = "/mnt/d/Starlng_paper/comparison_against_ground_truth/paga_cellrank_cluster_output.json"
with open(output_path, "w") as f:
    json.dump(results_dictionary, f)
