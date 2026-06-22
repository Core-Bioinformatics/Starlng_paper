import numpy as np
import pandas as pd
import sys

# Compatibility shims for legacy SERGIO code on modern NumPy
if not hasattr(np, "int"):
    np.int = int
if not hasattr(np, "float"):
    np.float = float

import os
os.getcwd()

import os
project_path = os.path.abspath(".")
if project_path not in sys.path:
    sys.path.append(project_path)
from SERGIO.SERGIO.sergio import sergio

def infer_ngenes(input_target, input_regs):
    max_id = -1
    with open(input_target) as f:
        for line in f:
            row = line.strip().split(",")
            if not row or row[0] == "":
                continue
            target_id = int(float(row[0]))
            n_regs = int(float(row[1]))
            reg_ids = [int(float(x)) for x in row[2:2 + n_regs]]
            if reg_ids:
                max_id = max(max_id, target_id, *reg_ids)
            else:
                max_id = max(max_id, target_id)

    with open(input_regs) as f:
        for line in f:
            row = line.strip().split(",")
            if row and row[0] != "":
                max_id = max(max_id, int(float(row[0])))

    return max_id + 1


n_tfs = [5, 10, 15, 20]
ncells_per_cluster = 400

for n_tf in n_tfs:
    soft_n_clusters = int(np.ceil(n_tf * 0.6))
    
    for prefix_regs in [f"hard_clusters_{n_tf}_tfs", f"soft_{soft_n_clusters}_clusters_{n_tf}_tfs"]:
        input_target = f"input_sergio_file_targets_{n_tf}_tfs.csv"
        input_regs = f"input_sergio_file_{prefix_regs}.csv"

        ngenes = infer_ngenes(input_target, input_regs)

        with open(input_regs) as f:
            nclusters = len(f.readline().strip().split(",")) - 1
        print(prefix_regs, ngenes, nclusters, ncells_per_cluster)

        # simulate
        sim = sergio(
            number_genes=ngenes,
            number_bins=nclusters,
            number_sc=ncells_per_cluster,
            noise_params=1, decays=0.8, sampling_state=15, noise_type='dpd'
        )
        sim.build_graph(
            input_file_taregts=input_target,
            input_file_regs=input_regs,
            shared_coop_state=2
        )
        sim.simulate()
        expr = sim.getExpressions()

        # save clean raw counts
        raw_counts = sim.convert_to_UMIcounts(expr)
        raw_counts = np.concatenate(raw_counts, axis=1)
        print(raw_counts.shape)
        np.savetxt(f"sergio_raw_counts_{prefix_regs}.csv", raw_counts, delimiter=",", fmt = "%d")

        # add noise
        expr_O = sim.outlier_effect(expr, outlier_prob=0.01, mean=0.8, scale=1)

        # library size effect
        libFactor, expr_O_L = sim.lib_size_effect(expr_O, mean = 4.6, scale = 0.4)

        # Add dropouts
        binary_ind = sim.dropout_indicator(expr_O_L, shape=6.5, percentile=82)
        expr_O_L_D = np.multiply(binary_ind, expr_O_L)

        # Convert to UMI count
        count_matrix = sim.convert_to_UMIcounts(expr_O_L_D)
        count_matrix = np.concatenate(count_matrix, axis=1)
        print(count_matrix.shape)
        np.savetxt(f"sergio_noisy_counts_{prefix_regs}.csv", count_matrix, delimiter=",", fmt = "%d")
