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
    soft_n_clusters = int(n_tf * 0.6)
    for prefix_regs in [f"soft_{soft_n_clusters}_clusters_{n_tf}_tfs"]:
        df = pd.read_csv(f"input_sergio_file_grn_{soft_n_clusters}_clusters_{n_tf}_tfs.tab", sep='\t', header=None, index_col=None)
        bMat = df.values
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
            dynamics = True,
            bifurcation_matrix= bMat,
            noise_params=0.2, noise_params_splice = 0.07, decays=0.8, sampling_state=1, noise_type='dpd'
        )
        sim.build_graph(
            input_file_taregts=input_target,
            input_file_regs=input_regs,
            shared_coop_state=2
        )
        sim.simulate_dynamics()
        exprU, exprS = sim.getExpressions_dynamics()

        # save clean raw counts
        raw_countsU, raw_countsS = sim.convert_to_UMIcounts_dynamics(exprU, exprS)
        raw_countsS = np.concatenate(raw_countsS, axis=1)
        print(raw_countsS.shape)
        np.savetxt(f"sergio_dynamic_ds12_raw_counts_{prefix_regs}.csv", raw_countsS, delimiter=",", fmt = "%d")


        # add noise
        exprU_O, exprS_O = sim.outlier_effect_dynamics(exprU, exprS, outlier_prob=0.01, mean=0.8, scale=1)

        # Library Size Effect
        libFactor, exprU_O_L, exprS_O_L = sim.lib_size_effect_dynamics(exprU_O, exprS_O, mean = 4.6, scale = 0.4)

        # Add dropouts
        binary_indU, binary_indS = sim.dropout_indicator_dynamics(exprU_O_L, exprS_O_L, shape = 6.5, percentile = 82)
        exprU_O_L_D = np.multiply(binary_indU, exprU_O_L)
        exprS_O_L_D = np.multiply(binary_indS, exprS_O_L)

        # Convert to UMI count
        count_matrix_U, count_matrix_S = sim.convert_to_UMIcounts_dynamics(exprU_O_L_D, exprS_O_L_D)
        count_matrix_U = np.concatenate(count_matrix_U, axis=1)
        count_matrix_S = np.concatenate(count_matrix_S, axis=1)
        print(count_matrix_U.shape, count_matrix_S.shape)
        np.savetxt(f"sergio_dynamic_s12_noisy_counts_{prefix_regs}.csv", count_matrix_S, delimiter=",", fmt = "%d")
