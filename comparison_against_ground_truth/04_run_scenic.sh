#!/bin/bash

data_folder="/mnt/d/Starlng_paper/data"
output_folder="/mnt/d/Starlng_paper/comparison_against_ground_truth/scenic"
if [ ! -d "$output_folder" ]; then
    mkdir -p "$output_folder"
fi

n_tfs=(5 10 15 20)
count_types=("raw" "noisy")
generation_type=("" "dynamic_ds12_")

for n_tf in "${n_tfs[@]}"; do
    for count_type in "${count_types[@]}"; do
        for gen_type in "${generation_type[@]}"; do
            echo "Processing n_tf=${n_tf}, count_type=${count_type}, generation_type=${gen_type}"
            input_file="${data_folder}/sergio_${gen_type}${count_type}_${n_tf}_tfs_seurat.loom"
            if [ ! -f "$input_file" ]; then
                echo "Input file $input_file does not exist. Skipping."
                continue
            fi

            echo "Step 1"
            pyscenic grn \
                --num_workers 10 \
                -o "${output_folder}/sergio_${gen_type}${count_type}_${n_tf}_tfs_expr_mat.adjacencies.tsv" \
                "$input_file" \
                "${output_folder}/allTFs_hg38.txt"
            
            echo "Step 2"

            pyscenic ctx \
                "${output_folder}/sergio_${gen_type}${count_type}_${n_tf}_tfs_expr_mat.adjacencies.tsv" \
                "${output_folder}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
                --annotations_fname "${output_folder}/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
                --expression_mtx_fname "$input_file" \
                --output "${output_folder}/sergio_${gen_type}${count_type}_${n_tf}_tfs_reg.csv" \
                --mask_dropouts \
                --num_workers 10
            
            echo "Step 3"
            pyscenic aucell \
                "$input_file" \
                "${output_folder}/sergio_${gen_type}${count_type}_${n_tf}_tfs_reg.csv" \
                --output "${output_folder}/sergio_${gen_type}${count_type}_${n_tf}_tfs_pyscenic_output.loom" \
                --num_workers 10
        done
    done
done
