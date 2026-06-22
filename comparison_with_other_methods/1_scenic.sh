#!/bin/bash

data_folder="/mnt/d/Starlng_paper/data"
output_folder="/mnt/d/Starlng_paper/comparison_with_other_methods/scenic_results"
if [ ! -d "$output_folder" ]; then
    mkdir -p "$output_folder"
fi

dts_name="cao_Liver_subset"
echo "Processing ${dts_name}"
input_file="${data_folder}/${dts_name}.loom"
if [ ! -f "$input_file" ]; then
    echo "Input file $input_file does not exist. Skipping."
    continue
fi

echo "Step 1"
pyscenic grn \
    --num_workers 10 \
    -o "${output_folder}/${dts_name}_expr_mat.adjacencies.tsv" \
    "$input_file" \
    "${output_folder}/allTFs_hg38.txt"

echo "Step 2"

pyscenic ctx \
    "${output_folder}/${dts_name}_expr_mat.adjacencies.tsv" \
    "${output_folder}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
    --annotations_fname "${output_folder}/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
    --expression_mtx_fname "$input_file" \
    --output "${output_folder}/${dts_name}_reg.csv" \
    --mask_dropouts \
    --num_workers 10

echo "Step 3"
pyscenic aucell \
    "$input_file" \
    "${output_folder}/${dts_name}_reg.csv" \
    --output "${output_folder}/${dts_name}_pyscenic_output.loom" \
    --num_workers 10

