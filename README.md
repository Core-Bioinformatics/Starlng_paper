# Starlng_paper
Contains code used to generate the figures and the analysis included in the Starlng manuscript.

The code relies on the usage of the [Starlng](https://github.com/Core-Bioinformatics/Starlng) package.

Methods that were investigated in the manuscript: Starlng, hdWGCNA, Hotspot, SCENIC, PAGA + CellRank.

Structure of the repository:
- `benchmarking`: test the runtime and memory usage of Starlng against immune-based synthetic datasets (10k-100k cells);
- `data`: preprocessing steps for the real datasets and generation of the synthetic datasets using Sergio;
- `comparison_against_ground_truth`: analysis of the similarity between the method results and the ground truth.
- `comparison_with_other_methods`: analysis of the output of each clustering method based on the 5 metrics described in the manuscript.
- `moran_i_impact`: exploratory analysis of the distribution of Moran's I values and its relationship to other cohesion metrics such as pseudotime iqr and median umap distance.
- `robustness_analysis`: analysis of the robustness of Starlng against subsampling, noise adding and seed changes.
- `starlng_run`: scripts to run Starlng on real datasets and generate the figures included in the manuscript.
