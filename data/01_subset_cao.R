library(Seurat)

if (basename(getwd()) != "data") {
    setwd("data")
}

so <- qs2::qs_read("cao_Liver_filtered_normalized.qs2", nthreads = 30)
umap_emb <- Embeddings(so, reduction = "umap")
mask <- umap_emb[,1] > -1
cells_emph <- rownames(umap_emb)[mask]

DimPlot(so, reduction = "umap", cells.highlight = cells_emph)
so <- subset(so, cells = cells_emph)

so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
so <- ScaleData(so, features = rownames(so), verbose = FALSE)
so <- RunPCA(so)
set.seed(42)
so <- RunUMAP(so, dims = 1:30, reduction = "pca")
so <- FindNeighbors(so, dims = 1:30, k.param = 30, reduction = "pca")
so <- FindClusters(so, random.seed = 42)

qs2::qs_save(so, "cao_Liver_subset_filtered_normalized.qs2", nthreads = 30)

DimPlot(so)
