# 02. scRNA-seq processing

Analysis used Seurat.

## 1. Create and quality-filter the Seurat object

```r
counts <- Read10X("filtered_feature_bc_matrix")
sty <- CreateSeuratObject(counts = counts, project = "StyPis_scRNA",
                          min.cells = 3, min.features = 200)
sty <- subset(sty, subset = nFeature_RNA > 400 &
                         nFeature_RNA < 3500 & nCount_RNA < 15000)
```

No mitochondrial percentage was calculated because the genome GTF used for this analysis did not include the mitochondrial chromosome.

## 2. Normalize, identify variable features, and calculate PCA

```r
sty <- NormalizeData(sty, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = TRUE)
sty <- FindVariableFeatures(sty, selection.method = "vst",
                            nfeatures = 3000, verbose = TRUE)
sty <- ScaleData(sty, features = VariableFeatures(sty), verbose = TRUE)
sty <- RunPCA(sty, features = VariableFeatures(sty), npcs = 50, verbose = TRUE)
```

## 3. Graph clustering and UMAP

The first 20 PCs were selected after inspecting an elbow plot.

```r
pcs_to_use <- 1:20
sty <- FindNeighbors(sty, dims = pcs_to_use, verbose = TRUE)
sty <- FindClusters(sty, resolution = 0.2, verbose = TRUE)
sty <- RunUMAP(sty, dims = pcs_to_use, verbose = TRUE)
```

The source comments call this Louvain clustering, although the explicit algorithm argument is not set; the Seurat default in the installed version cannot be reconstructed.

## 4. Identify cluster markers and annotate gene names

```r
markers <- FindAllMarkers(sty, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = TRUE)
markers_annotated <- markers %>%
  left_join(gene_map, by = c("gene" = "gene_id"))
```

Unmapped names were replaced by the gene ID. The 20 largest fold-change markers per cluster were written after accommodating either `avg_log2FC` or legacy `avg_logFC` output columns.

## 5. Generate all-gene dot-plot data and save object

```r
avg_exp_mat <- AverageExpression(sty, assays = "RNA", features = rownames(sty),
                                 slot = "data", verbose = FALSE)$RNA
data_matrix <- GetAssayData(sty, slot = "data", assay = "RNA")
# for each cluster: pct_expressed = rowSums(cluster_data > 0) / ncol(cluster_data) * 100
saveRDS(sty, "StyPis_scRNA_clustered_res0.2.rds")
```

The generated long table contains `gene_id`, `gene_name`, cluster, mean normalized expression, and percent-expressing cells (`dotplot_data_full.csv`).
