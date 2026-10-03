# 04. Neuron gene lists

## 1. Define neuronal cells

```r
sty <- readRDS("StyPis_scRNA_annotated.rds")
neuron_cells <- colnames(sty)[sty$cell_type == "neuron"]
non_neuron_cells <- colnames(sty)[sty$cell_type != "neuron"]
DefaultAssay(sty) <- "RNA"
```

## 2. Build the neuron-expression-fraction list

Using raw counts, the fraction was the fraction of each gene’s total counts contributed by cells annotated as neuron:

```r
count_data <- GetAssayData(sty, layer = "counts")
neuron_total_counts <- Matrix::rowSums(count_data[, neuron_cells])
all_total_counts <- Matrix::rowSums(count_data)
neuron_expr_fraction <- ifelse(all_total_counts > 0,
                               neuron_total_counts / all_total_counts, NA)
```

All genes were saved to `neuron_unique_expressed_all.tsv`; the selected list required `neuron_expr_fraction >= 0.80` and was saved as `neuron_unique_expressed_frac0.80.tsv`.

## 3. Build the neuron-high-expression list

The script first calculated descriptive statistics from log-normalized data, including `log2((mean_neuron + 1e-9)/(mean_non_neuron + 1e-9))`, then used a Wilcoxon marker test:

```r
sty$temp_ident <- ifelse(sty$cell_type == "neuron", "neuron", "other")
Idents(sty) <- "temp_ident"
neuron_markers <- FindMarkers(sty, ident.1 = "neuron", ident.2 = "other",
                              only.pos = TRUE, min.pct = 0.1,
                              logfc.threshold = 0.25, test.use = "wilcox")
```

The final selected table, `neuron_highexpr_log2FC1.0.tsv`, required calculated `log2FC >= 1.0` and `p_val_adj < 0.05`; all genes were retained in `neuron_highexpr_all.tsv`. Both tables were left-joined to the supplied gene annotation table (`gene_id`, `gene_name`, `desc`).
