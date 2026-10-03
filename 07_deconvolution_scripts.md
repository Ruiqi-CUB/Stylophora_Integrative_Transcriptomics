# 07. Deconvolution and neuron-specific analysis

Single-cell reference profiles were used to deconvolve bulk expression and test neuronal responses.

## 1. Prepare scRNA and bulk inputs

`03_prep_instaprism.R` extracts the RNA assay’s log-normalized `data` layer and cell-type labels from the annotated Seurat object:

```r
expr_matrix <- GetAssayData(sty, slot = "data", assay = "RNA")
cell_types <- sty$cell_type
cell_types[is.na(cell_types)] <- "unknown"
saveRDS(expr_matrix, "scRNA_expr_matrix.rds")
saveRDS(cell_types, "scRNA_cell_types.rds")
```

`05_prep_bulk_input.R` converts each gene-count matrix to CPM, retains only genes shared with the scRNA reference, and writes RDS bulk matrices:

```r
normalize_to_cpm <- function(counts_matrix) sweep(counts_matrix, 2, colSums(counts_matrix), "/") * 1e6
voolstra_filtered <- normalize_to_cpm(voolstra_bulk)[intersect(rownames(voolstra_bulk), scRNA_genes), ]
savary_filtered <- normalize_to_cpm(savary_bulk)[intersect(rownames(savary_bulk), scRNA_genes), ]
```

It also writes equivalent CSVs for EPIC-unmix, but no EPIC-unmix analysis script is included.

## 2. Run InstaPrism

`01_InstaPrism.R` builds the reference for Voolstra; `01b_Savary2021_InstaPrism.R` reuses it for Savary:

```r
refPhi_obj <- refPrepare(sc_Expr = as.matrix(scRNA_expr),
                         cell.type.labels = scRNA_cell_types,
                         cell.state.labels = scRNA_cell_types,
                         pseudo.min = 1e-8)
fit <- InstaPrism(bulk_Expr = as.matrix(bulk_expr), refPhi_cs = refPhi_obj,
                  filter = TRUE, outlier.cut = 0.01, outlier.fraction = 0.1,
                  pseudo.min = 1e-8, verbose = TRUE,
                  convergence.plot = FALSE, n.core = 1)
Z_array <- get_Z_array(fit)
```

The samples × cell-types fraction matrix was obtained by transposing the fitted object’s `theta` cell-type-fraction slot; `Z_array` is samples × genes × cell types. InstaPrism package version is not recorded.

## 3. Test cell-fraction shifts and neuron-specific expression

`02_ICN_36v30_neuron_analysis.R` compares ICN 36°C vs 30°C samples (source comments report n=7 and n=6); `04_3pop_36v30_neuron_analysis.R` applies the same function to AFR, ETR, and PTR. `02b_Savary2021_neuron_analysis.R` analogously analyzes Savary comparisons selected by sample-name patterns. For each cell type, the scripts calculate mean ratio and use a two-sided Wilcoxon rank-sum test with `exact = FALSE`, then Benjamini–Hochberg adjust across cell types.

Neuron expression is selected from `Z_array[, , "neuron"]`, transposed to genes × samples, and filtered:

```r
neuron_expr_filtered <- neuron_expr_t[rowSums(neuron_expr_t) > 1, ]
```

For each retained gene, 36°C versus 30°C (or 34°C versus 27°C) uses:

```r
log2FC <- log2((mean(expr_hot) + 1e-6) / (mean(expr_cold) + 1e-6))
pval <- suppressWarnings(wilcox.test(expr_hot, expr_cold, exact = FALSE)$p.value)
padj <- p.adjust(pval, method = "BH")
```

Neuron DEGs require `abs(log2FC) > 1` and `padj < 0.05`. These scripts explicitly treat the input as normalized, continuous CPM-like deconvolved expression—not raw counts—and therefore do not use DESeq2/edgeR for this stage.

## 4. Compare intersection and deconvolution methods

`05_method_comparison.R` compares each intersection set with a corresponding significant deconvolved-neuron DEG table. It uses gene-ID overlap, compares directions (`direction` versus sign of deconvolution `log2FC`), and writes shared genes with concordant direction. The script specifies four comparisons: ICN and Savary S_T1, each against both neuron-list definitions.
