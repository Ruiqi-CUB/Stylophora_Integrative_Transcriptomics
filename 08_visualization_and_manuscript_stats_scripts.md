# 08. Visualization and manuscript statistics

This section summarizes extension statistics, QC, candidate-gene plots, set overlaps, and cell-fraction tests.

## 1. Extension and bulk QC summaries

`01_gene_extension_stats.py` parses `gene` rows from the original and extended GTFs, matches common gene IDs, calculates coordinate-derived length differences (`end - start + 1`), and bins absolute changes as 0, 1–100, 101–500, 501–1,000, 1,001–5,000, and >5,000 bp.

The two bulk-QC scripts parse HISAT2 summary files for input-read count and overall alignment rate. `02_bulk_rnaseq_qc_summary.py` additionally counts reads in raw and trimmed FASTQs; `02b_bulk_rnaseq_qc_fast.py` is the faster alternative and records only HISAT2 input reads and alignment rate. Neither generates a MultiQC report.

## 2. Candidate-expression plots

`01_cluster_plots.R` and `02_dotplot.R` plot a hard-coded table of nine candidate gene IDs, names, and descriptions from the annotated Seurat object. Key calls are:

```r
FeaturePlot(sty, features = gene_id, reduction = "umap", order = TRUE, pt.size = 0.5)
DotPlot(sty, features = candidate_genes$gene_id, dot.scale = 8) + coord_flip()
```

`03_boxplots.R` reads Voolstra and Savary TMM matrices and makes per-gene boxplots with raw TMM values displayed on a log10 x-axis. `04_deconvolution_boxplots.R` similarly visualizes the neuron slice of the deconvolved expression arrays. `Neuro_DEG_data2boxplots.R` constructs annotated long TMM tables by joining sample, transcript-to-protein, and emapper annotation tables; `expr_boxplot.R` accepts `--gene` or `--protein` and creates three-panel Voolstra/Savary-S/Savary-L plots.

## 3. Set intersections and summaries

`03_upset_bulk_x_neuron.R` and `05_upset_neuron_lists.R` build binary membership tables from their named gene lists and render UpSetR plots. `04a_upset_voolstra_internal.R` splits 12 Voolstra pairwise DEG tables into positive and negative log2FC sets; `04b_upset_savary_internal.R` does the same for six Savary comparisons; and `04c_upset_cross_study_extreme.R` compares the four Voolstra 36 vs 30 sets with Savary S_T1 and L_T1 34 vs 27 sets. Each produces separate up/down membership tables and plots. `05_generate_stats_summary.R` reports values saved by these three scripts.

`04_upset_final_candidates.R` compares six sets: high-expression and unique-neuron intersections plus deconvolved-neuron DEGs for ICN, and the corresponding three sets for Savary S_T1. It labels genes found in all three sets for a dataset as high-confidence and genes in all six sets as ultra-high-confidence:

```r
icn_all <- Int_highexpr_ICN == 1 & Int_unique_ICN == 1 & Deconv_ICN == 1
s_t1_all <- Int_highexpr_S_T1 == 1 & Int_unique_S_T1 == 1 & Deconv_S_T1 == 1
all_six <- icn_all & s_t1_all
```

## 4. Cell-fraction statistics

`06_cell_fraction_stats.R` parses metadata from sample names, exports fraction tables, and for every cell type applies Wilcoxon rank-sum tests (`exact = FALSE`). It tests 36°C vs 30°C separately per Voolstra population, and 34°C vs 27°C separately per Savary condition. Benjamini–Hochberg correction is applied within population/condition:

```r
wilcox.test(group_hot, group_cold, exact = FALSE)
stats_voolstra <- stats_voolstra %>% group_by(population) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))
```

## 5. Candidate heatmap

`07_candidate_heatmap_v3.R` reads candidate IDs from `all4_shared_annotated.tsv`, confirms they are present in both TMM matrices, calculates mean TMM per temperature condition, and z-scores within each gene and study:

```r
safe_zscore <- function(x) if (is.na(sd(x)) || sd(x) == 0) rep(0, length(x)) else (x - mean(x)) / sd(x)
```

The color range is clipped to z = −2 to 2. Voolstra temperatures are ordered 30, 33, 36 and populations AFR, ETR, PTR, ICN; Savary temperatures are 27, 29, 32, 34 and groups S_T1, S_T2, L_T1, L_T2. Four PDF variants are written. The script’s comment says nine candidates, but candidate membership is determined by the external input table.
