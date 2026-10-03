# 06. Bulk DEG lists and neuron intersections

Bulk heat-response sets were intersected with the two neuronal gene sets.

## 1. Summarize the two scRNA neuron lists

`01_generate_scRNA_distributions.R` bins `neuron_expr_fraction` from 0 to 1 in increments of 0.1, and bins high-expression `log2FC` as `<0`, 0–1, …, 7–8, and ≥8. It counts genes meeting `p_val_adj < 0.05` in each log2FC bin and writes `scRNA_genelist_distributions.txt`.

## 2. Construct Voolstra 2020 heat-response lists

Input files are the four prefiltered 36°C vs 30°C lists: `ICN`, `AFR`, `ETR`, and `PTR`. Direction is `up` when `log2FC > 0`, otherwise `down`.

```r
nonpreloaded <- all_degs %>%
  group_by(transcript_id, direction) %>%
  summarize(n_populations = n(), populations_list = paste(population, collapse=","),
            mean_log2FC = mean(log2FC), .groups="drop") %>%
  filter(n_populations >= 3)

preloaded_gene_ids <- setdiff(ICN_DEGs, unique(c(AFR_DEGs, ETR_DEGs, PTR_DEGs)))
```

Outputs are `voolstra2020_heat_nonpreloaded.tsv` (same direction in ≥3/4 populations), `voolstra2020_heat_preloaded.tsv` (ICN-only DEGs), and `voolstra2020_heat_ICN36v30.tsv` (all ICN 36 vs 30 DEGs).

## 3. Construct Savary 2021 heat-response lists

`03_build_savary2021_DEG_lists.R` reads already-filtered S and L T1 comparisons, adds direction by the sign of `log2FC`, sorts by absolute fold change, and writes:

- `savary2021_heat_S_T1_34v27.tsv` from `S_T1_34_vs_S_T1_27`;
- `savary2021_heat_L_T1_34v27.tsv` from `L_T1_34_vs_L_T1_27`.

## 4. Conserved cross-study list: not recoverable

`05_intersect_neuron_DEG.R` expects `heat_conserved_ICN_S_T1_L_T1.tsv`, but the stated producing script, `04_build_conserved_DEG_list.R`, is zero bytes. Its selection logic and thresholds cannot be documented.

## 5. Intersect neuron and bulk-gene sets

The script uses `neuron_unique_expressed_frac0.80.tsv` and `neuron_highexpr_log2FC1.0.tsv`, then joins each with six bulk DEG lists on `gene_id = transcript_id`:

```r
intersection <- neuron_data %>%
  inner_join(deg_data, by = c("gene_id" = "transcript_id"))
```

This produces 12 result tables (2 neuron lists × 6 bulk lists) and counts up/down genes from the bulk `direction` field.

## 6. Add gene annotation columns

`06_add_annotations.R` left-joins the supplied `StyPis_geneID2gene.tsv` annotation (`gene_id`, `gene_name`, `desc`) onto all DEG and intersection tables, replacing any pre-existing `gene_name`/`desc` fields.
