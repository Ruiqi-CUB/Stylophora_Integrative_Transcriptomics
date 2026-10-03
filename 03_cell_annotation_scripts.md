# 03. Cell-type annotation

Orthology and keyword-based annotation supported the final cluster labels.

## 1. Prepare old and new genome protein sets

`deduplicate_old_genome.py` accepts `<old_proteins.faa> <old_longest_isoform.faa>` and groups RefSeq headers by `LOC` identifier (or accession if absent), retaining the longest protein per gene.

`extract_proteins.py` accepts `<new_transcripts.fa> <new_proteins.faa> [min_protein_len]`. It keeps the longest transcript per `g<digits>` gene, searches all six frames for start-to-stop ORFs, retains the longest ORF, and defaults to a minimum length of 30 amino acids.

## 2. Infer orthogroups

From `run_orthofinder.sh`:

```bash
orthofinder -f proteins/ -t 16 -a 16 -o results/ -n StyPis_old_vs_new
```

The OrthoFinder version was printed at runtime but its output is unavailable.

## 3. Transfer top old-genome markers through orthogroups

`step1_filter_markers.R` selected old gene IDs matching `^Spis_(XP_|YP_|LOC)`, ranked each old cell type by fold change, and retained the top 25 per type:

```r
fc_filtered <- fc_data[grep("^Spis_(XP_|YP_|LOC)", rownames(fc_data), value = TRUE), ]
top_genes <- head(names(sort(fc_filtered[, cell_type], decreasing = TRUE)), 25)
```

It stripped `Spis_`, joined `Spis_gene_combined.tsv`, and wrote `top_markers_per_celltype.tsv`.

`step2_find_orthologs.R` read `Orthogroups.tsv`, converted an isoform suffix such as `_1` to `.1` when searching old gene IDs, then collected all comma-separated IDs from the `StyPis_new_genome` column for each old-genome marker. The result was `top_markers_with_orthologs.tsv`.

## 4. Independent annotation support from keyword searches

`bottom_up_analysis.R` read cell-type keywords from `cell_keywords.csv` and searched case-insensitively in `gene_name` and `desc`. Keywords of five or fewer characters were whole-word searches; longer terms were substring searches:

```r
pattern <- if (nchar(keyword) <= 5) paste0("\\b", keyword, "\\b") else keyword
matches <- bind_rows(
  gene_annot[grepl(pattern, gene_annot$gene_name, ignore.case=TRUE, perl=TRUE), ],
  gene_annot[grepl(pattern, gene_annot$desc, ignore.case=TRUE, perl=TRUE), ]
) %>% distinct()
```

Only genes present in `StyPis_scRNA_clustered_res0.2.rds` were retained. Per cluster and cell type, the score was `sum(avg_expr * pct_expressed)` over matched genes. The related `create_annotations.R` instead used the mean of `avg_expr * pct_expressed` and classified confidence as high when top score >5 and >1.5× the runner-up, medium when >1, otherwise low. These are alternative/supporting procedures; the scripts do not establish which produced the final manual mapping.

## 5. Inspect supporting expression plots

`step4_create_plots.R` plotted ortholog-derived lists against the clustered Seurat object using `DotPlot`, `FeaturePlot`, and `VlnPlot`; it accepts either `StyPis_scRNA_clustered_res0.2.rds` or a fallback `res0.4` object. It makes pages of nine UMAP panels (3×3) and uses old-genome fold-change values for comparative dot plots. Plotting parameters are descriptive rather than an additional inferential analysis.

## 6. Apply the final cluster-to-cell-type mapping

The final operation in `01_add_cell_annotations.R` is explicit, but its mapping file is a required external/manual input:

```r
cell_annot <- read.csv("cell_annotation.csv") # required columns: cluster, type
cluster_to_type <- setNames(cell_annot$type, as.character(cell_annot$cluster))
sty$cell_type <- cluster_to_type[as.character(sty$seurat_clusters)]
saveRDS(sty, "StyPis_scRNA_annotated.rds")
```

Consequential ambiguity: the contents and derivation of `cell_annotation.csv` are not included, so final labels cannot be reproduced from this bundle alone.
