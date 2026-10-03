# Methods command guide

## Manuscript

**Integrating single-cell and bulk transcriptomes identifies neuronal heat-stress-associated candidate genes in a reef-building coral.**

This folder is a publication-safe, method-focused guide to the original analysis scripts. It records the substantive commands, function calls, settings, inputs, and outputs that could be verified in those scripts, while omitting local paths, scheduler configuration, environment activation, and routine file management.

These documents are **not complete, directly executable pipelines**. They are intentionally shortened records of the key analytical commands and parameters.

## Workflow overview

1. [Reference extension and scRNA reference](01_reference_extension_scripts.md)
2. [scRNA-seq processing](02_scrna_processing_scripts.md)
3. [Cell-type annotation](03_cell_annotation_scripts.md)
4. [Neuron gene lists](04_neuron_gene_lists_scripts.md)
5. [Bulk RNA-seq processing](05_bulk_rnaseq_processing_scripts.md)
6. [Bulk DEG lists and neuron intersections](06_deg_intersection_scripts.md)
7. [Deconvolution and neuron-specific analysis](07_deconvolution_scripts.md)
8. [Visualizations and manuscript statistics](08_visualization_and_manuscript_stats_scripts.md)

## Verified inputs and outputs

| Stage | Main verified inputs | Main verified outputs |
|---|---|---|
| 01 | Masked host genome FASTA, `braker.cellranger.gtf`, pooled 10X R2 reads, 3′ peaks | `StyPis_GeneExt.gtf`; combined host–symbiont Cell Ranger reference; 10X filtered feature-barcode matrix |
| 02 | Cell Ranger `filtered_feature_bc_matrix`; gene-ID mapping | `StyPis_scRNA_clustered_res0.2.rds`, cluster-marker tables, full dot-plot data |
| 03 | Clustered Seurat object, old-genome marker fold-change table, annotations, old/new protein sets | Ortholog/keyword-supported annotation materials and `StyPis_scRNA_annotated.rds` |
| 04 | Annotated Seurat object and gene annotation table | Neuron-enriched and neuron-expression-fraction gene tables |
| 05 | Paired bulk reads, combined HISAT2 index, extended GTF, sample/contrast files | Voolstra and Savary gene-count matrices, TMM matrices, DESeq2 result-derived DEG lists |
| 06 | Bulk DEG lists and two filtered neuron-gene tables | Study-specific/conserved DEG lists and 12 neuron × bulk intersections |
| 07 | Annotated scRNA object; Voolstra/Savary count matrices; deconvolution inputs | Cell-type fractions, deconvolved expression arrays, neuron-specific DEG results, method comparisons |
| 08 | GTFs, bulk QC summaries, candidate list, TMM matrices, Seurat/deconvolution objects, set tables | Extension/QC summaries, plots, UpSet membership tables, cell-fraction statistics, candidate heatmaps |

## Missing information or dependencies

- No combined host–symbiont HISAT2 index-build or FASTA/GTF concatenation script was recovered (`SCRIPT_BUNDLE_MANIFEST.md`).
- No MultiQC wrapper, final BLASTp-to-nr annotation script, or script that creates the Savary S/L subset matrices was recovered (`SCRIPT_BUNDLE_MANIFEST.md`; `Savary2021_gene_DEG.slurm`).
- `04_build_conserved_DEG_list.R`, `step3_build_gene_lists.R`, and `Voolstra2021_HISAT2.slurm` are empty in this bundle, so their implementation cannot be documented.
- Software versions are reported only where explicit: Cell Ranger 10.0.0 and HISAT2 2.2.0. Other environments/package versions were not captured reliably; several R scripts print `sessionInfo()` at runtime but the outputs are unavailable.
- Sample-sheet contents, candidate-list construction, and manual cluster-to-cell-type assignments are inputs rather than recoverable computational steps.
