# 05. Bulk RNA-seq processing

Sources: `Voolstra2020_Trim.slurm`, `Voolstra2021_FastQC.slurm`, `Voolstra2020_StyPis_featureCounts_gene.slurm`, `Voolstra2020_gene_DEG.slurm`, `Savary_Trim.slurm`, `Savary2021_FastQC.slurm`, `Savary2021_HISAT2_New.slurm`, `Savary2021_featureCounts_gene.slurm`, and `Savary2021_gene_DEG.slurm`. `Voolstra2021_HISAT2.slurm` is empty.

## 1. Trim paired-end reads and assess QC

The Voolstra (83-array task) and Savary (80-array task) scripts used the same Trimmomatic settings, with input sample names read from their respective sample lists:

```bash
trimmomatic PE -threads <threads> -phred33 \
  <sample>_1.fastq.gz <sample>_2.fastq.gz \
  <sample>_1P.qtrim.fq.gz <sample>_1U.qtrim.fq.gz \
  <sample>_2P.qtrim.fq.gz <sample>_2U.qtrim.fq.gz \
  ILLUMINACLIP:<TruSeq3-PE.fa>:2:30:10 \
  SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25

fastqc -t <threads> -o FastQC_raw <sample>_1.fastq.gz <sample>_2.fastq.gz
fastqc -t <threads> -o FastQC_trimmed <sample>_1P.qtrim.fq.gz <sample>_2P.qtrim.fq.gz
```

## 2. Map Savary reads to the combined reference

`Savary2021_HISAT2_New.slurm` specifies HISAT2 2.2.0 and maps each trimmed pair as follows:

```bash
hisat2 -p <threads> --dta -k 1 \
  -x StyPis_SymMic_hisat2_index \
  -1 <sample>_1P.qtrim.fq.gz -2 <sample>_2P.qtrim.fq.gz \
  -S <sample>.sam --summary-file <sample>.hisat2.summary.txt
samtools sort -@ <threads> -o <sample>.bam <sample>.sam
samtools index <sample>.bam
```

The index construction and host/symbiont-reference assembly are not recoverable. The Voolstra mapping source file is empty, so its mapping command/options are also not recoverable.

## 3. Generate gene-level count matrices

Both study-specific featureCounts scripts use the extended GTF and identical settings:

```bash
featureCounts -a StyPis_GeneExt.gtf -t exon -g gene_id \
  -s 2 -p -B -T <threads> \
  -o <study>_StyPis.gene.counts.txt *.bam
```

`-s 2` denotes reverse-strand counting. The scripts remove comment lines, retain the `gene_id` column plus columns 7 onward, rename `Geneid` to `gene_id`, and strip `.bam` from sample headers to produce `<study>_StyPis.gene.counts.matrix.tsv`.

## 4. Bulk QC, normalization, and differential expression

For Voolstra, `Voolstra2020_gene_DEG.slurm` ran the following on the count matrix and externally supplied sample/contrast files:

```bash
PtR --matrix Voolstra2020_StyPis.gene.counts.matrix.tsv \
  --min_rowSums 10 -s <Voolstra_sample_file> --log2 --CPM \
  --sample_cor_matrix --output Voolstra2020_gene
PtR --matrix Voolstra2020_StyPis.gene.counts.matrix.tsv \
  --min_rowSums 10 -s <Voolstra_sample_file> --log2 --CPM \
  --center_rows --prin_comp 3 --output Voolstra2020_gene
normalize_counts.R -in Voolstra2020_StyPis.gene.counts.matrix.tsv \
  -out Voolstra2020_gene_TMM.EXPR.matrix
run_DE_analysis.pl --matrix Voolstra2020_StyPis.gene.counts.matrix.tsv \
  --method DESeq2 --samples_file <Voolstra_sample_file> \
  --contrasts <Voolstra_contrasts_file> --output DESeq2_Voolstra2020_gene
```

For Savary, the same PtR/TMM/DESeq2 sequence was applied to a full matrix for global QC and separately to `Savary2021_S_gene` and `Savary2021_L_gene` subset matrices for DESeq2. The procedure that created these S/L matrices is absent.

For every `*.DE_results`, both scripts selected genes with absolute log2 fold change >1 and adjusted p-value <0.05 (columns 7 and 11 respectively):

```bash
awk -F'\t' 'NR>1 && $1!="" && (($7 > 1) || ($7 < -1)) && ($11 < 0.05) {print $1, $7, $11}' \
  <comparison>.DE_results > <comparison>_DE_FC1.txt
```

They then formed a unique-gene union across contrast lists. Exact sample assignments and contrasts are defined only in missing external files.
