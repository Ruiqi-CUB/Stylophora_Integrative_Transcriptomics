# 01. Reference extension and scRNA reference

Key commands for 3′ extension and the combined host–symbiont Cell Ranger reference.

## 1. Build a STAR genome index

A STAR index was made from the masked host genome and original annotation:

```bash
STAR --runThreadN <threads> --runMode genomeGenerate \
  --genomeDir StyPis_STAR_index \
  --genomeFastaFiles jaStyPist1.fa.masked \
  --sjdbGTFfile braker.cellranger.gtf \
  --sjdbOverhang 60 --genomeSAindexNbases 11 \
  --limitSjdbInsertNsj 2022663
```

`<threads>` was supplied by the scheduler and is not a fixed scientific parameter.

## 2. Map pooled 10X read 2 sequences for 3′ peak discovery

Pooled 10X read 2 FASTQs were mapped:

```bash
STAR --genomeDir StyPis_STAR_index --runThreadN <threads> \
  --readFilesIn <pooled_R2_FASTQs.csv> --readFilesCommand zcat \
  --outFileNamePrefix pool_StyPis. \
  --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
  --outFilterMultimapNmax 5 --outFilterMismatchNmax 3 \
  --alignIntronMax 5500 --limitBAMsortRAM 80000000000 \
  --genomeLoad NoSharedMemory --readNameSeparator ' ' \
  --outSAMattributes Standard
```

## 3. Call strand-aware broad peaks

Reads were separated by strand and peaks called using an effective genome size of 90% of the summed chromosome lengths:

```bash
sambamba view -f bam --filter "strand=='+'" -t <threads> \
  -o pool_StyPis.plus.bam pool_StyPis.Aligned.sortedByCoord.out.bam
sambamba view -f bam --filter "strand=='-'" -t <threads> \
  -o pool_StyPis.minus.bam pool_StyPis.Aligned.sortedByCoord.out.bam

EFF=$(awk '{s+=$1} END{print int(s*0.9)}' chrLength.txt)

macs2 callpeak -t pool_StyPis.plus.bam -f BAM -g "$EFF" \
  --keep-dup 20 -q 0.01 --shift 1 --extsize 20 --broad --nomodel \
  -n pool_StyPis_plus --outdir macs2
macs2 callpeak -t pool_StyPis.minus.bam -f BAM -g "$EFF" \
  --keep-dup 20 -q 0.01 --shift 1 --extsize 20 --broad --nomodel \
  -n pool_StyPis_minus --outdir macs2
```

The plus and minus `broadPeak` outputs were combined after setting column 6 to `+` or `-`, yielding `pool_StyPis_peaks_3p.broadPeak`.

## 4. Extend gene models to 3′ peaks

Genes were extended to qualifying 3′ peaks:

```bash
Rscript qsub_extend3p.R \
  -g braker.cellranger.gtf \
  -p pool_StyPis_peaks_3p.broadPeak \
  -o extend_StyPis_genes.gtf \
  -m 5000 -a -r StyPis -q 1e-3
```

From `qsub_extend3p.R`, peaks are imported as ranges, filtered at MACS2 q-value `<= 0.001`, and assigned strand-specifically to downstream gene ends when within 5,000 bp. Overlapping end peaks are unioned with a gene; the most distant qualifying downstream orphan peak is used per gene. `-a` retains unassigned peaks as orphan genes prefixed `StyPis`. The optional 5′-overlap trimming flag was not used.

The extension output is later used as `StyPis_GeneExt.gtf`.

## 5. Build the combined Cell Ranger reference

Cell Ranger 10.0.0 built the combined reference:

```bash
cellranger mkref --genome=StyPis_SymMic_CR \
  --fasta=StyPis_SymMic.fna --genes=StyPis_GeneExt.gtf \
  --nthreads=<threads>
```

The script establishes that `StyPis_SymMic.fna` is a combined host–symbiont FASTA, but the concatenation procedure is not recoverable.

## 6. Generate the single-cell count matrix

Two runs from the same GSM library were jointly processed:

```bash
cellranger count --id=StyPis_SymMic_GSM8793941 \
  --transcriptome=StyPis_SymMic_CR \
  --fastqs=<run1_fastqs>,<run2_fastqs> \
  --expect-cells=20000 --chemistry=auto --create-bam=true \
  --localcores=<threads> --localmem=150
```

The input runs were `SRR32332770` and `SRR32332771`.
