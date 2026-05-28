# rnaseq-flow — Output Guide

How to find and interpret every file the pipeline produces. All paths are
relative to `--outdir` (default `results/`).

See [README.md](README.md) for an overview and [USAGE.md](USAGE.md) for run
instructions.

---

## Output directory layout

```
results/
├── fastqc/             Raw-read QC (FastQC)
├── fastp/              Trimmed reads + trimming reports
├── star/   | hisat2/   Aligned BAMs, indices, alignment logs
├── salmon/ | kallisto/ Transcript-level quantification
├── rseqc/              Strandedness, read distribution, gene-body coverage
├── bam_to_bigwig/      CPM-normalised coverage tracks (.bw)
├── featurecounts/      Per-sample gene count tables (STAR/HISAT2)
├── featurecounts_exon/ Per-exon count tables for diffSplice (STAR/HISAT2)
├── tximport/           Salmon/Kallisto counts summarised to gene level
├── deseq2_output/      Differential expression — DESeq2
├── edger_output/       Differential expression — edgeR
├── gsea_output/        Gene-set enrichment (fgsea)
├── gprofiler_output/   GO / pathway over-representation (gprofiler2)
├── rmats_output/       Alternative splicing events
├── star_fusion/        Gene-fusion predictions
├── isoform_switch/     Isoform-switching analysis
├── dtu_output/         Differential transcript usage (DEXSeq)
├── diffsplice_output/  Differential splicing (edgeR diffSpliceDGE)
├── multiqc/            Aggregated MultiQC report
├── quarto_report/      Interactive Quarto analysis report
└── pipeline_info/      Nextflow execution timeline / trace / report
```

Which directories appear depends on the aligner and options you chose.

---

## Recommended reading order

1. **`multiqc/multiqc_report.html`** — start here; one page for whole-run QC.
2. **Alignment logs** — confirm mapping rates are acceptable.
3. **`rseqc/`** — confirm strandedness and library quality.
4. **`deseq2_output/pca_plot.png`** — do samples group by condition?
5. **DE tables** — the differentially expressed genes.
6. **Enrichment / splicing / fusions** — biological interpretation.

---

## 1. Read QC — `fastqc/`, `fastp/`

**`fastqc/`** — `*_fastqc.html` per sample. Check per-base quality (should stay
in the green), adapter content, and duplication. Raw RNA-seq normally shows
some duplication and a biased first ~12 bp — that is expected.

**`fastp/`** — trimmed FASTQ files plus `*.fastp.html` / `*.fastp.json`. The
report shows reads before/after filtering and adapter removal. A large drop in
read count means aggressive filtering — inspect the input quality.

> All QC metrics are also aggregated in `multiqc/` — usually easier to compare
> across samples there.

---

## 2. Alignment — `star/` or `hisat2/`

| File | Description |
|---|---|
| `*.Aligned.sortedByCoord.out.bam` / `*.bam` | Coordinate-sorted alignments |
| `*.bai` | BAM index (for IGV and downstream tools) |
| `*.Log.final.out` (STAR) | Alignment summary statistics |
| `*.hisat2.summary.log` (HISAT2) | Alignment summary statistics |
| `*.ReadsPerGene.out.tab` (STAR) | STAR's own gene counts (not used downstream) |

**Reading the STAR `Log.final.out`** — the key line is *Uniquely mapped reads
%*:

| Uniquely mapped % | Interpretation |
|---|---|
| > 85% | Excellent |
| 70–85% | Acceptable |
| < 70% | Investigate — contamination, wrong genome, or poor quality |

Also watch *% of reads mapped to multiple loci* (high values suggest rRNA or
repetitive contamination) and *% of reads unmapped: too short* (often adapter
or quality problems).

---

## 3. Transcript quantification — `salmon/`, `kallisto/`

One sub-directory per sample.

**Salmon** — `<sample>/quant.sf`, columns:

| Column | Meaning |
|---|---|
| `Name` | Transcript ID |
| `Length` / `EffectiveLength` | Transcript length / length corrected for fragment bias |
| `TPM` | Transcripts Per Million — normalised abundance, comparable across samples |
| `NumReads` | Estimated reads assigned to the transcript |

**Kallisto** — `<sample>/abundance.tsv` with `target_id`, `length`,
`eff_length`, `est_counts`, `tpm`.

Use `TPM` to compare expression directly. The estimated counts
(`NumReads`/`est_counts`) are summarised to gene level by the pipeline's
tximport step (see §6) and fed straight into DESeq2/edgeR.

---

## 4. Alignment QC — `rseqc/`

| File | What it tells you |
|---|---|
| `*.infer_experiment.txt` | **Library strandedness**, measured empirically |
| `*.strandedness.txt` | The one-word verdict (`forward`/`reverse`/`unstranded`) derived from the above |
| `*.read_distribution.txt` | Fraction of reads in CDS / UTR / intron / intergenic |
| `*.geneBodyCoverage.txt` + `.curves.pdf` | 5'–3' coverage evenness |

**`infer_experiment.txt`** reports two fractions. Roughly:

- Both near **0.5** → `unstranded`
- `"++,--"` fraction near **1.0** → `forward`
- `"+-,-+"` fraction near **1.0** → `reverse`

**Automatic feedback.** With `--strandedness auto` (the default), this
inference is not just reported — `*.strandedness.txt` records the per-sample
verdict and it is fed directly into that sample's featureCounts run, so the
correct `-s` flag is used without you having to know the library type. Setting
`--strandedness` explicitly overrides the inference; if you do, check it
against `*.strandedness.txt`, since a wrong strand setting roughly halves your
gene counts.

**`read_distribution.txt`** — most tags should fall in exonic regions (CDS +
UTRs). A high intronic/intergenic fraction suggests DNA contamination or
incomplete annotation.

**Gene-body coverage** — should be a flat plateau. A strong 3' skew indicates
RNA degradation; a 5' skew can indicate library-prep bias.

---

## 5. Coverage tracks — `bam_to_bigwig/`

`<sample>.bw` — CPM-normalised coverage in BigWig format. Load into IGV or the
UCSC Genome Browser to inspect coverage at specific loci. CPM normalisation
makes tracks comparable across samples of different depth.

---

## 6. Gene counts — `featurecounts/` and `tximport/`

Gene-level counts are the input to DESeq2 and edgeR. How they are produced
depends on the aligner.

**`featurecounts/` (STAR / HISAT2).** `<sample>.featureCounts.txt` holds
gene-level counts — columns `Geneid`, `Chr`, `Start`, `End`, `Strand`,
`Length`, and a final column of raw read counts.
`<sample>.featureCounts.txt.summary` breaks down assigned vs unassigned reads
(and *why* reads were unassigned).

**`tximport/` (Salmon / Kallisto).** `tximport_gene_counts.csv` is the
gene-by-sample count matrix obtained by summarising transcript-level
quantification to gene level (genes × samples). `txi.rds` is the full tximport
object — counts, abundances and a transcript-length matrix — which DESeq2
imports with `DESeqDataSetFromTximport` and edgeR imports as a length offset,
giving more accurate normalisation than counts alone.

Either way these counts are **not** normalised for library size — do not
compare them directly between samples; that is what DESeq2/edgeR do.

---

## 7. Differential expression — `deseq2_output/`

DESeq2 and edgeR (§8) both run, **in parallel**, on the *same* gene counts
(featureCounts for STAR/HISAT2, tximport for Salmon/Kallisto) as two independent
callers — neither is downstream of the other. This section covers DESeq2.

One results CSV per pairwise contrast, named
`deseq2_results_<A>_vs_<B>.csv` (fold change is A relative to B).

| Column | Meaning |
|---|---|
| *(row name)* | Gene ID |
| `gene_name` | Gene symbol, parsed from the GTF |
| `gene_biotype` | Gene biotype (e.g. protein_coding, lncRNA) |
| `baseMean` | Mean normalised count across all samples |
| `log2FoldChange` | Effect size — log2(A / B), **apeglm-shrunken** (see note) |
| `lfcSE` | Posterior SD of the shrunken log2 fold change (apeglm) |
| `stat` | Wald test statistic (from the unshrunken fit) |
| `pvalue` | Raw p-value |
| `padj` | **Benjamini–Hochberg adjusted p-value (FDR)** — use this |

**LFC shrinkage.** The `log2FoldChange` is shrunk with the apeglm estimator
(DESeq2 `lfcShrink`): low-count, high-variance genes are pulled toward zero,
giving more reliable effect sizes and cleaner MA / volcano plots. Shrinkage
changes only the effect-size estimate — the `stat`, `pvalue` and `padj` come
from the unshrunken Wald test, so significance calling is unaffected.

**Calling a gene significant** — a common threshold is `padj < 0.05` and
`|log2FoldChange| > 1` (a 2-fold change). `padj` may be `NA` when DESeq2
filters a gene for low counts or flags it as an outlier — that is normal.

Plots:

| File | Use |
|---|---|
| `pca_plot.png` | Sample clustering — replicates of a condition should group together; if not, suspect batch effects or mislabelling |
| `maplot_<A>_vs_<B>.png` | Fold change vs mean expression; significant genes highlighted |
| `volcano_<A>_vs_<B>.png` | Volcano plot — log2 fold change vs −log10 p-value; genes with `padj < 0.05` and `|log2FC| > 1` are coloured (up red, down blue) |
| `heatmap_top_var.png` | The 20 most variable genes across samples |

---

## 8. Differential expression — `edger_output/`

edgeR runs in parallel with DESeq2 on the same gene-count input — neither is
downstream of the other. It is an independent cross-check; its results are
**not** consumed by any later step (GSEA and gProfiler use the DESeq2 results
only). Results CSVs are `edger_results_<A>_vs_<B>.csv`.

| Column | Meaning |
|---|---|
| *(row name)* | Gene ID |
| `gene_name` | Gene symbol, parsed from the GTF |
| `gene_biotype` | Gene biotype (e.g. protein_coding, lncRNA) |
| `logFC` | log2 fold change |
| `logCPM` | Average log2 counts per million |
| `F` or `LR` | Test statistic (QL F-test or likelihood-ratio) |
| `PValue` | Raw p-value |
| `FDR` | **Adjusted p-value** — use this |

`mds_plot.png` is edgeR's sample-similarity plot (analogous to the PCA);
`smear_<A>_vs_<B>.png` plots logFC vs logCPM with significant genes marked;
`volcano_<A>_vs_<B>.png` plots logFC vs −log10 p-value, colouring genes with
`FDR < 0.05` and `|logFC| > 1` (up red, down blue).

> Genes called significant by **both** DESeq2 and edgeR are the most robust
> hits. Modest disagreement near the significance threshold is expected.

---

## 9. Gene-set enrichment — `gsea_output/`

fgsea is run on the **ranked** DESeq2 gene list for each contrast (genes ranked
by the `stat` column), so it captures coordinated, subtle shifts a hard cutoff
would miss.

`gsea_stats_<contrast>.csv` columns:

| Column | Meaning |
|---|---|
| `pathway` | Gene-set name |
| `pval` / `padj` | Enrichment p-value / FDR |
| `ES` | Enrichment score |
| `NES` | **Normalised enrichment score** — sign gives direction |
| `size` | Number of genes from the set found in the data |

Interpretation: `padj < 0.05` is significant; positive `NES` = the pathway is
up in condition A, negative `NES` = up in B. `gsea_plot_<contrast>.png` shows
the top up/down pathways.

---

## 10. GO / pathway over-representation — `gprofiler_output/`

gprofiler2 tests the **significant gene lists from the DESeq2 results** for
enriched GO terms and pathways, separately for up- and down-regulated genes.

| File | Contents |
|---|---|
| `gprofiler_UP_<contrast>.csv` | Terms enriched among up-regulated genes |
| `gprofiler_DOWN_<contrast>.csv` | Terms enriched among down-regulated genes |
| `gostplot_UP/DOWN_<contrast>.png` | Manhattan-style enrichment plot |

Key columns: `term_name`, `source` (GO:BP/MF/CC, KEGG, REACTOME…),
`p_value` (already multiple-testing corrected by g:SCS), `intersection_size`
(your genes in the term).

> GSEA (§9) uses the *whole ranked list*; gProfiler uses a *thresholded list*.
> They answer slightly different questions — agreement between them strengthens
> a conclusion.

---

## 11. Alternative splicing — `rmats_output/`

One sub-directory per condition pair. rMATS reports five event types:

| Code | Event |
|---|---|
| `SE` | Skipped exon |
| `MXE` | Mutually exclusive exons |
| `A3SS` | Alternative 3' splice site |
| `A5SS` | Alternative 5' splice site |
| `RI` | Retained intron |

For each type there are two files:

- `*.MATS.JC.txt` — junction counts only (reads spanning the splice junction).
- `*.MATS.JCEC.txt` — junction counts **plus** reads on the exon body.

Key columns:

| Column | Meaning |
|---|---|
| `GeneID`, `geneSymbol` | Gene |
| `IncLevel1`, `IncLevel2` | Inclusion levels (PSI) in condition 1 and 2 |
| `IncLevelDifference` | PSI difference (1 − 2) — the effect size |
| `PValue`, `FDR` | Significance of the difference |

A common cutoff: `FDR < 0.05` and `|IncLevelDifference| > 0.1`.

---

## 12. Gene fusions — `star_fusion/`

`<sample>.star-fusion.fusion_predictions.tsv` (full) and
`.abridged.tsv` (summary). Key columns:

| Column | Meaning |
|---|---|
| `#FusionName` | The two partner genes, e.g. `GENE1--GENE2` |
| `JunctionReadCount` | Reads spanning the fusion breakpoint |
| `SpanningFragCount` | Read pairs flanking the breakpoint |
| `LeftBreakpoint` / `RightBreakpoint` | Genomic coordinates |
| `LargeAnchorSupport` | Whether reads have long anchors (more reliable) |
| `FFPM` | Fusion fragments per million — normalised support |

Higher junction + spanning support and `YES` large-anchor support indicate
more confident calls. Always validate fusion candidates against known biology
and, ideally, orthogonal evidence.

---

## 13. Isoform switching — `isoform_switch/`

| File | Contents |
|---|---|
| `isoform_switches.csv` | Top genes with significant isoform switches |
| `switch_plot_<gene>.png` | Per-gene isoform usage across conditions |
| `switchList_analyzed.rds` | The full R object for custom downstream analysis |

A switch is a gene where the *dominant* transcript isoform changes between
conditions even if total gene expression does not — biologically important and
invisible to gene-level DE. Sort `isoform_switches.csv` by q-value.

---

## 14. Differential transcript usage — `dtu_output/`

Produced when `--dtu` is set (Salmon/Kallisto). DEXSeq tests whether the
*proportions* of a gene's transcript isoforms shift between conditions — a gene
can be significant here even when its total expression (gene-level DE) is flat.

| File | Contents |
|---|---|
| `dtu_transcript_results.csv` | Per-transcript test: `groupID` (gene), `gene_name`, `gene_biotype`, `featureID` (transcript), `exonBaseMean`, `dispersion`, `stat`, `pvalue`, `padj` |
| `dtu_gene_qvalues.csv` | Per-gene q-value: `gene`, `gene_name`, `gene_biotype`, `gene_qvalue` (`perGeneQValue`, aggregating the gene's transcripts) |

Use `dtu_gene_qvalues.csv` (`gene_qvalue < 0.05`) to call genes with significant
usage changes, then `dtu_transcript_results.csv` to see which transcripts of
that gene drive the switch. Only multi-transcript, expressed genes are tested.

> Complementary to §13: IsoformSwitchAnalyzeR highlights *which* switch and its
> functional consequence; DEXSeq DTU is the formal statistical test of
> transcript-usage change.

---

## 15. Differential splicing — `diffsplice_output/`

Produced when `--diffsplice` is set. edgeR's `diffSpliceDGE` tests, for each
feature, whether its log-fold-change between conditions departs from the
gene's overall log-fold-change — i.e. differential *usage*. The feature is an
**exon** for STAR/HISAT2 and a **transcript** for Salmon/Kallisto.

| File | Contents |
|---|---|
| `diffsplice_exon_results.csv` *or* `diffsplice_transcript_results.csv` | Per-feature test: `GeneID`, `gene_name`, `gene_biotype`, `FeatureID`, `logFC`, an exon-level statistic, `P.Value`, `FDR`, `comparison` |
| `diffsplice_gene_results.csv` | Per-gene test (Simes' method across a gene's features): `GeneID`, `gene_name`, `gene_biotype`, `NExons`, `P.Value`, `FDR`, `comparison` |

Use `diffsplice_gene_results.csv` (`FDR < 0.05`) to find genes with a splicing
change, then the per-feature table to see which exon/transcript drives it. The
`comparison` column names the condition contrast. Only multi-feature genes are
tested. For STAR/HISAT2, the per-exon counts feeding this test are also kept in
`featurecounts_exon/<sample>.exon.featureCounts.txt` (`featureCounts -f`).

> Complementary to §11 (rMATS) and §14 (DEXSeq DTU): rMATS classifies splicing
> *events*, DEXSeq DTU tests transcript proportions, and edgeR diffSplice tests
> each feature's fold change against its gene — three independent views of
> alternative splicing.

---

## 16. Aggregated reports — `multiqc/`, `quarto_report/`

- **`multiqc/multiqc_report.html`** — single interactive page combining FastQC,
  fastp, alignment, RSeQC and featureCounts metrics for every sample. The best
  starting point for a run-wide quality overview, and for spotting outlier
  samples.
- **`quarto_report/analysis_report.html`** — an interactive analysis report.
  Alongside a MultiQC general-statistics summary, it pulls together the run's
  differential-expression and enrichment results: per-contrast significant-gene
  counts, interactive (plotly) volcano plots for DESeq2 and edgeR, the PCA /
  MDS / heatmap panels, a DESeq2-vs-edgeR agreement table, and searchable
  (DT) DESeq2 / edgeR / GSEA / gProfiler result tables. Sections for stages
  that did not run are simply omitted, so the report adapts to each run.

---

## 17. Execution metadata — `pipeline_info/`

Nextflow writes an execution `timeline`, `report` and `trace` here. Use them to
see which processes ran, how long they took, peak memory, and — if a run failed
— exactly which task and why. `software_versions.yml` (via MultiQC) records the
exact version of every tool used, for reproducibility and methods sections.

- **`run_summary.html`** — an end-of-run summary written automatically when the
  pipeline finishes (whether it succeeded or failed). It shows the run status,
  duration and command line; links the MultiQC report and every key result
  directory that was produced; and tabulates, for each process, the number of
  tasks, total job time, peak memory and mean CPU usage (aggregated from the
  execution trace). A concise version is also printed to the console at the end
  of the run. This is the quickest place to confirm a run finished cleanly and
  to jump to its outputs.

---

## Quick interpretation checklist

- [ ] MultiQC: no outlier samples, adapter content low after trimming
- [ ] Alignment: uniquely-mapped rate acceptable for your organism
- [ ] RSeQC: measured strandedness matches `--strandedness`
- [ ] PCA/MDS: replicates cluster by condition
- [ ] DE: genes significant in **both** DESeq2 and edgeR are the high-confidence set
- [ ] Enrichment: GSEA and gProfiler tell a consistent biological story
