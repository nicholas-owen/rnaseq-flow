# rnaseq-flow Output Guide

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
├── deseq2_output/      Differential expression (DESeq2)
│   └── reproduce/      Scripts, data and objects that redraw every figure (§18)
├── edger_output/       Differential expression (edgeR)
│   └── reproduce/
├── gsea_output/        Gene-set enrichment (fgsea)
│   └── reproduce/
├── gprofiler_output/   GO / pathway over-representation (gprofiler2)
│   └── reproduce/
├── rmats_output/       Alternative splicing events, one subdirectory per contrast
├── star_fusion/        Gene-fusion predictions
├── isoform_switch/     Isoform-switching analysis
├── dtu_output/         Differential transcript usage (DEXSeq)
├── diffsplice_output/  Differential splicing (edgeR diffSpliceDGE)
├── multiqc/            Aggregated MultiQC report (rnaseq-flow_multiqc_report.html)
├── quarto_report/      Interactive Quarto analysis report
└── pipeline_info/      Run manifest, software versions, run summary, Nextflow trace
```

Which directories appear depends on the aligner and options you chose.

> **Count-matrix runs (`--counts`).** Starting from a pre-computed count matrix
> skips read QC, alignment and quantification, so those directories
> (`fastqc/`, `fastp/`, `star/`|`hisat2/`|`salmon/`|`kallisto/`, `rseqc/`,
> `bam_to_bigwig/`, `featurecounts/`, `tximport/`) and the MultiQC report are
> **not** produced. You get the differential-expression outputs onward
> (`deseq2_output/`, `edger_output/`, and any enrichment), plus a Quarto report
> whose QC section is omitted. See [USAGE.md §4.7](USAGE.md#47-starting-from-a-count-matrix---counts).

---

## Recommended reading order

1. **`quarto_report/analysis_report.html`, Run overview**: what was run, on
   what, against which reference and baseline. Read this first if the results
   directory is not your own.
2. **`multiqc/rnaseq-flow_multiqc_report.html`**: one page for whole-run QC.
3. **Alignment logs**: confirm mapping rates are acceptable.
4. **`rseqc/`**: confirm strandedness and library quality.
5. **`deseq2_output/pca_plot.svg`**: do samples group by condition?
6. **DE tables**: the differentially expressed genes.
7. **Enrichment / splicing / fusions**: biological interpretation.

---

## 1. Read QC: `fastqc/`, `fastp/`

**`fastqc/`**: `*_fastqc.html` per sample. Check per-base quality (should stay
in the green), adapter content, and duplication. Raw RNA-seq normally shows
some duplication and a biased first ~12 bp. That is expected.

**`fastp/`**: trimmed FASTQ files plus `*.fastp.html` / `*.fastp.json`. The
report shows reads before/after filtering and adapter removal. A large drop in
read count means aggressive filtering: inspect the input quality.

> All QC metrics are also aggregated in `multiqc/`, which is usually easier to compare
> across samples there.

---

## 2. Alignment: `star/` or `hisat2/`

| File | Description |
|---|---|
| `*.Aligned.sortedByCoord.out.bam` / `*.bam` | Coordinate-sorted alignments |
| `*.bai` | BAM index (for IGV and downstream tools) |
| `*.Log.final.out` (STAR) | Alignment summary statistics |
| `*.hisat2.summary.log` (HISAT2) | Alignment summary statistics |
| `*.ReadsPerGene.out.tab` (STAR) | STAR's own gene counts (not used downstream) |

**Reading the STAR `Log.final.out`**: the key line is *Uniquely mapped reads
%*:

| Uniquely mapped % | Interpretation |
|---|---|
| > 85% | Excellent |
| 70–85% | Acceptable |
| < 70% | Investigate: contamination, wrong genome, or poor quality |

Also watch *% of reads mapped to multiple loci* (high values suggest rRNA or
repetitive contamination) and *% of reads unmapped: too short* (often adapter
or quality problems).

---

## 3. Transcript quantification: `salmon/`, `kallisto/`

One sub-directory per sample.

**Salmon**: `<sample>/quant.sf`, columns:

| Column | Meaning |
|---|---|
| `Name` | Transcript ID |
| `Length` / `EffectiveLength` | Transcript length / length corrected for fragment bias |
| `TPM` | Transcripts Per Million: normalised abundance, comparable across samples |
| `NumReads` | Estimated reads assigned to the transcript |

**Kallisto**: `<sample>/abundance.tsv` with `target_id`, `length`,
`eff_length`, `est_counts`, `tpm`.

Use `TPM` to compare expression directly. The estimated counts
(`NumReads`/`est_counts`) are summarised to gene level by the pipeline's
tximport step (see §6) and fed straight into DESeq2/edgeR.

---

## 4. Alignment QC: `rseqc/`

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
inference is not just reported: `*.strandedness.txt` records the per-sample
verdict and it is fed directly into that sample's featureCounts run, so the
correct `-s` flag is used without you having to know the library type. Setting
`--strandedness` explicitly overrides the inference; if you do, check it
against `*.strandedness.txt`, since a wrong strand setting roughly halves your
gene counts.

**`read_distribution.txt`**: most tags should fall in exonic regions (CDS +
UTRs). A high intronic/intergenic fraction suggests DNA contamination or
incomplete annotation.

**Gene-body coverage**: should be a flat plateau. A strong 3' skew indicates
RNA degradation; a 5' skew can indicate library-prep bias.

---

## 5. Coverage tracks: `bam_to_bigwig/`

`<sample>.bw`: CPM-normalised coverage in BigWig format. Load into IGV or the
UCSC Genome Browser to inspect coverage at specific loci. CPM normalisation
makes tracks comparable across samples of different depth.

---

## 6. Gene counts: `featurecounts/` and `tximport/`

Gene-level counts are the input to DESeq2 and edgeR. How they are produced
depends on the aligner.

**`featurecounts/` (STAR / HISAT2).** `<sample>.featureCounts.txt` holds
gene-level counts: columns `Geneid`, `Chr`, `Start`, `End`, `Strand`,
`Length`, and a final column of raw read counts.
`<sample>.featureCounts.txt.summary` breaks down assigned vs unassigned reads
(and *why* reads were unassigned).

**`tximport/` (Salmon / Kallisto).** `tximport_gene_counts.csv` is the
gene-by-sample count matrix obtained by summarising transcript-level
quantification to gene level (genes × samples). `txi.rds` is the full tximport
object (counts, abundances and a transcript-length matrix), which DESeq2
imports with `DESeqDataSetFromTximport` and edgeR imports as a length offset,
giving more accurate normalisation than counts alone.

Either way these counts are **not** normalised for library size: do not
compare them directly between samples; that is what DESeq2/edgeR do.

---

## 7. Differential expression: `deseq2_output/`

DESeq2 and edgeR (§8) both run, **in parallel**, on the *same* gene counts
(featureCounts for STAR/HISAT2, tximport for Salmon/Kallisto) as two independent
callers: neither is downstream of the other. This section covers DESeq2.

One results CSV per pairwise contrast, named
`deseq2_results_<A>_vs_<B>.csv`, where B is the baseline
([USAGE.md §4.5](USAGE.md#45-which-condition-results-are-measured-against---reference_level))
and the fold change is A relative to B.

| Column | Meaning |
|---|---|
| `gene_id` | The stable gene identifier from the GTF. Symbols are neither unique nor stable across annotation releases, so join on this |
| `gene_name` | Gene symbol, parsed from the GTF |
| `gene_biotype` | Gene biotype (e.g. protein_coding, lncRNA) |
| `baseMean` | Mean normalised count across all samples |
| `log2FoldChange` | Effect size, log2(A / B), **apeglm-shrunken** (see note) |
| `log2FoldChange_MLE` | The unshrunken maximum-likelihood fold change, kept beside the shrunken one so the size of the shrinkage can be seen rather than inferred |
| `lfcSE` | Posterior SD of the shrunken log2 fold change (apeglm) |
| `stat` | Wald test statistic (from the unshrunken fit) |
| `pvalue` | Raw p-value |
| `padj` | **Benjamini–Hochberg adjusted p-value (FDR)**: use this |

**LFC shrinkage.** The `log2FoldChange` is shrunk with the apeglm estimator
(DESeq2 `lfcShrink`): low-count, high-variance genes are pulled toward zero,
giving more reliable effect sizes and cleaner MA / volcano plots. Shrinkage
changes only the effect-size estimate: the `stat`, `pvalue` and `padj` come
from the unshrunken Wald test, so significance calling is unaffected.

**Calling a gene significant**: a common threshold is `padj < 0.05` and
`|log2FoldChange| > 1` (a 2-fold change). `padj` may be `NA` when DESeq2
filters a gene for low counts or flags it as an outlier. That is normal.

Plots. Every figure is written twice, as `.svg` and as a 300 dpi `.png`: the
vector copy for publication and figure assembly, the raster for slides and
email. Each can be redrawn and adapted from `reproduce/` (§18).

| File | Use |
|---|---|
| `pca_plot` | Sample clustering: replicates of a condition should group together; if not, suspect batch effects or mislabelling |
| `maplot_<A>_vs_<B>` | Fold change vs mean expression; significant genes highlighted |
| `volcano_<A>_vs_<B>` | Volcano plot: log2 fold change vs −log10 p-value; genes with `padj < 0.05` and `|log2FC| > 1` are coloured (up red, down blue) |
| `shrinkage_ma_<A>_vs_<B>` | The same contrast before and after apeglm, so you can see how far shrinkage moved the fold changes rather than take it on faith |
| `heatmap_top_var` | The 20 most variable genes across samples, variance-stabilised |

---

## 8. Differential expression: `edger_output/`

edgeR runs in parallel with DESeq2 on the same gene-count input: neither is
downstream of the other. It is an independent cross-check; its results are
**not** consumed by any later step (GSEA and gProfiler use the DESeq2 results
only). Results CSVs are `edger_results_<A>_vs_<B>.csv`.

| Column | Meaning |
|---|---|
| `gene_id` | The stable gene identifier from the GTF |
| `gene_name` | Gene symbol, parsed from the GTF |
| `gene_biotype` | Gene biotype (e.g. protein_coding, lncRNA) |
| `logFC` | log2 fold change |
| `logCPM` | Average log2 counts per million |
| `F` | Quasi-likelihood F-test statistic |
| `PValue` | Raw p-value |
| `FDR` | **Adjusted p-value**: use this |

`mds_plot` is edgeR's sample-similarity plot (analogous to the PCA);
`smear_<A>_vs_<B>` plots logFC vs logCPM with significant genes marked;
`volcano_<A>_vs_<B>` plots logFC vs −log10 p-value, colouring genes with
`FDR < 0.05` and `|logFC| > 1` (up red, down blue). Each is written as `.svg`
and `.png`, and each can be redrawn from `reproduce/` (§18).

> Genes called significant by **both** DESeq2 and edgeR are the most robust
> hits. Modest disagreement near the significance threshold is expected.

---

## 9. Gene-set enrichment: `gsea_output/`

fgsea is run on the **ranked** DESeq2 gene list for each contrast (genes ranked
by the `stat` column), so it captures coordinated, subtle shifts a hard cutoff
would miss.

`gsea_stats_<contrast>.csv` columns:

| Column | Meaning |
|---|---|
| `pathway` | Gene-set name |
| `pval` / `padj` | Enrichment p-value / FDR |
| `ES` | Enrichment score |
| `NES` | **Normalised enrichment score**: sign gives direction |
| `size` | Number of genes from the set found in the data |

Interpretation: `padj < 0.05` is significant; positive `NES` = the pathway is
up in condition A, negative `NES` = up in B.

| File | Contents |
|---|---|
| `gsea_stats_<contrast>.csv` | The table above |
| `gsea_ranks_<contrast>.csv` | The ranked gene list the test was run on |
| `gsea_plot_<contrast>` | Enrichment-table plot of the top up and down pathways (`.svg` and `.png`) |
| `gsea_dotplot_<contrast>` | The most significant gene sets by NES and padj (`.svg` and `.png`) |
| `reproduce/` | Redraws both figures (§18). GSEA is the case that most needs it: `plotGseaTable` requires gene-set membership from the GMT, a pipeline input that is never published with the results, so the plot object is saved here |

---

## 10. GO / pathway over-representation: `gprofiler_output/`

gprofiler2 tests the **significant gene lists from the DESeq2 results** for
enriched GO terms and pathways, separately for up- and down-regulated genes.

| File | Contents |
|---|---|
| `gprofiler_UP_<contrast>.csv` | Terms enriched among up-regulated genes |
| `gprofiler_DOWN_<contrast>.csv` | Terms enriched among down-regulated genes |
| `gostplot_UP/DOWN_<contrast>` | Manhattan-style enrichment plot (`.svg` and `.png`) |
| `reproduce/` | Redraws the plots (§18). `gostplot` needs the full gost result object, not just the table, so that is saved as `.rds` beside the gzipped CSVs |

Key columns: `term_name`, `source` (GO:BP/MF/CC, KEGG, REACTOME…),
`p_value` (already multiple-testing corrected by g:SCS), `intersection_size`
(your genes in the term).

> GSEA (§9) uses the *whole ranked list*; gProfiler uses a *thresholded list*.
> They answer slightly different questions. Agreement between them strengthens
> a conclusion.

---

## 11. Alternative splicing: `rmats_output/`

One sub-directory per condition pair, named `<condition>_vs_REF` to match the
differential-expression tables, with `REF` as the denominator.

> **Changed after v1.5.1.** Directories were previously named `REF_vs_<condition>`
> and `IncLevelDifference` carried the opposite sign, so a positive value meant
> higher inclusion in `REF`, the reverse of what a positive `log2FoldChange`
> means in the DE tables beside it. **Splicing results produced before this change
> point the other way**, and the files give no indication of which convention they
> follow beyond the directory name. Do not compare results across that boundary
> without checking.

rMATS reports five event types:

| Code | Event |
|---|---|
| `SE` | Skipped exon |
| `MXE` | Mutually exclusive exons |
| `A3SS` | Alternative 3' splice site |
| `A5SS` | Alternative 5' splice site |
| `RI` | Retained intron |

For each type there are two files:

- `*.MATS.JC.txt`: junction counts only (reads spanning the splice junction).
- `*.MATS.JCEC.txt`: junction counts **plus** reads on the exon body.

Key columns:

| Column | Meaning |
|---|---|
| `GeneID`, `geneSymbol` | Gene |
| `IncLevel1`, `IncLevel2` | Inclusion levels (PSI) in condition 1 and 2 |
| `IncLevelDifference` | PSI difference (1 − 2), the effect size |
| `PValue`, `FDR` | Significance of the difference |

**Reading the sign.** Condition 1 is the first name in the directory, condition 2
the second. In `NaCl_vs_REF`, `IncLevelDifference` is PSI(NaCl) − PSI(REF), so a
**positive value means higher inclusion in NaCl**: the same direction a positive
`log2FoldChange` means in `deseq2_output/deseq2_results_NaCl_vs_REF.csv`.

PSI is a proportion, so the difference is bounded at ±1. A value of 0.2 means
twenty percentage points of the transcript pool shifted. It is not a fold change
and should not be compared to one.

A common cutoff: `FDR < 0.05` and `|IncLevelDifference| > 0.1`.

> **`summary.txt` uses different criteria.** rMATS builds it with an
> `|IncLevelDifference|` cutoff of **0**, so its `SignificantEvents*` columns
> count every FDR-significant event regardless of effect size. Those numbers will
> be larger, often much larger, than the cutoff above produces. Recompute from
> the `*.MATS.JC.txt` files rather than quoting `summary.txt`.

---

## 12. Gene fusions: `star_fusion/`

`<sample>.star-fusion.fusion_predictions.tsv` (full) and
`.abridged.tsv` (summary). Key columns:

| Column | Meaning |
|---|---|
| `#FusionName` | The two partner genes, e.g. `GENE1--GENE2` |
| `JunctionReadCount` | Reads spanning the fusion breakpoint |
| `SpanningFragCount` | Read pairs flanking the breakpoint |
| `LeftBreakpoint` / `RightBreakpoint` | Genomic coordinates |
| `LargeAnchorSupport` | Whether reads have long anchors (more reliable) |
| `FFPM` | Fusion fragments per million, normalised support |

Higher junction + spanning support and `YES` large-anchor support indicate
more confident calls. Always validate fusion candidates against known biology
and, ideally, orthogonal evidence.

---

## 13. Isoform switching: `isoform_switch/`

| File | Contents |
|---|---|
| `isoform_switches.csv` | Top genes with significant isoform switches |
| `switch_plot_<gene>.png` | Per-gene isoform usage across conditions |
| `switchList_analyzed.rds` | The full R object for custom downstream analysis |

A switch is a gene where the *dominant* transcript isoform changes between
conditions even if total gene expression does not: biologically important and
invisible to gene-level DE. Sort `isoform_switches.csv` by q-value.

---

## 14. Differential transcript usage: `dtu_output/`

Produced when `--dtu` is set (Salmon/Kallisto). DEXSeq tests whether the
*proportions* of a gene's transcript isoforms shift between conditions. A gene
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

## 15. Differential splicing: `diffsplice_output/`

Produced when `--diffsplice` is set. edgeR's `diffSpliceDGE` tests, for each
feature, whether its log-fold-change between conditions departs from the
gene's overall log-fold-change, i.e. differential *usage*. The feature is an
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
> each feature's fold change against its gene: three independent views of
> alternative splicing.

---

## 16. Aggregated reports: `multiqc/`, `quarto_report/`

- **`multiqc/rnaseq-flow_multiqc_report.html`**: single interactive page
  combining FastQC, fastp, alignment, RSeQC and featureCounts metrics for every
  sample. The best starting point for a run-wide quality overview, and for
  spotting outlier samples. Samples are shown by samplesheet id, with a Source
  FASTQ column mapping each back to the file it came from, and a Software
  Versions section lists every tool that ran. The filename carries the
  `rnaseq-flow_` prefix because MultiQC prefixes its output with the `title:`
  from `assets/multiqc_config.yml`.
- **`quarto_report/analysis_report.html`**: an interactive analysis report
  that brings the run together in one document. It opens with a **Run
  overview** (what ran it, what went in, the library layout, the aligner, the
  baseline, and the reference set with its source and release, read from
  `pipeline_info/run_manifest.json`), a **Sample design** table, and a
  **Requested but not performed** notice that renders only when an analysis was
  asked for and skipped. Then the results: a MultiQC general-statistics summary,
  per-contrast significant-gene counts, a shrinkage before/after view,
  interactive (plotly) volcano and MA plots for DESeq2 and edgeR, PCA, MDS and
  the top-variable-gene heatmap with dendrograms, a DESeq2-vs-edgeR agreement
  table, searchable (DT) result tables, a GSEA dot plot and a gProfiler
  Manhattan plot. It closes with **Files and locations**, an index of where the
  run wrote everything. Sections for stages that did not run are omitted.

  Its **Alternative splicing** section reads `rmats_output/` directly and covers
  the ground described in §11: an events-per-type count recomputed from the
  `MATS.JC.txt` files rather than taken from `summary.txt`, a JC-against-JCEC
  agreement table, a PSI-difference volcano coloured by event type, a stacked
  bar of significant events, per-replicate PSI for the top 16 events, a
  JC-against-JCEC concordance scatter, and a searchable results table. It ends
  by counting the genes whose splicing changed while their total expression did
  not, which the DE tables miss by construction. Two things the result files do
  not show are stated above the results: rMATS runs without `--libType`, so a
  stranded run is flagged, and a design with fewer than three replicates per
  group is flagged as well.

---

## 17. Execution metadata: `pipeline_info/`

| File | Contents |
|---|---|
| `run_manifest.json` | How the run was configured, machine-readable: pipeline version and commit, the command line, the aligner, the baseline (`--reference_level`), the design, the samples per condition, the notices, and the reference set with its release where known. Written before any process runs, so it records what was asked for; the report derives what was produced from the result directories that exist |
| `software_versions_mqc.yml` | Every tool and version collected at run time, the same table MultiQC shows. Use it for methods sections |
| `run_summary.html` | See below |
| `execution_report_*.html`, `execution_timeline_*.html`, `execution_trace_*.txt` | Nextflow's own per-run resource report, timeline and trace, timestamped per launch. Use them to see which processes ran, how long they took, peak memory, and if a run failed, exactly which task and why |

The reference set's own provenance is not here: `--download_refs` writes
`reference_metadata.json` beside the reference files
([USAGE.md §2](USAGE.md#2-mode-1-download-references)), and the manifest
reads it at launch to fill in the annotation source and release. Keep it with
the references.

- **`run_summary.html`**: an end-of-run summary written automatically when the
  pipeline finishes (whether it succeeded or failed). It shows the run status,
  duration and command line; links the MultiQC report and every key result
  directory that was produced; and tabulates, for each process, the number of
  tasks, total job time, peak memory and mean CPU usage (aggregated from the
  execution trace). A concise version is also printed to the console at the end
  of the run. This is the quickest place to confirm a run finished cleanly and
  to jump to its outputs.

---

## 18. Redrawing and adapting figures: `reproduce/`

Every figure the pipeline publishes can be redrawn, and adapted, from the
results directory alone: no rerunning, no access to the original compute, no
digging through `work/`. Each of `deseq2_output/`, `edger_output/`,
`gsea_output/` and `gprofiler_output/` carries a `reproduce/` subfolder holding
a standalone R script per figure type, gzipped copies of the tables behind each
figure, and the R objects a plot genuinely needs, always the smallest that will
do.

The split is that the parent folder is for reading and `reproduce/` is for
running. Parent CSVs stay uncompressed so they open by double-click; the
`reproduce/` copies are gzipped, which R reads transparently.

```
deseq2_output/reproduce/
├── pca_plot.R                              redraws pca_plot from pca_data.csv.gz
├── volcano.R                               redraws volcano_<A>_vs_<B> from the results CSV
├── maplot.R                                redraws maplot_<A>_vs_<B> from the .rds
├── shrinkage_ma.R                          redraws shrinkage_ma_<A>_vs_<B> from the results CSV, which carries both fold changes
├── heatmap_top_var.R                       redraws heatmap_top_var from the two heatmap CSVs
├── pca_data.csv.gz                         PCA coordinates and variance explained
├── heatmap_top_var.csv.gz                  variance-stabilised matrix, top 20 genes
├── heatmap_top_var_annotation.csv.gz       sample-to-condition map for the colour bar
├── deseq2_results_<A>_vs_<B>.csv.gz        the results table, gzipped
└── deseq2_results_<A>_vs_<B>.rds           the DESeq2 results object the MA plot needs, with $res_mle, the unshrunken fit

edger_output/reproduce/
├── mds_plot.R                              redraws mds_plot from mds_plot.rds
├── smear.R                                 redraws smear_<A>_vs_<B> from the .rds and the results CSV
├── volcano.R                               redraws volcano_<A>_vs_<B> from the results CSV
├── mds_data.csv.gz                         MDS coordinates
├── mds_plot.rds                            the small MDS object, not the whole DGEList
├── edger_results_<A>_vs_<B>.csv.gz         the results table, gzipped
└── edger_results_<A>_vs_<B>.rds            the test result object the smear plot needs

gsea_output/reproduce/
├── gsea_plot.R                             redraws gsea_plot_<A>_vs_<B> from the .rds
├── gsea_dotplot.R                          redraws gsea_dotplot_<A>_vs_<B> from the stats CSV
├── gsea_plot_<A>_vs_<B>.rds                the plotGseaTable object; gene-set membership lives only here
├── gsea_stats_<A>_vs_<B>.csv.gz            the results table, gzipped
└── gsea_ranks_<A>_vs_<B>.csv.gz            the ranked gene list the test was run on

gprofiler_output/reproduce/
├── gostplot.R                              redraws gostplot_UP and gostplot_DOWN from the .rds files
├── gprofiler_UP_<A>_vs_<B>.csv.gz          the results tables, gzipped
├── gprofiler_DOWN_<A>_vs_<B>.csv.gz
├── gprofiler_UP_<A>_vs_<B>.rds             the full gost result objects, which gostplot needs
└── gprofiler_DOWN_<A>_vs_<B>.rds
```

One script per figure type, one gzipped CSV per published table, and an
`.rds` only where a plot cannot be rebuilt from a CSV: the MDS object rather
than the whole `DGEList`, which is 28 times larger even at 2,000 genes; the
fgsea plot object, because `plotGseaTable` needs gene-set membership from the
GMT, a pipeline input never published with the results; and the gost result,
because `gostplot` takes the object rather than the table.

The scripts are written to be changed, not merely run. Contrast, cutoffs,
colours and figure dimensions are named constants in the first lines, so
adjusting a legend does not mean reading a `ggplot` chain:

```r
CONTRAST      <- "NaCl_vs_REF"
PADJ_CUT      <- 0.05
LFC_CUT       <- 1
WIDTH         <- 8
HEIGHT        <- 6.5
FIG_RES       <- 300
OUTPUT_PREFIX <- "repro_"
```

Run one from a terminal:

```bash
cd results/deseq2_output/reproduce
Rscript volcano.R
```

or in RStudio: open the file, then Session > Set Working Directory > To Source
File Location, and Source. Regenerated figures are prefixed `repro_` and get a
matching `.info.txt` recording the settings used, the environment that produced
them and where they came from, so a regenerated figure is never mistaken for
the pipeline's own once it leaves the folder.

> `isoform_switch/` is outside these conventions. IsoformSwitchAnalyzeR writes
> its own figures, so they are PNG-only and have no `reproduce/` folder.

---

## Quick interpretation checklist

- [ ] MultiQC: no outlier samples, adapter content low after trimming
- [ ] Alignment: uniquely-mapped rate acceptable for your organism
- [ ] RSeQC: measured strandedness matches `--strandedness`
- [ ] PCA/MDS: replicates cluster by condition
- [ ] DE: genes significant in **both** DESeq2 and edgeR are the high-confidence set
- [ ] Enrichment: GSEA and gProfiler tell a consistent biological story
