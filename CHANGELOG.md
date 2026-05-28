# Changelog

All notable changes to **rnaseq-flow** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

No unreleased changes. Roadmap items are tracked in
[future_improvements.md](future_improvements.md); current candidates include
contamination / rRNA screening, a `--contrasts` parameter, an Arriba fusion
caller and a bundled CI test profile.

## [1.1.0] - 2026-05-25

### Added

- **Differential splicing — `--diffsplice`.** An edgeR `diffSpliceDGE` test of
  differential feature usage: exon-level on the STAR/HISAT2 route (via a new
  `featureCounts -f` per-exon count run) and transcript-level on the
  Salmon/Kallisto route. Adds the `DIFFSPLICE` and `FEATURECOUNTS_EXON`
  processes and a `diffsplice_output/` results directory.
- **Batch / covariate modelling.** The DESeq2/edgeR design is no longer fixed
  to `~ condition`: an optional `batch` samplesheet column is auto-detected and
  makes the model `~ batch + condition`, and a `--design` parameter accepts any
  formula (e.g. `~ sex + batch + condition`). Samplesheet validation checks that
  every design covariate is a real samplesheet column.
- **Gene-symbol annotation.** A new `GTF2GENEINFO` step parses the GTF into a
  gene_id/gene_name/gene_biotype table; the DESeq2, edgeR, DTU and diffSplice
  result tables now carry `gene_name` and `gene_biotype` columns, so they are
  readable without a separate ID-mapping step.
- **Per-contrast volcano plots.** DESeq2 and edgeR now write a
  `volcano_<A>_vs_<B>.png` for every contrast — log2 fold change vs
  −log10 p-value, with genes passing `padj`/`FDR < 0.05` and `|log2FC| > 1`
  coloured (up red, down blue) — alongside the existing MA / smear plots.
- **apeglm LFC shrinkage (DESeq2).** DESeq2 log2 fold changes are now shrunk
  with the apeglm estimator (`lfcShrink`), pulling low-count / high-variance
  estimates toward zero for better gene ranking and cleaner MA / volcano
  plots; the Wald `stat` / p-value / FDR are kept from the unshrunken fit so
  significance calling is unchanged. For non-reference contrasts the condition
  factor is releveled and the GLM refitted so the effect is a single
  coefficient apeglm can shrink. `DESEQ2` is now Conda/Wave-provisioned, as
  apeglm is not carried by the deseq2-only biocontainer.
- **Custom MultiQC report.** `assets/multiqc_config.yml` sets the report title,
  orders the modules in pipeline order and cleans sample names; a custom CSS
  theme and logo restyle the report in the rnaseq-flow identity.
- **Expanded Quarto analysis report.** The Quarto report
  (`quarto_report/analysis_report.html`, formerly QC-only `qc_report.html`)
  grew into a full analysis report: per-contrast significant-gene counts,
  interactive plotly volcano plots for DESeq2 and edgeR, the PCA / MDS /
  heatmap panels, a DESeq2-vs-edgeR agreement table, and searchable (DT)
  DESeq2 / edgeR / GSEA / gProfiler result tables — each section rendered only
  when its data is present. `QUARTO_REPORT` now takes the DE/enrichment result
  directories as inputs and is Conda/Wave-provisioned (it needs plotly + DT,
  which the previous `rocker/verse` image lacked).
- **Run-completion summary.** A `workflow.onComplete` handler now writes
  `pipeline_info/run_summary.html` at the end of every run (success or
  failure) — run status, duration and command line, links to the MultiQC
  report and every key result directory that was produced, and a per-process
  table of task count, total job time, peak memory and mean CPU usage
  aggregated from the execution trace (`trace.raw = true` keeps that trace
  machine-readable). A concise version is also printed to the console.
- **`CITATIONS.md`.** Every tool with its publication (verified against PubMed),
  grouped by pipeline stage, plus a ready-to-paste methods-paragraph template.
- **rnaseq-flow logo.** A vector logo (`assets/rnaseq-flow_logo.svg`, with a
  PNG companion) now appears on the MultiQC report, the `overview.html`
  workflow diagram and the Word user guide.

### Changed

- **edgeR contrasts now use the quasi-likelihood GLM throughout.** Every
  pairwise contrast is tested with `glmQLFTest`; the previous `exactTest` path
  for non-reference pairs has been removed. `exactTest` ignores the design
  matrix, so the GLM path is what makes the new covariate modelling correct and
  keeps all contrasts mutually consistent.
- **Kallisto now honours `--strandedness`.** `forward` / `reverse` are mapped
  to kallisto's `--fr-stranded` / `--rf-stranded`. Kallisto produces no BAM for
  RSeQC, so strandedness cannot be auto-inferred — `auto` runs
  library-type-agnostic and logs a warning.
- **`overview.html`** gained the diffSplice node and metro-map stations, and
  the metro map was reworked with compact, evenly-paired interchanges.

## [1.0.0] - 2026-05-23

First release of the modular bulk RNA-seq pipeline.

### Added

#### Core workflow

- Three run modes from a single entry point: `--download_refs` (fetch genome
  FASTA + GTF from Ensembl), `--build_indices` (STAR / HISAT2 / Salmon /
  Kallisto), and the full `RNASEQ` analysis workflow.
- Read QC and trimming with FastQC and fastp.
- Four interchangeable aligners selected with `--aligner`: STAR and HISAT2
  (genome alignment) and Salmon and Kallisto (pseudo-alignment).
- Post-alignment QC with RSeQC, CPM-normalised BigWig coverage tracks
  (deepTools) and gene-level quantification with featureCounts.
- Differential expression with DESeq2 and edgeR, run in parallel as two
  independent callers on the same gene counts.
- Alternative splicing (rMATS), gene-fusion detection (STAR-Fusion) and
  functional enrichment (fgsea GSEA and gprofiler2).
- Aggregated reporting with MultiQC and a Quarto QC report.
- `--stop_at` staging (`preQC`, `postQC`, `DE`, `GSEA`).

#### Transcript-level analysis

- `tximport` step that summarises Salmon/Kallisto transcript quantification to
  gene level, so the pseudo-aligners feed DESeq2/edgeR with length-aware
  normalisation.
- IsoformSwitchAnalyzeR for transcript isoform-switch detection (Salmon).
- DEXSeq differential transcript usage (DTU) test, opt-in via `--dtu`
  (Salmon/Kallisto), flagging genes whose isoform proportions shift between
  conditions.

#### Reproducibility & robustness

- Per-process resource configuration (`conf/base.config`) with `check_max`
  ceilings and a one-retry-with-doubled-resources policy.
- `conf/modules.config` `publishDir` rules and a `stub` block in every process
  for fast `-stub-run` dry-runs.
- Ensembl release pinning (`--download_release`), writing references into a
  versioned `v<release>` subfolder.
- Per-sample strandedness auto-detection (`--strandedness auto`): RSeQC
  `infer_experiment` results are fed straight into featureCounts.
- Fail-fast samplesheet validation (required columns, unique sample ids, FASTQ
  existence, >= 2 conditions and >= 2 replicates per condition).
- `nextflow_schema.json` enabling `--help` and parameter typo-detection.
- Wave on-the-fly container provisioning for Conda-declared processes.

#### Documentation

- `README.md`, `USAGE.md` and `OUTPUTS.md`.
- `overview.html` — an interactive workflow diagram (DAG and metro-map views).
- `rnaseq-flow_User_Guide.docx` — a formatted user guide.
