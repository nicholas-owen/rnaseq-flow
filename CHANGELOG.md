# Changelog

All notable changes to **rnaseq-flow** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]


### Fixed

- **GSEA now ranks genes by the identifier the gene sets actually use.** `gsea.R`
  ranked genes by their Ensembl gene IDs, but MSigDB `--gmt` files (and the
  pipeline's own downloaded gene sets) are keyed by gene *symbol*, so the ranked
  list and the pathways overlapped at essentially zero genes — fgsea returned an
  empty table and the report's GSEA section silently showed "no data". Ranking
  now uses whichever identifier overlaps the gene sets (the annotated gene
  symbol, falling back to the Ensembl ID for an Ensembl-keyed GMT), collapses
  duplicate symbols to their most-extreme statistic, and warns (rather than
  failing silently) when the overlap is too low to trust.


Roadmap items are tracked in
[future_improvements.md](future_improvements.md); current candidates include
contamination / rRNA screening, a `--contrasts` parameter, an Arriba fusion
caller and a bundled CI test profile.

## [1.4.0] - 2026-07-17

### Added

- **Count-matrix entry point (`--counts`).** The pipeline can now start from a
  pre-computed gene count matrix, skipping QC and alignment and entering directly
  at differential expression (DESeq2, edgeR, GSEA, gProfiler, Quarto report). This
  supports re-analysing published or collaborator count tables (e.g. GEO
  supplementary files). Two layouts are auto-detected: featureCounts wide output
  and a plain gene-id + samples matrix (tab or comma delimited), including
  gzipped and header-less matrices. Requires `--input` (for the design) and
  `--gtf` (for gene annotation); the `R1`/`R2` samplesheet columns become
  optional. Non-integer input (normalised / TPM / FPKM or Salmon/Kallisto
  estimates) and Entrez/RefSeq-keyed matrices are rejected with clear guidance.
  See [USAGE.md §4.6](USAGE.md#46-starting-from-a-count-matrix---counts).

## [1.3.0] - 2026-07-08

Nextflow 26 support. Nextflow's strict language parser — the default from
Nextflow 25.10 and used by `nextflow lint` — rejected several constructs the
pipeline relied on, so it would not run on current Nextflow. The config and
scripts are migrated to the strict language and now pass `nextflow lint` with
zero errors and zero warnings; the full workflow was verified end-to-end on
Nextflow 26.04 (STAR, Salmon and reference-download paths).

**Minimum Nextflow is now 25.10** (`nextflowVersion = '!>=25.10.0'`, was
`22.10.1`). The strict-language constructs used here — `process.resourceLimits`,
the entry-workflow `onComplete:` section, and no top-level script statements —
require the new parser, which is the default from 25.10.

### Changed

- **Resource ceilings now use `process.resourceLimits`** (native, Nextflow
  24.04+) instead of the `check_max()` helper. The strict config parser forbids
  function definitions in a config file, so `check_max()` (and the
  `hisat2_build_mem()` helper added in 1.2.0) could not be defined there.
  `--max_cpus` / `--max_memory` / `--max_time` still cap every request exactly as
  before; the HISAT2 index-build memory estimate is now an inline `ext`
  closure in `conf/base.config`. No change to the resources any process requests.

- **Entry-workflow `onComplete:` section.** The run-completion summary
  (`pipeline_info/run_summary.html`) was written by a top-level
  `workflow.onComplete { }` handler, which the strict parser does not allow. It
  is now the entry workflow's `onComplete:` section; behaviour is unchanged
  (verified: the summary is written on both success and failure).

- **Trace-file timestamps are inlined.** The shared `def trace_timestamp`
  variable is gone (the strict parser forbids config variable declarations); each
  of the timeline / report / trace filenames computes the timestamp inline.

### Fixed

- **Top-level script statements moved into the workflow.** The schema-driven
  `--help`, parameter typo-detection and required-input checks ran as top-level
  statements in `main.nf`, which the strict parser rejects. They are now a
  `checkParameters()` function called from the entry workflow's `main:` section.
  `--help` and typo-detection behave exactly as before.

- **Strict-syntax script cleanups.** Removed a `while` loop (no longer supported
  — replaced with a range iterator), C-style multi-variable declarations, and
  `;`-joined statements in the run-summary code; converted the run-summary
  formatting closures to top-level functions (the strict parser does not resolve
  a closure variable called from inside another closure).

- **Deprecation warnings cleared.** `Channel.*` factory access → `channel.*`;
  implicit closure parameters (`it`) → explicit parameters; unused closure
  parameters prefixed with `_`; single workflow emit given as an unnamed
  expression. `nextflow lint .` now reports no warnings.

## [1.2.0] - 2026-07-08

A reference-download and index-building release. `--download_refs` selected the
wrong GTF and its "current release" path had stopped working against Ensembl;
both are fixed, the workflow now fetches the transcriptome it was always missing,
and the HISAT2 index build no longer runs out of memory on mammalian genomes.

**Re-download your reference set.** Any run whose references came from
`--download_refs` before this version used a patch/haplotype annotation
(`*.chr_patch_hapl_scaff.gtf.gz`) against a primary-assembly genome, and should
be regarded as provisional. See *Fixed* below for what that did and did not
affect.

### Fixed

- **`--download_refs` downloaded the wrong GTF.** `assets/download_refs.py`
  chose the annotation by *excluding* known variants. Ensembl lists them
  alphabetically — `abinitio`, `chr`, `chr_patch_hapl_scaff`, then the canonical
  file — so the first surviving candidate was always
  `*.chr_patch_hapl_scaff.gtf.gz`, a patch/haplotype annotation that does not
  describe the `dna.primary_assembly` genome downloaded alongside it. On Ensembl
  release-116 human that GTF spans 528 contigs against the FASTA's 70 (458 with
  no sequence at all), carries 86,411 gene records instead of 78,941, and raises
  duplicated gene symbols from 484 to 3,384 (`HLA-A` appears 8× instead of
  once). Gene `gene_id`s stay unique and the surplus genes receive no reads, so
  gene-level counts were not themselves corrupted — but STAR's index disagreed
  with its annotation, every annotation-derived table carried phantom genes, and
  symbol-keyed steps saw colliding names.

  The GTF and the transcriptome are now selected by **the assembly named in the
  genome FASTA**, and anything that does not match exactly one file is a hard
  error rather than a guess. Two simpler rules were tried and rejected: the
  trailing number in a GTF filename is not always the Ensembl release
  (release-116 ships `Saccharomyces_cerevisiae.R64-1-1.63.gtf.gz`), and a single
  release directory can hold GTFs for more than one assembly (release-110 has
  both `Drosophila_melanogaster.BDGP6.32.110.gtf.gz` and `...BDGP6.46.110.gtf.gz`
  while the genome is `BDGP6.46`).

- **`--download_refs` no longer works against `pub/current_gtf/`.** Ensembl has
  removed that alias — it now returns 404 over both HTTP and HTTPS, although
  `pub/current_fasta/` still resolves. A `current` download therefore failed to
  find any annotation at all. `current` is now resolved to a concrete release
  number up front (Ensembl REST `/info/data`, falling back to scraping `pub/`),
  and the genome, GTF and transcriptome are all fetched from that same
  `release-<N>/` directory.

- **Partial reference sets were reported as success.** A failed or missing
  download left `references/` holding whatever had arrived, and the process
  exited 0. `download_refs.py` now exits non-zero unless all three files are
  found, downloaded and verified, and `DOWNLOAD_REFS` declares `fasta` and `gtf`
  as required rather than `optional` outputs.

- **`--download_source ncbi` was a silent no-op**, producing a directory that
  contained only a log file. It is now a hard error naming the supported
  alternatives.

- **`HISAT2_BUILD` ran out of memory on mammalian genomes.** A splice-aware
  HISAT2 index (`hisat2-build --ss --exon`) needs on the order of 200 GB for
  human, but the process only requested the 48 GB of the `process_high` label,
  so the build was OOM-killed, retried once, and failed. It now requests memory
  estimated from the GTF size — roughly `8.GB + 45.GB` per GB of uncompressed
  annotation (about 204 GB for human, 158 GB for mouse, 9 GB for yeast), still
  capped by `--max_memory`. When the granted memory is below what a splice-aware
  build needs, `HISAT2_BUILD` now drops `--ss`/`--exon` and builds a
  non-splice-aware index (with a clear warning) instead of being killed. A new
  `--hisat2_build_memory` parameter overrides the estimate.

### Added

- **`--hisat2_build_memory` parameter.** Sets the memory a splice-aware HISAT2
  index build is assumed to need (e.g. `200.GB`); default is estimated from the
  GTF size. Doubles as the threshold below which the build degrades to a
  non-splice-aware index.

- **The reference download now fetches the transcriptome** (`*.cdna.all.fa.gz`,

- **The reference download now fetches the transcriptome** (`*.cdna.all.fa.gz`,
  matched to the genome's assembly, never the `cdna.abinitio` prediction set).
  `--download_refs` previously produced no transcript FASTA at all, while
  `--build_indices` hard-requires one for Salmon and Kallisto — the two helper
  workflows could not feed each other. `DOWNLOAD_REFS` gained a matching
  `transcript_fasta` output.

  This is cDNA only. A total-RNA / rRNA-depleted library should also index the
  Ensembl ncRNA FASTA, or lncRNAs go unquantified.

- **Download integrity checking.** Every file is fetched over HTTPS with retries
  and backoff, checked against the advertised `Content-Length`, and streamed
  through a full gzip decode to catch truncation before the pipeline uses it.

- **Reference provenance.** `references/download_log.txt` now records the
  resolved Ensembl release, the assembly name, all three filenames, and the
  exact `--download_release <N>` needed to reproduce the set — so a `current`
  download stays interpretable after Ensembl moves on.

### Changed

- **`DOWNLOAD_REFS` output contract.** The `fasta` output glob is now anchored
  on `.dna.` so the soft-masked (`dna_sm`) and repeat-masked (`dna_rm`) genome
  variants can never be emitted in its place; `gtf` and `fasta` are no longer
  `optional`; and the `stub` block uses realistic Ensembl filenames so the
  output globs are actually exercised by `-stub-run`.

- **Genome FASTA selection** prefers a whole-genome `dna.primary_assembly` file
  and falls back to `dna.toplevel` only when the species publishes none — which
  correctly handles Drosophila, whose `primary_assembly` files are per-chromosome
  (`...dna.primary_assembly.2L.fa.gz`) rather than whole-genome.

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
