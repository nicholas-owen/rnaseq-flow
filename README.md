# rnaseq-flow

A modular **bulk RNA-seq analysis pipeline** built with Nextflow DSL2. It takes
raw FASTQ files all the way through quality control, alignment/quantification,
differential expression, alternative splicing, fusion detection and functional
enrichment, and aggregates everything into interactive HTML reports.

The pipeline is *modular*: you choose the aligner, choose which analysis stages
run, and supply references either directly or via built-in helper workflows
that download and index them for you.

---

## Contents

- [Features](#features)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [The three run modes](#the-three-run-modes)
- [Pipeline stages](#pipeline-stages)
- [Samplesheet format](#samplesheet-format)
- [Parameter reference](#parameter-reference)
- [Choosing an aligner](#choosing-an-aligner)
- [Documentation](#documentation)
- [Known limitations](#known-limitations)

For step-by-step run instructions see **[USAGE.md](USAGE.md)**.
For how to read every output file see **[OUTPUTS.md](OUTPUTS.md)**.
For a visual map of the workflow, open **[overview.html](overview.html)** in a
browser — it shows the pipeline as both a DAG and an nf-core-style metro map.

---

## Features

| Stage | Tools | Notes |
|---|---|---|
| Read QC | FastQC, fastp | Raw QC + adapter/quality trimming |
| Genome alignment | STAR, HISAT2 | Produces sorted, indexed BAMs |
| Pseudo-alignment | Salmon, Kallisto | Fast transcript-level quantification |
| Alignment QC | RSeQC | Strandedness (auto-fed to featureCounts), read distribution, gene-body coverage |
| Coverage tracks | deepTools | CPM-normalised BigWig files for IGV/UCSC |
| Gene quantification | featureCounts | Gene-level count matrix (STAR/HISAT2) |
| Transcript&rarr;gene counts | tximport | Summarises Salmon/Kallisto to gene level for DE |
| Differential expression | DESeq2, edgeR | Two independent callers, run in parallel on the same gene counts |
| Alternative splicing | rMATS | SE, MXE, A3SS, A5SS, RI events |
| Fusion detection | STAR-Fusion | Requires a CTAT genome library |
| Isoform switching | IsoformSwitchAnalyzeR | Salmon transcript-level analysis |
| Differential transcript usage | DEXSeq | Transcript-usage changes (Salmon/Kallisto) |
| Differential splicing | edgeR `diffSpliceDGE` | Exon-usage (STAR/HISAT2) and transcript-usage (Salmon/Kallisto) shifts vs the gene |
| Functional enrichment | fgsea, gprofiler2 | GSEA + GO/pathway over-representation, from the DESeq2 results |
| Reporting | MultiQC, Quarto | Aggregated HTML reports |

Helper workflows download reference genomes/annotation from Ensembl and build
aligner indices, so you can go from nothing to results with three commands.

---

## Requirements

- **Nextflow** `>=22.10.1` (`curl -s https://get.nextflow.io | bash`)
- **Java** 11–21 (required by Nextflow)
- A **container engine** — Docker *or* Singularity/Apptainer — or **Conda**.
  Every process declares its own container, so nothing else needs installing.

Run with one of the bundled profiles: `-profile docker`, `-profile singularity`
or `-profile conda`.

---

## Quick start

```bash
# 1. Download a reference genome + annotation, pinned to an Ensembl release
nextflow run main.nf \
    --download_refs \
    --download_species homo_sapiens \
    --download_release 102 \
    --download_gmt \
    --organism hsapiens \
    --outdir references/human \
    -profile docker
# genome FASTA + GTF land in references/human/v102/

# 2. Build the index for your chosen aligner
nextflow run main.nf \
    --build_indices \
    --aligner star \
    --genome_fasta references/human/v102/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    --gtf         references/human/v102/Homo_sapiens.GRCh38.102.gtf.gz \
    --outdir indices/human \
    -profile docker

# 3. Run the analysis
nextflow run main.nf \
    --input samplesheet.csv \
    --aligner star \
    --star_index indices/human/star_index \
    --gtf references/human/Homo_sapiens.GRCh38.110.gtf.gz \
    --strandedness reverse \
    --outdir results \
    -profile docker
```

Reference FASTA/GTF files may be supplied **gzipped or uncompressed** — the
pipeline decompresses them internally where a tool requires it.

---

## The three run modes

`main.nf` has a single entry point that dispatches to one of three workflows
based on the flags you pass:

| Flag | Workflow | Purpose |
|---|---|---|
| `--download_refs` | `DOWNLOAD` | Fetch genome FASTA + GTF (and optionally GMT gene sets) from Ensembl |
| `--build_indices` | `BUILD_INDICES` | Build STAR / HISAT2 / Salmon / Kallisto indices from references |
| *(neither; `--input` given)* | `RNASEQ` | The full analysis pipeline |

Only one mode runs per invocation.

---

## Pipeline stages

The main `RNASEQ` workflow runs these stages in order. Use `--stop_at` to halt
early (see [Parameter reference](#parameter-reference)).

1. **Read QC & trimming** — FastQC on raw reads; fastp adapter/quality trimming.
2. **Alignment / quantification** — STAR or HISAT2 (genome) *or* Salmon or
   Kallisto (pseudo-alignment).
3. **Post-alignment QC** *(genome aligners only)* — BAM indexing, RSeQC, BigWig
   coverage tracks.
4. **Gene quantification** — featureCounts gene matrix (STAR/HISAT2), or
   tximport summarising Salmon/Kallisto transcript quantification to gene level.
5. **Differential expression** — DESeq2 and edgeR. Both run, in parallel, on the
   same gene counts as two independent callers; neither feeds the other. DESeq2
   log2 fold changes are apeglm-shrunken (`lfcShrink`) for better gene ranking.
6. **Alternative splicing** — rMATS (genome aligners, CSV samplesheet).
7. **Fusion detection** — STAR-Fusion (STAR only, needs `--ctat_lib`).
8. **Isoform switching** — IsoformSwitchAnalyzeR (Salmon only, needs
   `--isoform_switch` and `--transcript_fasta`).
9. **Differential transcript usage** — DEXSeq (Salmon/Kallisto, opt-in with
   `--dtu`): tests whether transcript-isoform proportions shift between
   conditions, complementing gene-level DE.
10. **Differential splicing** — edgeR `diffSpliceDGE` (opt-in with
    `--diffsplice`): an exon-usage test on STAR/HISAT2 and a transcript-usage
    test on Salmon/Kallisto, each comparing a feature's fold change against its
    gene's overall fold change.
11. **Functional enrichment** — fgsea (needs `--gmt`) and gprofiler2, run on the
    **DESeq2** results only (edgeR results are not used downstream).
12. **Reporting** — MultiQC and an interactive Quarto analysis report (QC
    summary plus DE tables, plotly volcano plots and enrichment summaries).

Every differential-expression result table (DESeq2, edgeR, DTU and diffSplice)
is annotated with `gene_name` and `gene_biotype` columns parsed from the GTF, so
the outputs are readable without a separate gene-ID lookup.

---

## Samplesheet format

`--input` takes a CSV file with at least these columns:

```csv
sample,R1,R2,condition
ctrl_rep1,data/ctrl_rep1_R1.fastq.gz,data/ctrl_rep1_R2.fastq.gz,REF
ctrl_rep2,data/ctrl_rep2_R1.fastq.gz,data/ctrl_rep2_R2.fastq.gz,REF
treat_rep1,data/treat_rep1_R1.fastq.gz,data/treat_rep1_R2.fastq.gz,treatment
treat_rep2,data/treat_rep2_R1.fastq.gz,data/treat_rep2_R2.fastq.gz,treatment
```

| Column | Description |
|---|---|
| `sample` | Unique sample identifier (used to name all outputs) |
| `R1` | Path to the R1 (or single-end) FASTQ; absolute or relative to the launch directory |
| `R2` | Path to the R2 FASTQ. **Leave empty for single-end data** |
| `condition` | Experimental group — drives the differential-expression design |
| `batch` *(optional)* | Batch / covariate label; if present, the DE model becomes `~ batch + condition` automatically |

> **Important — the `REF` convention.** The differential-expression scripts
> treat a condition literally named `REF` as the baseline/denominator of every
> contrast, so a positive log2 fold change means "up relative to `REF`". Name
> your control group `REF` to get correctly-oriented results. If no group is
> called `REF`, conditions are ordered alphabetically and the comparison
> direction is arbitrary.

Single-end and paired-end samples are detected automatically from whether `R2`
is present.

> **Validation.** The samplesheet is checked before any work starts, and the
> run aborts immediately if there is a problem. The checks: required columns
> present, sample ids unique, every R1/R2 FASTQ exists, and — because DESeq2 and
> edgeR need them — **at least 2 conditions and at least 2 replicates per
> condition**. (With `--stop_at preQC`/`postQC`, where no differential
> expression runs, the condition/replicate rules are downgraded to warnings.)

> **Batch effects & covariates.** Extra columns beyond the four above are kept
> and may be used as model covariates. A `batch` column is auto-detected and
> turns the gene-level DESeq2/edgeR model into `~ batch + condition`. For other
> confounders, add the column(s) and pass an explicit formula with `--design`
> (e.g. `--design "~ sex + batch + condition"`); `condition` must always be
> included, as it is the term being contrasted.

---

## Parameter reference

Run `nextflow run main.nf --help` to print every parameter, grouped, with its
default. All parameters are declared in `nextflow_schema.json`; an unrecognised
`--parameter` (e.g. a typo) **aborts the run immediately** rather than being
silently ignored.

### Input / output

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--input` | Yes (RNASEQ mode) | — | Path to the CSV samplesheet |
| `--outdir` | No | `results` | Output directory |
| `--publish_dir_mode` | No | `copy` | How results are published (`copy`, `symlink`, `link`) |
| `--stop_at` | No | — | Stop after a stage: `preQC`, `postQC`, `DE`, `GSEA` |

### Alignment & references

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--aligner` | No | `star` | `star`, `hisat2`, `salmon` or `kallisto` |
| `--strandedness` | No | `auto` | `auto` (RSeQC infers it for STAR/HISAT2), `unstranded`, `forward`, `reverse` — `forward`/`reverse` also set kallisto `--fr-stranded`/`--rf-stranded` |
| `--gtf` | Yes* | — | Gene annotation (GTF, may be gzipped) |
| `--star_index` | Yes* | — | STAR index directory (if `--aligner star`) |
| `--hisat2_index` | Yes* | — | HISAT2 index directory (if `--aligner hisat2`) |
| `--salmon_index` | Yes* | — | Salmon index directory (if `--aligner salmon`) |
| `--kallisto_index` | Yes* | — | Kallisto index file (if `--aligner kallisto`) |
| `--transcript_fasta` | Yes* | — | Transcript FASTA (isoform switching) |

### Analysis options

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--read_length` | No | `100` | Read length passed to rMATS |
| `--ctat_lib` | No | — | CTAT genome library directory — enables STAR-Fusion |
| `--isoform_switch` | No | `false` | Enable IsoformSwitchAnalyzeR (needs `--aligner salmon`) |
| `--dtu` | No | `false` | Enable DEXSeq differential transcript usage (needs `--aligner salmon`/`kallisto` and `--gtf`) |
| `--diffsplice` | No | `false` | Enable edgeR `diffSpliceDGE` — exon-level usage (STAR/HISAT2) or transcript-level usage (Salmon/Kallisto); needs `--gtf` |
| `--design` | No | auto | DESeq2/edgeR model formula, e.g. `~ batch + condition`. Default: `~ batch + condition` if the samplesheet has a `batch` column, else `~ condition` |
| `--gmt` | No | — | GMT gene-set file — enables GSEA |
| `--organism` | No | `hsapiens` | Organism ID for gProfiler / GMT download |

### Helper workflows

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--download_refs` | No | `false` | Run the reference-download workflow |
| `--download_species` | Yes* | `homo_sapiens` | Ensembl species name |
| `--download_source` | No | `ensembl` | Reference source |
| `--download_release` | No | `current` | Pin an Ensembl release (e.g. `102`) for a reproducible download; output goes to `<outdir>/v<release>` |
| `--download_gmt` | No | `false` | Also download MSigDB GMT gene sets |
| `--build_indices` | No | `false` | Run the index-building workflow |
| `--genome_fasta` | Yes* | — | Genome FASTA for index building |

### Resources

| Parameter | Default | Description |
|---|---|---|
| `--max_cpus` | `16` | Per-process CPU ceiling |
| `--max_memory` | `128.GB` | Per-process memory ceiling |
| `--max_time` | `240.h` | Per-process wall-time ceiling |

\* "Yes*" = required only when the relevant tool or mode is used.

---

## Choosing an aligner

| Aligner | Type | Best for | Downstream stages available |
|---|---|---|---|
| **STAR** | Genome | The default. Full analysis incl. fusions & splicing | QC, BigWig, featureCounts, **DESeq2/edgeR**, rMATS, STAR-Fusion, GSEA, gProfiler, diffSplice (exon) |
| **HISAT2** | Genome | Lighter-weight genome alignment | QC, BigWig, featureCounts, **DESeq2/edgeR**, rMATS, GSEA, gProfiler, diffSplice (exon) |
| **Salmon** | Pseudo | Fast transcript quantification, isoform switching | Transcript quant, **DESeq2/edgeR** (via tximport), GSEA, gProfiler, IsoformSwitchAnalyzeR, DEXSeq DTU, diffSplice (transcript) |
| **Kallisto** | Pseudo | Fastest transcript quantification | Transcript quant, **DESeq2/edgeR** (via tximport), GSEA, gProfiler, DEXSeq DTU, diffSplice (transcript) |

> **Differential expression** runs with **every aligner**. STAR/HISAT2 feed
> DESeq2/edgeR from the featureCounts gene matrix; Salmon/Kallisto feed them via
> **tximport**, which summarises transcript-level quantification to gene level
> (DESeq2 uses the transcript-length offset for accurate normalisation). The
> tximport path needs `--gtf` for the transcript-to-gene map.

---

## Documentation

| File | Contents |
|---|---|
| [USAGE.md](USAGE.md) | Detailed run instructions, worked scenarios, troubleshooting |
| [OUTPUTS.md](OUTPUTS.md) | Every output directory and how to interpret each file |
| [overview.html](overview.html) | Interactive workflow diagram — DAG and metro-map views (open in a browser) |
| [CHANGELOG.md](CHANGELOG.md) | Release history |
| [CITATIONS.md](CITATIONS.md) | Every tool with its publication, plus a methods-paragraph template |
| [CITATIONS.html](CITATIONS.html) | The citations page as a styled HTML report (open in a browser) |

---

## Known limitations

- Differential expression is **gene-level** by default (DESeq2/edgeR). A
  transcript-level differential transcript usage (DTU) test using DEXSeq is
  available opt-in via `--dtu` for the Salmon/Kallisto routes.
- The tximport DE path requires `--gtf` so a transcript-to-gene map can be
  built; without it, DE is skipped for Salmon/Kallisto.
- The Ensembl reference downloader scrapes the Ensembl FTP listing; for
  unusual species, download references manually and pass them directly.
- `--ctat_lib` must point to an already-extracted CTAT genome library
  directory (the pipeline does not untar archives).

---

## Pipeline configuration files

| File | Purpose |
|---|---|
| `nextflow.config` | Parameters, profiles, manifest, `check_max` resource function |
| `conf/base.config` | Default CPU/memory/time per process label |
| `conf/modules.config` | `publishDir` rules — controls what lands in `--outdir` |
| `nextflow_schema.json` | Parameter schema — powers `--help`, typo detection, and nf-core tooling |
| `assets/multiqc_config.yml` | MultiQC report config — title, module order, sample-name cleaning, and a custom CSS/logo theme |
| `assets/analysis_report.qmd` | Quarto template for the interactive analysis report (QC summary, DE tables and plotly volcano plots, enrichment) |

---

## License

rnaseq-flow is released under the **GNU General Public License v3.0**
(GPL-3.0). See the [LICENSE](LICENSE) file for the full licence text.
