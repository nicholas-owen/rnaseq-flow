# rnaseq-flow — Usage Guide

This guide walks through installing the pipeline, preparing references, and
running each of the three workflows. For output interpretation see
[OUTPUTS.md](OUTPUTS.md); for an overview see [README.md](README.md).

---

## 1. Installation

```bash
# Install Nextflow (needs Java 11-21)
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/      # or add it to your PATH

nextflow -version                     # confirm it runs
```

You also need **Docker**, **Singularity/Apptainer** or **Conda**. You do not
need to install any bioinformatics tools yourself — each process pulls its own
container.

Get the pipeline:

```bash
git clone <your-repo-url> rnaseq-flow
cd rnaseq-flow
```

### Validate the install without real data

Every process ships a `stub` block, so you can dry-run the whole DAG in seconds
without containers or data:

```bash
nextflow run main.nf --input samplesheet_test.csv --aligner star \
    --star_index fake --gtf fake.gtf -stub-run
```

`-stub-run` executes the lightweight stub of each process (it just creates
placeholder output files), which confirms the channel wiring is sound.

### Discovering parameters

```bash
nextflow run main.nf --help
```

prints every parameter — grouped, with defaults — from `nextflow_schema.json`.
Parameters are also validated on every run: an unrecognised `--option` (a typo)
aborts immediately with a clear message, instead of being silently ignored.

---

## 2. Mode 1 — Download references

Fetches the genome FASTA and GTF annotation for any Ensembl species, and
(optionally) MSigDB gene sets for GSEA.

```bash
nextflow run main.nf \
    --download_refs \
    --download_species homo_sapiens \
    --download_source ensembl \
    --download_release 102 \
    --download_gmt \
    --organism hsapiens \
    --outdir references/human \
    -profile docker
```

| Parameter | Description |
|---|---|
| `--download_species` | Ensembl species, lower-case with underscores (`homo_sapiens`, `mus_musculus`) |
| `--download_release` | Ensembl release to pin, e.g. `102` (a leading `v` is accepted). Omit it for the rolling `current` release. **Pin a release for reproducible references.** |
| `--download_gmt` | Also download Hallmark / C2 / C5 gene sets |
| `--organism` | gProfiler-style organism ID; mapped to a scientific name for MSigDB |

**Reproducibility — the versioned subfolder.** When `--download_release` is
set, the genome FASTA + GTF are written to a `v<release>` subfolder *under*
`--outdir`, so each Ensembl build is kept separate. With
`--outdir references/human --download_release 102` you get:

```
references/human/
├── v102/                                             # pinned Ensembl release 102
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz   # genome FASTA
│   ├── Homo_sapiens.GRCh38.102.gtf.gz                   # annotation (note the .102.)
│   └── download_log.txt                                 # provenance log (records the release)
└── gmt/                                              # if --download_gmt
    ├── hallmark.gmt
    ├── c2_curated.gmt   c2_kegg.gmt   c2_reactome.gmt
    └── c5_go.gmt        c5_go_bp.gmt
```

Without `--download_release`, the latest release is downloaded straight into
`--outdir` (no `v<release>` subfolder). The `download_log.txt` always records
the exact release used.

> For a transcript FASTA (needed for Salmon/Kallisto or isoform switching),
> download the Ensembl `cdna.all.fa.gz` file for your species separately.

---

## 3. Mode 2 — Build indices

Builds the index for one aligner, or all four at once with `--aligner all`.

```bash
# STAR (needs genome FASTA + GTF)
nextflow run main.nf --build_indices --aligner star \
    --genome_fasta references/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    --gtf         references/human/Homo_sapiens.GRCh38.110.gtf.gz \
    --outdir indices/human -profile docker

# HISAT2 (genome FASTA + GTF; splice sites and exons are extracted from the GTF)
nextflow run main.nf --build_indices --aligner hisat2 \
    --genome_fasta references/human/genome.fa.gz \
    --gtf         references/human/annotation.gtf.gz \
    --outdir indices/human -profile docker

# Salmon (decoy-aware: needs transcript FASTA + genome FASTA)
nextflow run main.nf --build_indices --aligner salmon \
    --transcript_fasta references/human/transcripts.cdna.all.fa.gz \
    --genome_fasta     references/human/genome.fa.gz \
    --outdir indices/human -profile docker

# Kallisto (transcript FASTA only)
nextflow run main.nf --build_indices --aligner kallisto \
    --transcript_fasta references/human/transcripts.cdna.all.fa.gz \
    --outdir indices/human -profile docker
```

| Aligner | Requires | Produces |
|---|---|---|
| `star` | `--genome_fasta`, `--gtf` | `star_index/` |
| `hisat2` | `--genome_fasta`, `--gtf` | `hisat2_index/` |
| `salmon` | `--transcript_fasta`, `--genome_fasta` | `salmon_index/` |
| `kallisto` | `--transcript_fasta` | `kallisto_index` |
| `all` | all of the above | all four |

Indices are written under `--outdir`. FASTA/GTF inputs may be gzipped.

---

## 4. Mode 3 — Run the analysis

### 4.1 Prepare the samplesheet

```csv
sample,R1,R2,condition
ctrl_1,/data/ctrl_1_R1.fastq.gz,/data/ctrl_1_R2.fastq.gz,REF
ctrl_2,/data/ctrl_2_R1.fastq.gz,/data/ctrl_2_R2.fastq.gz,REF
treat_1,/data/treat_1_R1.fastq.gz,/data/treat_1_R2.fastq.gz,treatment
treat_2,/data/treat_2_R1.fastq.gz,/data/treat_2_R2.fastq.gz,treatment
```

- Name the control/baseline group **`REF`** so fold changes are oriented as
  "treatment vs control" (see README → Samplesheet format).
- Leave `R2` empty for single-end reads.
- Provide **at least 2 conditions and at least 2 replicates per condition** —
  this is enforced (DESeq2/edgeR cannot estimate dispersion otherwise).
- Optionally add a **`batch`** column (or other covariate columns) to model
  confounders in the differential-expression design — see §4.5.

The samplesheet is validated before the run starts: if a required column is
missing, a sample id is duplicated, a FASTQ path does not exist, or the
condition/replicate minimums are not met, the run aborts immediately with a
list of every problem found. (When `--stop_at preQC`/`postQC` skips
differential expression, the condition/replicate rules become warnings.)

### 4.2 Worked scenarios

**Scenario A — comprehensive STAR run** (fusions + splicing + DE + enrichment):

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --aligner star \
    --star_index indices/human/star_index \
    --gtf references/human/Homo_sapiens.GRCh38.110.gtf.gz \
    --strandedness reverse \
    --ctat_lib /refs/GRCh38_gencode_v37_CTAT_lib/ctat_genome_lib_build_dir \
    --read_length 150 \
    --diffsplice \
    --gmt references/human/gmt/c2_kegg.gmt \
    --organism hsapiens \
    --outdir results \
    -profile docker
```

**Scenario B — Salmon transcript-level analysis (isoform switching + DTU + diffSplice):**

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --aligner salmon \
    --salmon_index indices/human/salmon_index \
    --gtf references/human/Homo_sapiens.GRCh38.110.gtf.gz \
    --transcript_fasta references/human/transcripts.cdna.all.fa.gz \
    --isoform_switch \
    --dtu \
    --diffsplice \
    --outdir results \
    -profile docker
```

> When `--gtf` is supplied, Salmon/Kallisto transcript quantification is
> summarised to gene level by **tximport** and run through DESeq2/edgeR and the
> enrichment stages — the same downstream analysis as the genome aligners.
> Salmon additionally supports isoform-switch analysis. `--dtu` adds a DEXSeq
> **differential transcript usage** test (Salmon or Kallisto) that flags genes
> whose isoform proportions shift between conditions — it is opt-in because
> DEXSeq is compute-heavy. `--diffsplice` adds edgeR's `diffSpliceDGE` test:
> transcript-level usage on the Salmon/Kallisto route, exon-level usage on the
> STAR/HISAT2 route (it also needs `--gtf`). Without `--gtf`, DE, DTU and
> diffSplice are skipped (no transcript-to-gene map can be built). Every
> DE / DTU / diffSplice result table carries `gene_name` and `gene_biotype`
> columns parsed from the GTF, so the outputs are readable at a glance.

**Scenario C — standard HISAT2 run:**

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --aligner hisat2 \
    --hisat2_index indices/human/hisat2_index \
    --gtf references/human/Homo_sapiens.GRCh38.110.gtf.gz \
    --strandedness reverse \
    --outdir results \
    -profile docker
```

### 4.3 Strandedness

`--strandedness` sets the library strand setting for featureCounts (the genome
aligners) and for kallisto.

| Setting | Behaviour |
|---|---|
| `auto` *(default)* | STAR/HISAT2: RSeQC `infer_experiment` measures the strandedness of every sample and feeds the verdict straight into featureCounts. Kallisto: no BAM exists to inspect, so kallisto runs library-type-agnostic (see below). |
| `unstranded` | Force unstranded (`featureCounts -s 0`; no kallisto strand flag). |
| `reverse` | Force reverse-stranded (`featureCounts -s 2`; `kallisto --rf-stranded`) — Illumina TruSeq Stranded mRNA and most dUTP kits. |
| `forward` | Force forward-stranded (`featureCounts -s 1`; `kallisto --fr-stranded`) — ligation-based / some older kits. |

With the default `auto`, the pipeline detects strandedness per sample for the
genome aligners (STAR/HISAT2): each sample's verdict is written to
`results/rseqc/<sample>.strandedness.txt` and applied to its featureCounts run.
Set `--strandedness` explicitly only to override the inference.

**Salmon** auto-detects library type internally (`--libType A`), so it needs no
setting. **Kallisto** has no equivalent and produces no BAM for RSeQC, so it
*cannot* auto-detect: with `auto` (or `unstranded`) it runs without a strand
flag. If your library is stranded and you use kallisto, pass `--strandedness
forward` or `reverse` explicitly so the `--fr-stranded` / `--rf-stranded` flag
is applied (the pipeline logs a warning when kallisto runs with `auto`).

### 4.4 Staging with `--stop_at`

| Value | Pipeline stops after |
|---|---|
| `preQC` | FastQC + fastp trimming |
| `postQC` | Alignment, BAM indexing, RSeQC, BigWig |
| `DE` | featureCounts + DESeq2 + edgeR |
| `GSEA` *(or unset)* | Everything, incl. enrichment, fusions, splicing |

```bash
nextflow run main.nf --input samplesheet.csv --aligner star \
    --star_index indices/human/star_index \
    --gtf references/human/annotation.gtf --stop_at postQC -profile docker
```

### 4.5 Batch effects & covariates

The gene-level differential-expression model (DESeq2 and edgeR) defaults to
`~ condition`. To adjust for a confounder there are two ways:

- **Quick way — a `batch` column.** Add a `batch` column to the samplesheet;
  the pipeline detects it and the model automatically becomes
  `~ batch + condition` for both DESeq2 and edgeR.
- **Full control — `--design`.** Pass an explicit model formula, e.g.
  `--design "~ sex + batch + condition"`. Every variable in the formula must be
  a samplesheet column, and `condition` must be included — it stays the
  variable contrasted, so the result tables and the downstream GSEA / gProfiler
  steps are unchanged. An explicit `--design` overrides the automatic `batch`
  behaviour.

```csv
sample,R1,R2,condition,batch
ctrl_1,/data/ctrl_1_R1.fastq.gz,/data/ctrl_1_R2.fastq.gz,REF,b1
ctrl_2,/data/ctrl_2_R1.fastq.gz,/data/ctrl_2_R2.fastq.gz,REF,b2
treat_1,/data/treat_1_R1.fastq.gz,/data/treat_1_R2.fastq.gz,treatment,b1
treat_2,/data/treat_2_R1.fastq.gz,/data/treat_2_R2.fastq.gz,treatment,b2
```

Covariates are modelled as categorical factors. Avoid a covariate that is
fully confounded with `condition` (e.g. every `REF` sample in one batch and
every treated sample in another): the model would be rank-deficient and
DESeq2/edgeR would abort. `--design` covers the gene-level DE tests; rMATS, DTU
and diffSplice keep their own designs.

---

## 5. Profiles & execution

| Profile | Effect |
|---|---|
| `-profile docker` | Run every process in its Docker container |
| `-profile singularity` | Use Singularity/Apptainer (HPC-friendly) |
| `-profile conda` | Create per-process Conda environments |

Add `-resume` to continue from cached results after an interruption or a
parameter tweak:

```bash
nextflow run main.nf --input samplesheet.csv ... -profile docker -resume
```

To submit to an HPC scheduler, add an executor block (SLURM example) to a
custom config and pass it with `-c`:

```groovy
process.executor = 'slurm'
process.queue    = 'normal'
```

---

## 6. Resource tuning

Per-process CPU/memory/time come from `conf/base.config`, keyed off process
labels (`process_low/medium/high`). They are automatically capped by:

```bash
--max_cpus 32 --max_memory 256.GB --max_time 120.h
```

A failed task is retried **once** with doubled memory/time. STAR genome
indexing of a mammalian genome needs ~38 GB RAM — make sure `--max_memory`
allows it.

---

## 7. Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| `Please specify an input samplesheet ... with --input` | No `--input`, `--download_refs` or `--build_indices` given |
| `STAR index not provided via --star_index` | Pass the matching `--<aligner>_index` for your `--aligner` |
| `Cannot find any reads matching` | FASTQ paths in the samplesheet are wrong or relative to the wrong directory |
| STAR-Fusion produced nothing | `--ctat_lib` not set, or not pointing at the extracted `ctat_genome_lib_build_dir` |
| No `deseq2_output/` | Needs a CSV samplesheet; for Salmon/Kallisto also needs `--gtf` (the tximport transcript-to-gene map) |
| No `diffsplice_output/` | `--diffsplice` not set, or `--gtf` missing (needed for exon counting and the transcript-to-gene map) |
| Splicing/fusion missing | These run only at full depth — don't combine with `--stop_at DE` |
| Out-of-memory kills | Raise `--max_memory`, or lower concurrency |
| Container pull failures | Check internet access / Docker daemon; on HPC prefer `-profile singularity` |

Execution reports are written to `<outdir>/pipeline_info/` (timeline, trace,
resource report) — inspect these to find which process failed and why.

To validate pipeline structure after editing it, re-run with `-stub-run`
(see §1).
