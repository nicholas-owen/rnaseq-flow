# rnaseq-flow Example Analysis

A complete, end-to-end worked example that takes **GSE151251** (human hepatic
stellate cells, ± TGF-β) from raw FASTQ files through HISAT2 alignment,
DESeq2 / edgeR differential expression, functional enrichment, and the run
reports. It runs in well under an hour on a modern workstation and produces a
textbook fibrosis-signature response to TGF-β, so it doubles as a verification
that your **rnaseq-flow** setup is working.

For the pipeline itself see [README.md](README.md); every parameter used below
is documented in [USAGE.md](USAGE.md), and every file it produces in
[OUTPUTS.md](OUTPUTS.md).

> **Why this dataset?** GSE151251 (BioProject **PRJNA635294**) is a small,
> modern, publicly available human bulk RNA-seq study with a clean
> treated-vs-control design (3 control + 3 TGF-β-treated hepatic stellate
> cells). TGF-β activation of HSCs drives an unmistakable upregulation of
> fibrosis / myofibroblast genes: COL1A1, COL3A1, ACTA2 (α-SMA), SERPINE1,
> IGFBP3. That gives you an unambiguous positive control to check against.

---

## 1. What you'll need

- **Nextflow 25.10.0+** (strict language; verified on 26.04) and **Java 17–26**.
- A **container engine**: Docker, Singularity/Apptainer, or Conda.
- Command-line tools: `curl`, `seqtk` for read subsampling, and `md5sum` for
  the download checks.
- **Disk space**: ~40 GB total. About 15 GB for the GRCh38 reference and
  HISAT2 index, ~10 GB for the FASTQs as downloaded, ~3 GB for the subsampled
  copies, the rest for results. The fetch script downloads full depth and
  subsamples afterwards, so the 10 GB is unavoidable even though the run only
  reads the 3 GB. Delete `hsc_data/raw/` once `hsc_data/subset/` exists if you
  are short of space and do not intend to try the splicing example.
- **Memory**: 16 GB RAM is comfortable. The HISAT2 index build is the heaviest
  step (~8 GB).

Install Nextflow if you don't already have it:

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version
```

Create a clean working directory and `cd` into it. All commands in this guide
assume you are at its top level and that `main.nf` is the **rnaseq-flow**
checkout (clone the repo separately or point at it with an absolute path):

```bash
mkdir rnaseq-flow_example && cd rnaseq-flow_example
```

---

## 2. Step 1: Fetch the reads and write the samplesheet

One script does all three jobs: it downloads the six runs, subsamples them,
and writes the samplesheet. `<rnaseq-flow>` below is your checkout.

```bash
<rnaseq-flow>/examples/human/fetch_hsc_tgfb_data.sh hsc_data 2000000
```

The first argument is the output directory, the second the number of read
pairs to keep per sample. Two million is plenty for differential expression
and leaves roughly 250 MB per file, about 3 GB for all six. Omit the second
argument and you get full depth, about 10 GB, which is what the splicing
example in [examples/human/README.md](examples/human/README.md) needs and this
one does not.

When it finishes you have:

```
hsc_data/raw/                 the downloaded FASTQs, checksum-verified
hsc_hsc_data/subset/              2 M pairs per sample, what the run below uses
hsc_data/samplesheet_hsc.csv  six rows, absolute paths, controls as REF
```

The dataset is BioProject **PRJNA635294** (GEO accession **GSE151251**), and
the six accessions are fixed in the script: `SRR11854548` to `SRR11854550` are
the controls, `SRR11854551` to `SRR11854553` the TGF-β-treated. FASTQ URLs are
resolved from ENA's `filereport` API at run time rather than hard-coded, so no
SRA Toolkit is needed and the script keeps working if ENA reorganises its
paths.

Three details worth knowing:

**Downloads are idempotent.** Each file lands as `.part`, has its MD5 checked
against ENA, and is renamed only on a match. An interrupted run resumes, and a
truncated file is never mistaken for a finished one. Re-running the script
costs nothing once the data is there.

**Subsampling needs `seqtk`.** It uses one fixed seed for both mates so pairs
stay in step. Without `seqtk` the script falls back to taking the first N
reads, says so loudly, and that is not a random sample: fine for checking the
pipeline runs, not for drawing conclusions.

**The controls are named `REF`.** That is the pipeline's baseline convention,
so a positive log2 fold change means up on TGF-β. The treated condition is
`TGFb`, which is why every contrast below reads `TGFb_vs_REF`.

The samplesheet validator aborts early if anything is wrong: missing columns,
missing FASTQ files, fewer than two conditions, or fewer than two replicates
per condition.

---

## 3. Step 2: Download references and build the HISAT2 index

These steps run **once** and are reused for every subsequent pipeline
invocation.

### 3.1 Download the GRCh38 genome, GTF and MSigDB gene sets

Pin a specific Ensembl release for reproducibility (release **110** used
here; any modern release works):

```bash
nextflow run main.nf \
    --download_refs \
    --download_species homo_sapiens \
    --download_source ensembl \
    --download_release 110 \
    --download_gmt \
    --organism hsapiens \
    --outdir references/human \
    -profile docker
```

This writes:

- `references/human/v110/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`: the genome FASTA.
- `references/human/v110/Homo_sapiens.GRCh38.110.gtf.gz`: the annotation.
- `references/human/v110/download_log.txt`: provenance (which release, when).
- `references/human/gmt/hallmark.gmt`, `c2_curated.gmt`, `c2_kegg.gmt`, `c2_reactome.gmt`, `c5_go.gmt`, `c5_go_bp.gmt`: MSigDB gene sets for GSEA.

### 3.2 Build the HISAT2 index

The HISAT2 build extracts splice sites and exons from the GTF automatically:

```bash
nextflow run main.nf \
    --build_indices \
    --aligner hisat2 \
    --genome_fasta references/human/v110/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    --gtf         references/human/v110/Homo_sapiens.GRCh38.110.gtf.gz \
    --outdir indices/human \
    -profile docker
```

The index lands in `indices/human/hisat2_index/`. This is the slow step
(typically 20–60 minutes and ~8 GB RAM) but you only do it once.

---

## 4. Step 3: Run the analysis

With the samplesheet, references and index in place, run the full RNASEQ
workflow with HISAT2:

```bash
nextflow run main.nf \
    --input hsc_data/samplesheet_hsc.csv \
    --aligner hisat2 \
    --hisat2_index indices/human/hisat2_index \
    --gtf references/human/v110/Homo_sapiens.GRCh38.110.gtf.gz \
    --strandedness auto \
    --gmt references/human/gmt/hallmark.gmt \
    --organism hsapiens \
    --read_length 100 \
    --outdir results \
    -profile docker
```

What happens, in order:

- **FastQC + fastp** QC and trim the subsampled reads.
- **HISAT2** aligns each sample to GRCh38; **samtools** indexes the BAMs.
- **RSeQC** infers per-sample strandedness (with `--strandedness auto` the
  verdict is fed straight into featureCounts), and checks read distribution
  and gene-body coverage.
- **deepTools** writes CPM-normalised BigWig coverage tracks.
- **featureCounts** produces the gene-level count matrix.
- **DESeq2** and **edgeR** run in parallel on those counts and emit one
  results CSV per contrast, here `TGFb_vs_REF`. DESeq2 log2 fold changes
  are apeglm-shrunken; each method also writes volcano, MA / smear and
  (DESeq2) PCA / (edgeR) MDS plots.
- **fgsea** (because you supplied `--gmt`) runs GSEA on the ranked DESeq2
  gene list, and **gprofiler2** tests the up- and down-regulated significant
  gene lists against GO terms and pathways.
- **rMATS** runs an alternative-splicing analysis on the BAMs (this is
  automatic on the HISAT2 route).
- **MultiQC** aggregates all the QC and the **Quarto analysis report**
  stitches everything (QC summary, plotly volcano plots, searchable DT result
  tables) into one HTML page.
- A **run-completion summary** (`results/pipeline_info/run_summary.html`) is
  written at the end with links to all of the above plus a per-process table
  of job time, peak memory and CPU.

> **Resume after an interruption.** Re-run the same command with `-resume`:
> Nextflow only re-executes the tasks affected by anything that changed.

---

## 5. Step 4: What was produced and how to read it

Start with the two reports:

- **`results/quarto_report/analysis_report.html`**: the headline interactive
  report. It opens with a Run overview (aligner, reference set and release,
  baseline) and the sample design. The DESeq2 and edgeR sections each show:
  significant-gene counts per contrast, an interactive plotly volcano (hover
  for gene name and stats), the embedded PCA / MDS / heatmap, and a searchable
  results table. A DESeq2-vs-edgeR agreement section lists how many genes were
  significant in both methods, and an Alternative splicing section summarises
  the rMATS events. It ends with a Files and locations index.
- **`results/multiqc/rnaseq-flow_multiqc_report.html`**: per-sample QC across
  FastQC, fastp, HISAT2 alignment, RSeQC and featureCounts. Use this to confirm
  alignment rates are healthy and strandedness was inferred consistently.

Then the result tables, all under `results/`:

| Path | What's in it |
|---|---|
| `deseq2_output/deseq2_results_TGFb_vs_REF.csv` | DESeq2 table: apeglm-shrunken log2FoldChange (and the unshrunken `log2FoldChange_MLE` beside it) + Wald stat / pvalue / padj, with `gene_id`, `gene_name` and `gene_biotype` columns |
| `edger_output/edger_results_TGFb_vs_REF.csv` | edgeR counterpart: logFC, logCPM, F-stat, PValue, FDR |
| `gsea_output/gsea_stats_TGFb_vs_REF.csv` + `gsea_plot_TGFb_vs_REF.svg` | Hallmark pathways enriched in TGF-β (positive NES = enriched in TGFb); every figure is also written as `.png` |
| `gprofiler_output/gprofiler_UP_TGFb_vs_REF.csv`, `gprofiler_DOWN_TGFb_vs_REF.csv` | GO terms / pathways over-represented in the up- and down-regulated gene lists, with their `gostplot_*` images |
| `*_output/reproduce/` | A script, data and objects that redraw each figure from the results alone ([OUTPUTS.md §18](OUTPUTS.md#18-redrawing-and-adapting-figures-reproduce)) |
| `featurecounts/`, `rseqc/`, `bam_to_bigwig/`, `hisat2/` | Upstream counts, QC and alignment outputs |
| `rmats_output/TGFb_vs_REF/` | Alternative-splicing events (SE / MXE / A3SS / A5SS / RI); `IncLevelDifference` is PSI(TGFb) minus PSI(REF), the same direction as the fold changes |
| `pipeline_info/run_manifest.json` | How this run was configured: aligner, references, baseline, design |
| `pipeline_info/run_summary.html` | Run status, duration, links to all outputs, per-process resource table |

---

## 6. Expected result: a built-in positive control

TGF-β activation of hepatic stellate cells drives an unmistakable, textbook
fibrosis signature. In `deseq2_results_TGFb_vs_REF.csv` you should see strong
upregulation (large positive `log2FoldChange`, very small `padj`) of:

- **COL1A1, COL3A1**: type I/III collagens.
- **ACTA2** (α-smooth-muscle actin): myofibroblast marker.
- **SERPINE1** (PAI-1): TGF-β's classical immediate target.
- **IGFBP3**: the dataset's own published top finding.

In the Hallmark GSEA output, the gene sets `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`
and `HALLMARK_TGF_BETA_SIGNALING` should be highly NES-positive. gProfiler
should likewise return extracellular-matrix / collagen / fibrosis GO terms
in the up-regulated list.

If you don't see this signature, the likeliest causes are: a wrong
samplesheet `condition` assignment (controls and treated swapped), a wrong
`--strandedness` (keep it at `auto` and check
`results/rseqc/<sample>.strandedness.txt` against your library type), or a
reference / index mismatch.

---

## 7. Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| `Samplesheet validation failed` | Re-check column names; `R1` paths must exist; conditions need ≥ 2 levels and ≥ 2 replicates each. |
| Very few reads `Assigned` by featureCounts | Strandedness mismatch. Check `results/rseqc/<sample>.strandedness.txt` and either leave `--strandedness auto` or set the matching value explicitly. |
| HISAT2 index build runs out of memory | Raise `--max_memory` (default `128.GB`; you need ~8 GB free for the build). |
| No `gsea_output/` | `--gmt` not supplied, or `--stop_at` was set before the enrichment stage. |
| Container pull failures | On HPC prefer `-profile singularity`; locally check Docker is running. |
| The DE table has very few significant genes | You subsampled too aggressively. Try increasing `N` (e.g. to 5 M pairs) and re-running. |

---

## 8. Reusing the example

The subsampled FASTQs in `hsc_data/subset/`, the references in
`references/human/v110/` and the HISAT2 index in
`indices/human/hisat2_index/` are all one-time setup. Once they exist,
re-running the pipeline against the same samplesheet takes only minutes,
which makes it easy to iterate on options, for example swapping to STAR with
`--aligner star --star_index ...`, or trying a different `--gmt` file
(`c2_kegg.gmt`, `c5_go_bp.gmt`).

For the isoform-level analyses, `--dtu`, `--diffsplice` and
`--isoform_switch`, go to [examples/human/README.md](examples/human/README.md).
It uses this same dataset and the same fetch script, so the download is already
done, but at full depth rather than subsampled: junction and isoform evidence
is a small fraction of gene-level counts, and 2 M pairs leaves too little to
test. Keep `hsc_data/raw/` and you can start there without downloading
anything again.

The two examples share one samplesheet, and the fetch script rewrites it each
run to point at the depth you asked for. If you go to the splicing example and
come back, re-run the command in section 2 first. Nothing is downloaded or
subsampled twice, so it takes seconds.

If instead of raw reads you already have a **gene count matrix** (from another
pipeline, or a public dataset such as a GEO supplementary file), you can skip the
FASTQ, alignment and quantification steps entirely and enter at differential
expression with `--counts`. You still supply the samplesheet (for `condition`)
and `--gtf` (for gene annotation); see
[USAGE.md §4.7](USAGE.md#47-starting-from-a-count-matrix---counts).

The whole example is deliberately small, reproducible (subsample seed,
pinned Ensembl release) and self-verifying through the TGF-β signature,
so it serves as both a tutorial and an ongoing setup check.
