# rnaseq-flow — Example Analysis

A complete, end-to-end worked example that takes **GSE151251** (human hepatic
stellate cells, ± TGF-β) from raw FASTQ files through HISAT2 alignment,
DESeq2 / edgeR differential expression, functional enrichment, and the run
reports. It runs in well under an hour on a modern workstation and produces a
textbook fibrosis-signature response to TGF-β, so it doubles as a verification
that your **rnaseq-flow** setup is working.

> **Why this dataset?** GSE151251 (BioProject **PRJNA635294**) is a small,
> modern, publicly available human bulk RNA-seq study with a clean
> treated-vs-control design (3 control + 3 TGF-β-treated hepatic stellate
> cells). TGF-β activation of HSCs drives an unmistakable upregulation of
> fibrosis / myofibroblast genes — COL1A1, COL3A1, ACTA2 (α-SMA), SERPINE1,
> IGFBP3 — which gives you an unambiguous positive control to check against.

---

## 1. What you'll need

- **Nextflow 22.10.1+** and **Java 11–21**.
- A **container engine**: Docker, Singularity/Apptainer, or Conda.
- Command-line tools: `curl` (or `wget`), `seqtk` (for read subsampling),
  `column` and `awk` (for reading the run table).
- **Disk space**: ~30 GB total — about 15 GB for the GRCh38 reference and
  HISAT2 index, ~2 GB for the subsampled FASTQs, the rest for results.
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

## 2. Step 1 — Fetch the FASTQ files

The dataset is BioProject **PRJNA635294** (GEO accession **GSE151251**). The
easiest way to list and download its FASTQs is via ENA's `filereport` API —
no SRA Toolkit needed, FASTQs come straight from the ENA HTTP/FTP server.

Fetch the run table:

```bash
mkdir -p data
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA635294&result=read_run&fields=run_accession,sample_title,library_layout,fastq_ftp,fastq_md5,read_count&format=tsv" \
    > data/runinfo.tsv

column -t -s $'\t' data/runinfo.tsv | less -S    # browse the table
```

The `sample_title` column identifies each sample. For this example you want
the **six samples that are TGF-β-treated vs control without any other
perturbation** — three controls and three TGF-β-treated. Note the
`run_accession` of each.

Once you've picked the six runs, download their FASTQs from the URLs in the
`fastq_ftp` column. A small loop will do it (the `fastq_ftp` cell contains
one URL for single-end or two semicolon-separated URLs for paired-end):

```bash
# Put the six run_accession values you picked here, controls first.
SELECTED=(SRRxxxxxxx1 SRRxxxxxxx2 SRRxxxxxxx3 SRRxxxxxxx4 SRRxxxxxxx5 SRRxxxxxxx6)

mkdir -p data/raw
for run in "${SELECTED[@]}"; do
    urls=$(awk -F'\t' -v r="$run" 'NR>1 && $1==r {print $4}' data/runinfo.tsv | tr ';' '\n')
    for url in $urls; do
        echo "Downloading $url"
        curl -L -o "data/raw/$(basename "$url")" "ftp://${url}"
    done
done
```

> **Alternative.** If you prefer the SRA Toolkit, the equivalent is
> `prefetch <SRR>` followed by `fasterq-dump --split-files <SRR>`. ENA's
> direct FASTQ downloads are usually faster and need no toolkit setup.

---

## 3. Step 2 — Subset the reads to a reusable size

Modern RNA-seq depths (20–40 M reads / sample) are more than enough for a
worked example. Subsample each FASTQ to **2 million read pairs** with `seqtk`,
using the **same fixed seed for R1 and R2** of each sample so the pairs stay
in sync:

```bash
mkdir -p data/subset
SEED=42
N=2000000

ls data/raw/*_1.fastq.gz | sed 's|.*/||; s|_1\.fastq\.gz$||' | while read run; do
    echo "Subsampling $run"
    seqtk sample -s $SEED "data/raw/${run}_1.fastq.gz" $N | gzip > "data/subset/${run}_1.fastq.gz"
    if [ -f "data/raw/${run}_2.fastq.gz" ]; then
        seqtk sample -s $SEED "data/raw/${run}_2.fastq.gz" $N | gzip > "data/subset/${run}_2.fastq.gz"
    fi
done
```

`data/subset/` now contains the reusable example data — roughly ~250 MB per
file (so ~3 GB total for six paired-end samples). You can keep these and
re-run the pipeline against them as many times as you like.

---

## 4. Step 3 — Write the samplesheet

rnaseq-flow's samplesheet is a CSV with one row per sample. Naming the
controls `REF` makes positive log2 fold changes mean "up on TGF-β" (the
pipeline's REF convention). Create `samplesheet.csv` next to `data/`:

```csv
sample,R1,R2,condition
hsc_ctrl_1,data/subset/<run1>_1.fastq.gz,data/subset/<run1>_2.fastq.gz,REF
hsc_ctrl_2,data/subset/<run2>_1.fastq.gz,data/subset/<run2>_2.fastq.gz,REF
hsc_ctrl_3,data/subset/<run3>_1.fastq.gz,data/subset/<run3>_2.fastq.gz,REF
hsc_tgfb_1,data/subset/<run4>_1.fastq.gz,data/subset/<run4>_2.fastq.gz,tgfb
hsc_tgfb_2,data/subset/<run5>_1.fastq.gz,data/subset/<run5>_2.fastq.gz,tgfb
hsc_tgfb_3,data/subset/<run6>_1.fastq.gz,data/subset/<run6>_2.fastq.gz,tgfb
```

Substitute your actual run accessions for `<run1>…<run6>`. Use the
`sample_title` column in `data/runinfo.tsv` to confirm which three are
controls and which three are TGF-β-treated. If the runs are single-end (no
`_2.fastq.gz`), leave the `R2` column empty — the pipeline auto-detects
single-end vs paired-end from whether `R2` is present.

The pipeline's samplesheet validator will abort early if anything is wrong
(missing columns, missing FASTQ files, fewer than two conditions, fewer than
two replicates per condition).

---

## 5. Step 4 — Download references and build the HISAT2 index

These steps run **once** and are reused for every subsequent pipeline
invocation.

### 5.1 Download the GRCh38 genome, GTF and MSigDB gene sets

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

- `references/human/v110/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz` — the genome FASTA.
- `references/human/v110/Homo_sapiens.GRCh38.110.gtf.gz` — the annotation.
- `references/human/v110/download_log.txt` — provenance (which release, when).
- `references/human/gmt/hallmark.gmt`, `c2_curated.gmt`, `c2_kegg.gmt`, `c2_reactome.gmt`, `c5_go.gmt`, `c5_go_bp.gmt` — MSigDB gene sets for GSEA.

### 5.2 Build the HISAT2 index

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

## 6. Step 5 — Run the analysis

With the samplesheet, references and index in place, run the full RNASEQ
workflow with HISAT2:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
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
  results CSV per contrast — here, `tgfb_vs_REF`. DESeq2 log2 fold changes
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

> **Resume after an interruption.** Re-run the same command with `-resume` —
> Nextflow only re-executes the tasks affected by anything that changed.

---

## 7. Step 6 — What was produced and how to read it

Start with the two reports:

- **`results/quarto_report/analysis_report.html`** — the headline interactive
  report. The DESeq2 and edgeR sections each show: significant-gene counts
  per contrast, an interactive plotly volcano (hover for gene name and
  stats), the embedded PCA / MDS / heatmap, and a searchable results table.
  A DESeq2-vs-edgeR agreement section lists how many genes were significant
  in both methods.
- **`results/multiqc/multiqc_report.html`** — per-sample QC across FastQC,
  fastp, HISAT2 alignment, RSeQC and featureCounts. Use this to confirm
  alignment rates are healthy and strandedness was inferred consistently.

Then the result tables, all under `results/`:

| Path | What's in it |
|---|---|
| `deseq2_output/deseq2_results_tgfb_vs_REF.csv` | DESeq2 table — apeglm-shrunken log2FoldChange + Wald stat / pvalue / padj, with `gene_name` and `gene_biotype` columns |
| `edger_output/edger_results_tgfb_vs_REF.csv` | edgeR counterpart — logFC, logCPM, F-stat, PValue, FDR |
| `gsea_output/gsea_stats_tgfb_vs_REF.csv` + `gsea_plot_tgfb_vs_REF.png` | Hallmark pathways enriched in TGF-β (positive NES = enriched in tgfb) |
| `gprofiler_output/gprofiler_UP_tgfb_vs_REF.csv`, `gprofiler_DOWN_tgfb_vs_REF.csv` | GO terms / pathways over-represented in the up- and down-regulated gene lists, with their `gostplot_*.png` images |
| `featurecounts/`, `rseqc/`, `bam_to_bigwig/`, `hisat2/` | Upstream counts, QC and alignment outputs |
| `rmats_output/` | Alternative-splicing events (SE / MXE / A3SS / A5SS / RI) |
| `pipeline_info/run_summary.html` | Run status, duration, links to all outputs, per-process resource table |

---

## 8. Expected result — a built-in positive control

TGF-β activation of hepatic stellate cells drives an unmistakable, textbook
fibrosis signature. In `deseq2_results_tgfb_vs_REF.csv` you should see strong
upregulation (large positive `log2FoldChange`, very small `padj`) of:

- **COL1A1, COL3A1** — type I/III collagens.
- **ACTA2** (α-smooth-muscle actin) — myofibroblast marker.
- **SERPINE1** (PAI-1) — TGF-β's classical immediate target.
- **IGFBP3** — the dataset's own published top finding.

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

## 9. Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| `Samplesheet validation failed` | Re-check column names; `R1` paths must exist; conditions need ≥ 2 levels and ≥ 2 replicates each. |
| Very few reads `Assigned` by featureCounts | Strandedness mismatch. Check `results/rseqc/<sample>.strandedness.txt` and either leave `--strandedness auto` or set the matching value explicitly. |
| HISAT2 index build runs out of memory | Raise `--max_memory` (default `128.GB`; you need ~8 GB free for the build). |
| No `gsea_output/` | `--gmt` not supplied, or `--stop_at` was set before the enrichment stage. |
| Container pull failures | On HPC prefer `-profile singularity`; locally check Docker is running. |
| The DE table has very few significant genes | You subsampled too aggressively — try increasing `N` (e.g. to 5 M pairs) and re-running. |

---

## 10. Reusing the example

The subsampled FASTQs in `data/subset/`, the references in
`references/human/v110/` and the HISAT2 index in
`indices/human/hisat2_index/` are all one-time setup. Once they exist,
re-running the pipeline against the same samplesheet takes only minutes —
which makes it easy to iterate on options, for example adding
`--diffsplice` to enable differential splicing, swapping to STAR with
`--aligner star --star_index ...`, or trying a different `--gmt` file
(`c2_kegg.gmt`, `c5_go_bp.gmt`).

The whole example is deliberately small, reproducible (subsample seed,
pinned Ensembl release) and self-verifying through the TGF-β signature —
so it serves as both a tutorial and an ongoing setup check.
