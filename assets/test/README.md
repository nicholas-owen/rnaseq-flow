# Yeast test dataset

A small, real *S. cerevisiae* RNA-seq dataset for exercising the whole pipeline
end to end on a laptop or a memory-limited WSL instance.

## Why yeast

The genome is 12 Mb, so a STAR index builds in about 2 GB of RAM and a minute or
two. Human and mouse need far more — mostly for index building, where a
splice-aware HISAT2 index alone needs ~204 GB for human. Every downstream stage
supports yeast: `--organism scerevisiae` is valid for gProfiler, and MSigDB
gene sets are available for *Saccharomyces cerevisiae* via `--download_gmt`
(see the caveat under [GSEA](#gsea-optional) — they are human sets mapped by
orthology, not yeast-native ones).

## The data

**BioProject [PRJDB13901](https://www.ebi.ac.uk/ena/browser/view/PRJDB13901)** —
BY4741 carrying the empty `pTOWug2836` vector, grown in YPD or in NaCl.

| sample | run | condition |
| --- | --- | --- |
| `ypd_rep1..3` | DRR392077, DRR392078, DRR392079 | `REF` (control) |
| `nacl_rep1..3` | DRR392092, DRR392093, DRR392094 | `NaCl` |

Two conditions × three replicates — comfortably above the pipeline's minimum of
two conditions and two replicates each, so DESeq2/edgeR, GSEA and gProfiler all
run. Naming the control `REF` orients every fold change as NaCl-vs-control.

Full runs are 3.5–9.8 M read pairs. The fetch script subsamples to 1 M pairs by
default, which is enough to exercise every process while keeping the download
and the run small.

## Fetch it

```bash
assets/test/fetch_yeast_test_data.sh                    # test_data/, 1M reads
assets/test/fetch_yeast_test_data.sh test_data 250000   # smaller/faster
```

Uses `seqtk` for a proper random subsample if it is installed; otherwise it takes
the first N reads. That fallback is fine for a smoke test — every process still
sees real data — but it is **not** a random sample, so don't draw biological
conclusions from it. `seqtk` is in Ubuntu's universe repo if you want the real
thing:

```bash
sudo apt install seqtk      # or: conda install -c bioconda seqtk
```

then delete `test_data/subset/` and re-run; the downloads in `raw/` are kept.

The script is idempotent: a file is only given its final name once it is
complete (downloads and subsamples are written to `.part` first), so re-running
skips whatever is already there and an interrupted run cannot be mistaken for a
finished one. FASTQ URLs are resolved from the ENA API rather than hard-coded.

It writes `test_data/samplesheet_yeast.csv` with absolute paths.

## Run it

**Pass `-profile test_yeast` on every command.** Without it the defaults in
`conf/base.config` apply — those are sized for mammalian genomes (48 GB for
`process_high`, 38 GB for STAR indexing) and the local executor refuses to
schedule a task larger than the machine:

```
Process requirement exceeds available memory -- req: 38 GB; avail: 15.5 GB
```

The profile caps CPUs/memory/time, and `process.resourceLimits` clamps every
request down to fit. Raise the caps to whatever the machine actually has. A
quick way to tell whether it loaded: the STAR command should show
`--runThreadN 4`, not `--runThreadN 12`.

> **Why a profile rather than `-c`?** `conf/base.config` builds `resourceLimits`
> from `params.max_*` *at parse time*. Profiles are merged before it is
> included, so they take effect; a file passed with `-c` is merged *afterwards*,
> so it can change `params.max_memory` without changing the already-computed
> limits — the caps silently do nothing. `conf/test_yeast.config` now restates
> `resourceLimits` itself so `-c` works too, but the profile is the reliable
> route. Explicit CLI flags (`--max_memory 12.GB`) also work, because Nextflow
> injects those before config parsing.

The flags are repeated in full on each command below rather than held in a shell
variable, so that copying any single command still works.

```bash
# 1. references  (yeast has no dna.primary_assembly; the downloader
#    correctly falls back to dna.toplevel)
nextflow run main.nf --download_refs \
    --download_species saccharomyces_cerevisiae --download_release 116 \
    --outdir refs/yeast \
    -profile test_yeast,docker

# 2. index
nextflow run main.nf --build_indices --aligner star \
    --genome_fasta refs/yeast/v116/*.dna.toplevel.fa.gz \
    --gtf          refs/yeast/v116/*.gtf.gz \
    --outdir idx/yeast \
    -profile test_yeast,docker

# 3. analysis
nextflow run main.nf \
    --input      test_data/samplesheet_yeast.csv \
    --aligner    star \
    --star_index idx/yeast/star_index \
    --gtf        refs/yeast/v116/*.gtf.gz \
    --outdir     results \
    -profile test_yeast,docker
```

If you would rather not pass the config, the same caps as explicit flags:

```bash
--max_cpus 4 --max_memory 12.GB --max_time 4.h
```

### GSEA (optional)

Add `--download_gmt` to step 1 and `--gmt refs/yeast/gmt/c5_go_bp.gmt` to
step 3. Step 1 already runs under `-profile test_yeast`, which sets
`organism = 'scerevisiae'`, so the gene sets come back for yeast without any
extra flag. `--download_gmt` on its own does nothing — it is handled inside the
download workflow, which only runs when `--download_refs` is given.

`--download_gmt` needs outbound HTTPS: from msigdbr 24 the gene sets are
fetched at run time rather than bundled.

**Use `c5_go_bp.gmt` (GO biological process), not `hallmark.gmt`.** Hallmark is
50 human-defined sets; projected onto yeast orthologs, few of them retain enough
genes to test. `gsea.R` warns below 15 overlapping genes and skips a contrast
entirely at zero, so hallmark can produce an empty GSEA step that looks like a
failure. GO:BP has far more sets and much better overlap. Either way, check the
run log for the

```
ranking by <gene symbol|gene ID> - N of M genes overlap the gene sets
```

line: that number tells you whether the GMT is usable before you look at any
result.

> **These are not yeast-native gene sets.** msigdbr maps the human MSigDB
> collections to other species by orthology (`db_species = "HS"`, the default —
> see `assets/download_gmt.R`). That is fine for exercising the GSEA step, which
> is what this test dataset is for, but the resulting enrichments are a weak
> basis for biology: many human sets have no meaningful yeast counterpart, and
> ortholog mapping loses genes in both directions. For yeast GSEA that is meant
> to be interpreted, build a GMT from SGD's own GO annotations
> (`org.Sc.sgd.db`, or the SGD GO slim) instead.

## What to check

- `results/multiqc/multiqc_report.html` — every sample should appear under
  FastQC, fastp, STAR, RSeQC and featureCounts.
- `results/quarto_report/analysis_report.html` — PCA should separate YPD from
  NaCl, and the DESeq2/edgeR volcano plots should be populated.
- `results/deseq2_output/deseq2_results_NaCl_vs_REF.csv` — a salt-stress
  response, so the usual osmotic-stress genes (`GRE2`, `HSP12`, `CTT1`, `STL1`)
  are a reasonable sanity check for the up-regulated set.
- The `STAR_GENOME_GENERATE` log line reporting the genome length and the
  `--genomeSAindexNbases` it derived — 10 for this genome, not STAR's mammalian
  default of 14.
