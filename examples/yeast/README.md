# Yeast test dataset

A small, real *S. cerevisiae* RNA-seq dataset for exercising the whole pipeline
end to end on a laptop or a memory-limited WSL instance.

## Why yeast

The genome is 12 Mb, so a STAR index builds in about 2 GB of RAM and a minute or
two. Human and mouse need far more, mostly for index building, where a
splice-aware HISAT2 index alone needs ~204 GB for human. Every downstream stage
supports yeast: `--organism scerevisiae` is valid for gProfiler, and MSigDB
gene sets are available for *Saccharomyces cerevisiae* via `--download_gmt`
(see the caveat under [GSEA](#gsea-optional): they are human sets mapped by
orthology, not yeast-native ones).

## The data

**BioProject [PRJDB13901](https://www.ebi.ac.uk/ena/browser/view/PRJDB13901)**:
BY4741 carrying the empty `pTOWug2836` vector, grown in YPD or in NaCl.

| sample | run | condition |
| --- | --- | --- |
| `ypd_rep1..3` | DRR392077, DRR392078, DRR392079 | `REF` (control) |
| `nacl_rep1..3` | DRR392092, DRR392093, DRR392094 | `NaCl` |

Two conditions × three replicates, comfortably above the pipeline's minimum of
two conditions and two replicates each, so DESeq2/edgeR, GSEA and gProfiler all
run. Naming the control `REF` orients every fold change, and every rMATS PSI
difference, as NaCl-vs-control.

That is why none of the commands below pass `--reference_level`: it defaults to
`REF`, which is what this dataset calls its control, so the flag would be a
no-op. Adapting these commands to data whose control group has another name is
the case that needs it: `--reference_level Untreated`, say. See
[USAGE.md §4.5](../USAGE.md#45-which-condition-results-are-measured-against---reference_level).

Full runs are 3.5–9.8 M read pairs. The fetch script subsamples to 1 M pairs by
default, which is enough to exercise every process while keeping the download
and the run small.

## Fetch it

```bash
examples/yeast/fetch_yeast_test_data.sh                    # test_data/, 1M reads
examples/yeast/fetch_yeast_test_data.sh test_data 250000   # smaller/faster
```

Uses `seqtk` for a proper random subsample if it is installed; otherwise it takes
the first N reads. That fallback is fine for a smoke test (every process still
sees real data), but it is **not** a random sample, so don't draw biological
conclusions from it.

The seqtk path is seeded (`-s42`, the same seed for both mates so pairs stay in
step), so the subsample is reproducible: the same accessions at the same read
count give the same reads on any machine. Which of the two paths was taken is
worth noting before you compare against anyone else's numbers: the script
prints it, and the two produce different data from the same input. `seqtk` is in Ubuntu's universe repo if you want the real
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
`conf/base.config` apply: those are sized for mammalian genomes (48 GB for
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
> limits: the caps silently do nothing. `conf/test_yeast.config` now restates
> `resourceLimits` itself so `-c` works too, but the profile is the reliable
> route. Explicit CLI flags (`--max_memory 12.GB`) also work, because Nextflow
> injects those before config parsing.

The flags are repeated in full on each command below rather than held in a shell
variable, so that copying any single command still works.

Every route is three steps: download the references once, build the index for the
aligner you want, then run the analysis. Step 1 is shared by all four.

```bash
# 1. references  (yeast has no dna.primary_assembly; the downloader
#    correctly falls back to dna.toplevel). --download_gmt also fetches
#    the gene sets -- see GSEA below.
nextflow run main.nf --download_refs \
    --download_species saccharomyces_cerevisiae --download_release 116 \
    --download_gmt \
    --outdir refs/yeast \
    -profile test_yeast,docker
```

### All four aligner paths

The pipeline supports two genome aligners, which produce BAMs and count genes
with featureCounts, and two pseudo-aligners, which quantify transcripts and are
summarised to gene level by tximport:

| `--aligner` | index built from | index param | count route |
| --- | --- | --- | --- |
| `star` | genome + GTF | `--star_index` | BAM → featureCounts |
| `hisat2` | genome + GTF | `--hisat2_index` | BAM → featureCounts |
| `salmon` | transcripts + genome (decoys) | `--salmon_index` | tximport |
| `kallisto` | transcripts | `--kallisto_index` | tximport |

Only the BAM routes produce RSeQC, bigWig and rMATS output; the pseudo-aligner
routes skip straight from quantification to differential expression. Everything
downstream of the count matrix (DESeq2, edgeR, GSEA, gProfiler, MultiQC, the
Quarto report and the `reproduce/` folders) is identical on all four.

**STAR**

```bash
nextflow run main.nf --build_indices --aligner star \
    --genome_fasta refs/yeast/v116/*.dna.toplevel.fa.gz \
    --gtf          refs/yeast/v116/*.gtf.gz \
    --outdir idx/yeast_star \
    -profile test_yeast,docker

nextflow run main.nf \
    --input      test_data/samplesheet_yeast.csv \
    --aligner    star \
    --star_index idx/yeast_star/star_index \
    --gtf        refs/yeast/v116/*.gtf.gz \
    --gmt        refs/yeast/gmt/c5_go_bp.gmt \
    --outdir     results_star \
    -profile test_yeast,docker
```

**HISAT2**

```bash
nextflow run main.nf --build_indices --aligner hisat2 \
    --genome_fasta refs/yeast/v116/*.dna.toplevel.fa.gz \
    --gtf          refs/yeast/v116/*.gtf.gz \
    --outdir idx/yeast_hisat2 \
    -profile test_yeast,docker

nextflow run main.nf \
    --input        test_data/samplesheet_yeast.csv \
    --aligner      hisat2 \
    --hisat2_index idx/yeast_hisat2/hisat2_index \
    --strandedness reverse \
    --gtf          refs/yeast/v116/*.gtf.gz \
    --gmt          refs/yeast/gmt/c5_go_bp.gmt \
    --outdir       results_hisat2 \
    -profile test_yeast,docker
```

**Salmon**: the transcriptome is cut from the genome with `gffread` using the
GTF, so the index and the annotation describe the same transcripts. The genome
is needed twice over, for that and for the decoy set. It also writes
`idx/yeast_salmon/transcripts.fa.gz`, which is the file to pass to
`--transcript_fasta` on any later run that needs one.

```bash
nextflow run main.nf --build_indices --aligner salmon \
    --genome_fasta refs/yeast/v116/*.dna.toplevel.fa.gz \
    --gtf          refs/yeast/v116/*.gtf.gz \
    --outdir idx/yeast_salmon \
    -profile test_yeast,docker

nextflow run main.nf \
    --input        test_data/samplesheet_yeast.csv \
    --aligner      salmon \
    --salmon_index idx/yeast_salmon/salmon_index \
    --gtf          refs/yeast/v116/*.gtf.gz \
    --gmt          refs/yeast/gmt/c5_go_bp.gmt \
    --outdir       results_salmon \
    -profile test_yeast,docker
```

**Kallisto**: transcripts only, no decoys, but the genome is still needed to
cut the transcriptome from.

```bash
nextflow run main.nf --build_indices --aligner kallisto \
    --genome_fasta refs/yeast/v116/*.dna.toplevel.fa.gz \
    --gtf          refs/yeast/v116/*.gtf.gz \
    --outdir idx/yeast_kallisto \
    -profile test_yeast,docker

nextflow run main.nf \
    --input          test_data/samplesheet_yeast.csv \
    --aligner        kallisto \
    --kallisto_index idx/yeast_kallisto/kallisto_index \
    --strandedness   reverse \
    --gtf            refs/yeast/v116/*.gtf.gz \
    --gmt            refs/yeast/gmt/c5_go_bp.gmt \
    --outdir         results_kallisto \
    -profile test_yeast,docker
```

> **Why `--strandedness reverse` on HISAT2 and Kallisto, and nothing on Salmon.**
> This library is reverse-stranded. RSeQC puts 0.011 of reads on the sense
> strand against 0.896 antisense, the same to within 0.01 across all six
> samples, which is the signature of a dUTP protocol. STAR is left on `auto`
> above, where RSeQC infers exactly that, so the inference path gets exercised
> at least once.
>
> `auto` needs a BAM, which the pseudo-aligners do not produce, so the value has
> to be stated for them. Kallisto takes `reverse` as `--rf-stranded`. Salmon is
> given no value because it would ignore one: it is invoked with `--libType A`
> and detects the library type itself.
>
> These three carried `--strandedness unstranded` until 2026-10-30, on the
> mistaken premise that the dataset was unstranded. featureCounts therefore ran
> `-s 0` on a stranded library and Kallisto quantified without a strand flag.
> The four-route table below was regenerated after the correction.

The `--gtf` is required on the pseudo-aligner runs too: tximport needs it to
build the transcript-to-gene map.

If you would rather not pass the config, the same caps as explicit flags:

```bash
--max_cpus 4 --max_memory 12.GB --max_time 4.h
```

### GSEA (optional)

`--download_gmt` is already on the reference command above, and `--gmt
refs/yeast/gmt/c5_go_bp.gmt` on each of the four run commands. The download runs
under `-profile test_yeast`, which sets `organism = 'scerevisiae'`, so the gene
sets come back for yeast without any extra flag. `--download_gmt` on its own
does nothing: it is handled inside the download workflow, which only runs when
`--download_refs` is given.

> **If you downloaded references before v1.5.0**, the gene sets were published
> one level too deep, at `refs/yeast/gmt/gmt/c5_go_bp.gmt`. The path above is
> the corrected one. A re-download writes the un-nested copy alongside the old
> directory rather than replacing it, so check which one you are pointing at:
> both will exist and both are readable, which makes the mistake quiet.

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
> collections to other species by orthology (`db_species = "HS"`, the default:
> see `assets/download_gmt.R`). That is fine for exercising the GSEA step, which
> is what this test dataset is for, but the resulting enrichments are a weak
> basis for biology: many human sets have no meaningful yeast counterpart, and
> ortholog mapping loses genes in both directions. For yeast GSEA that is meant
> to be interpreted, build a GMT from SGD's own GO annotations
> (`org.Sc.sgd.db`, or the SGD GO slim) instead.

## What to check

Replace `results_*` below with whichever output directory you used.

- `results_*/multiqc/rnaseq-flow_multiqc_report.html`: every sample should appear, under
  FastQC, fastp and the modules for the route you ran: STAR/HISAT2, RSeQC and
  featureCounts for the genome aligners, Salmon for the Salmon run. Kallisto has
  no MultiQC module, so that run shows only FastQC and fastp. All four should
  also show the `sample_provenance` custom table mapping sample IDs back to run
  accessions.
- `results_*/quarto_report/analysis_report.html`: PCA should separate YPD from
  NaCl, and the DESeq2/edgeR volcano plots should be populated.
- `results_*/deseq2_output/deseq2_results_NaCl_vs_REF.csv`: a salt-stress
  response, so the usual osmotic-stress genes (`GRE2`, `HSP12`, `CTT1`, `STL1`)
  are a reasonable sanity check for the up-regulated set.
- `results_*/deseq2_output/reproduce/`: each of the four routes should produce
  the same set of standalone redraw scripts (5 for DESeq2, 3 for edgeR, 2 for
  GSEA, 1 for gProfiler).
- STAR only: the `STAR_GENOME_GENERATE` log line reporting the genome length and
  the `--genomeSAindexNbases` it derived: 10 for this genome, not STAR's
  mammalian default of 14.
- HISAT2 only: the per-sample `results_hisat2/hisat2/*.summary.log` overall
  alignment rate, ~96–98% on this data.

### Cross-aligner concordance

Running more than one route is the cheapest real test of the pipeline, because
the four are independent implementations that should reach the same biology.

Counted from `deseq2_output/deseq2_results_NaCl_vs_REF.csv`: a row is tested if
it is in the file, and significant if `padj < 0.05` and `|log2FoldChange| > 1`,
taking the shrunk `log2FoldChange` column rather than `log2FoldChange_MLE`. Up
and down split that set by the sign of the same column. Stating the column
matters: the two differ for genes sitting on the threshold, which is enough to
move a count by one.

Count these with a CSV reader that respects quoting, not by splitting on commas.
The yeast annotation contains `ADE5,7`, whose gene name holds a comma, so a
naive split shifts every column after it by one and reads `baseMean` where it
wants `log2FoldChange`. That turns one strongly down-regulated gene into an
up-regulated one, in three rows of 5579, and it is exactly how the first
regeneration of this table got the up/down split wrong by one per route.

| route | genes tested | significant | up | down |
| --- | --- | --- | --- | --- |
| STAR | 5579 | 852 | 152 | 700 |
| HISAT2 | 5587 | 847 | 153 | 694 |
| Salmon | 5508 | 778 | 130 | 648 |
| Kallisto | 5532 | 854 | 152 | 702 |

735 genes were called significant by all four, no pairwise overlap fell below
95% of the smaller set, and all four reproduced the same strong down-regulation
bias.

These are the DESeq2 calls. edgeR on the same data is consistently a little more
generous, having gained power in version 4.x: 924, 930, 895 and 944 significant
for STAR, HISAT2, Salmon and Kallisto respectively, 814 called by all four, and
no pairwise overlap below 91%. Its columns are `FDR` and `logFC` rather than
`padj` and `log2FoldChange`, at the same thresholds. The direction of the bias
is the same, though less lopsided: edgeR calls roughly twice as many genes up.

> **These numbers were regenerated on 2026-10-30 and the previous ones were
> wrong twice over.** HISAT2 and Kallisto ran `--strandedness unstranded` on a
> reverse-stranded library, and the Salmon and Kallisto transcriptome changed on
> 2026-10-24 without the table being redone. Correcting both raised the genes
> called by all four from 668 to 735 and the worst pairwise overlap from 83% to
> 95%, which is the direction a real fix should move a concordance figure.
>
> STAR did not move. It was always on `auto`, which infers `reverse`, and it
> uses no transcriptome, so its row reproduces the 2026-10-18 figures exactly:
> 5579 tested, 852 significant, 152 up, 700 down, and 924 for edgeR. Salmon's
> significant counts are unchanged too, because `--libType A` ignored the
> strandedness flag; only its `tested` moved, by two genes, from the
> transcriptome change. That is the control for everything else in the table.

> **These numbers come from a subsample, so treat them as a shape to match, not
> targets to hit.** They were produced from the fetch script's default output:
> 1 M read pairs per sample, drawn by `seqtk sample -s42` (seqtk 1.4) from raw
> libraries of 3.5–9.8 M pairs. Because the seed is fixed, seqtk users at the
> default depth should land very close to the table, but you will *not* match
> it if seqtk was missing when you fetched (the script falls back to the first
> N reads, a different set of reads entirely, and says so in its output), if you
> passed a different read count, or if container versions have moved. Any of
> those shifts the absolute counts without meaning anything is wrong.

What matters is the agreement between routes, not the absolute counts: all four
moving together is the expected result, whereas one route disagreeing sharply
with the other three is worth investigating.

The overlaps split along the genome-aligner / pseudo-aligner line, as you would
hope: STAR and HISAT2 agree on 99% of the smaller set, while every pair that
crosses the line sits at 95% or 96%. The routes differ in multimapper handling
and length correction, so some gap across that line is expected.

The previous version of this section reported the opposite, with STAR closer to
Salmon than to HISAT2, and excused it. That was the strandedness bug: HISAT2
was being counted unstranded on a reverse-stranded library, which pushed it away
from the aligner it most resembles. An anomaly that needs excusing is worth
checking first.
