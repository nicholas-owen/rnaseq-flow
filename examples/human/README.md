# Human splicing, transcript usage and isoform switching

Six samples of primary human hepatic stellate cells, three untreated and three
treated with TGF-beta, run through Salmon to exercise `--dtu`, `--diffsplice`
and `--isoform_switch`. The expected result is published and checkable.

Dataset: GSE151251 / PRJNA635294, paired-end, 22.9 to 30.9 million pairs per
sample at about 101 bp.

## Why this example exists

The yeast example cannot test any of this. Ensembl R64-1-1 has 7127 genes and
7127 transcripts, one each, so no gene has the two or more isoforms that
transcript usage and isoform switching are defined for. Those analyses run
there and correctly produce empty results.

So without a human dataset, three of the pipeline's analyses had no worked
example and no expected answer anywhere in the repository.

## What you will need

About 60 GB of disk: 10 GB of FASTQ, 20 GB of reference and Salmon index, and
the rest for results and Nextflow's work directory. 22 GB of RAM is enough for
everything here. A machine that can spare 10 cores finishes in a few hours.

The commands below cap resources explicitly. Drop `--max_cpus` and
`--max_memory` on a cluster, or raise them if you have more.

## 1. Fetch the reads

```bash
examples/human/fetch_hsc_tgfb_data.sh
```

This downloads all six runs at full depth into `hsc_data/`, verifies each
file's MD5 against ENA, and writes `hsc_data/samplesheet_hsc.csv` with absolute
paths and the three controls as `REF`. It resumes cleanly if interrupted, and a
file is given its final name only once its checksum matches.

If you already ran [the expression walkthrough](../../rnaseq-flow_Example_Analysis.md),
it used this same script on the same dataset, so `hsc_data/raw/` is already
populated and nothing is downloaded twice.

Note that the two share one samplesheet, and the script rewrites it on every
run to point at whichever depth you asked for. So run the command above before
starting here, and re-run it with `hsc_data 2000000` before going back. Neither
re-download nor re-subsample anything that already exists, so the switch costs
seconds.

**Do not subsample this one.** The yeast example takes 1 million pairs per
sample and that is the right call for gene-level counts. Junction and isoform
evidence is a small fraction of those counts, so subsampling removes the signal
before the test sees it. The script accepts a read count as its second argument
if you want to anyway, which is useful for checking the pipeline runs and
useless for drawing conclusions.

## 2. References

```bash
nextflow run main.nf --download_refs --download_species homo_sapiens --download_release 116 --download_gmt --outdir refs/human --max_cpus 10 --max_memory 22.GB -profile docker
```

Pinning `--download_release 116` matters here: the expected results below were
produced against that annotation, and isoform-level results move with the
transcript set more than gene-level ones do.

Check what landed before using the paths below, since the file names carry the
assembly:

```bash
ls refs/human/v116/
```

## 3. Salmon index

```bash
nextflow run main.nf --build_indices --aligner salmon --transcript_fasta refs/human/v116/Homo_sapiens.GRCh38.cdna.all.fa.gz --genome_fasta refs/human/v116/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --outdir idx/human_salmon --max_cpus 10 --max_memory 22.GB -profile docker
```

A decoy-aware index, so it needs the genome as well as the transcriptome. It
lands at `idx/human_salmon/salmon_index`.

## 4. Run

```bash
nextflow run main.nf --input hsc_data/samplesheet_hsc.csv --aligner salmon --salmon_index idx/human_salmon/salmon_index --gtf refs/human/v116/Homo_sapiens.GRCh38.116.gtf.gz --transcript_fasta refs/human/v116/Homo_sapiens.GRCh38.cdna.all.fa.gz --gmt refs/human/gmt/c5_go_bp.gmt --outdir results_hsc_salmon --dtu --diffsplice --isoform_switch --max_cpus 10 --max_memory 22.GB -profile docker
```

No `--strandedness` flag. Salmon is invoked with `--libType A` and detects the
library type per sample, so the parameter never reaches it. The yeast example
passes `--strandedness unstranded` because its genome-aligner routes need it.

`--transcript_fasta` is required twice over: once to build the index, and again
at run time because `--isoform_switch` is skipped without it.

## 5. What you should see

Counts move a little with Salmon's version and the annotation release. The
named genes are the durable part.

| Analysis | Result |
|---|---|
| DTU | 2131 transcripts at `padj < 0.05` |
| Isoform switching | 490 genes, 734 switches across 596 isoforms; DEXSeq puts genes with DTU at 313 to 522 |
| diffSplice | Per-gene results for every multi-transcript gene tested |

The ten highest-ranked switching genes by q-value:

```
AHCYL1  MYLK  COL3A1  NR1H2  TXN  ZBTB38  DIDO1  PLEKHA8  SRI  C3
```

Two of those are worth recognising. COL3A1 is one of the five fibrosis markers
in this dataset's own published expression signature, now moving at the level
of isoform usage as well as abundance. MYLK is myosin light chain kinase,
central to the contractile apparatus that defines an activated stellate cell.
Neither was looked for.

## 6. The published expectation

Alli-Oke et al. 2025 profiled TGF-beta-induced myofibroblast differentiation in
primary skin fibroblasts, validated a five-event splicing signature by ddPCR
and AS-PCR, and recovered it in independent public RNA-seq of lung and skin
fibroblasts. The five are ACTN1 19A/19B, COL5A1 64A/64B, COL6A3 exon 4, FLNA
exon 30 and TPM1 6a/6b.

Three of the five carry significant transcript usage changes here, in hepatic
stellate cells rather than fibroblasts:

| Gene | Isoforms at `padj < 0.05` | Best `padj` |
|---|---|---|
| FLNA | 2 | 3.9e-05 |
| COL6A3 | 3 | 3.1e-04 |
| TPM1 | 3 | 1.1e-03 |

Check them with:

```bash
grep -iE '"(ACTN1|COL5A1|COL6A3|FLNA|TPM1)"' results_hsc_salmon/dtu_output/dtu_transcript_results.csv | awk -F, '$NF!="NA" && $NF+0<0.05' | cut -d, -f2,4,9
```

ACTN1 and COL5A1 do not reach significance. They were tested rather than
filtered out, so this is a real negative, but a weak one for ACTN1: its best
isoform reaches `padj` 0.139 across 35 tested isoforms, and usage divided 35
ways leaves DEXSeq very little power per isoform. A two-exon swap need not shift
transcript proportions enough to register.

**What this is not.** The published signature is event-level, naming particular
exons. DTU is transcript-level. A gene can change isoform usage without the
named exon moving, and the reverse. These results are consistent with the
signature rather than a confirmation of it. Confirming the events themselves
needs rMATS, which means the STAR route.

## 7. The rMATS half

Not covered here, because it has not been run. `STAR_GENOME_GENERATE` requests
38 GB for a human index, which most workstations cannot give it. Two routes
work:

On a cluster:

```bash
nextflow run main.nf --build_indices --aligner star --genome_fasta refs/human/v116/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --gtf refs/human/v116/Homo_sapiens.GRCh38.116.gtf.gz --outdir idx/human_star -profile sge,singularity
```

Or locally with a sparse index, which roughly halves the memory at the cost of
slower alignment:

```bash
cat > star_sparse.config <<'CFG'
process {
    withName: STAR_GENOME_GENERATE {
        ext.args = '--genomeSAsparseD 2 --limitGenomeGenerateRAM 20000000000'
    }
}
CFG
```

```bash
nextflow run main.nf --build_indices --aligner star --genome_fasta refs/human/v116/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --gtf refs/human/v116/Homo_sapiens.GRCh38.116.gtf.gz --outdir idx/human_star -c star_sparse.config --max_cpus 10 --max_memory 22.GB -profile docker
```

Then the run, where `--read_length 101` matters because rMATS uses it for event
detection and the default is 100:

```bash
nextflow run main.nf --input hsc_data/samplesheet_hsc.csv --aligner star --star_index idx/human_star/star_index --gtf refs/human/v116/Homo_sapiens.GRCh38.116.gtf.gz --gmt refs/human/gmt/c5_go_bp.gmt --outdir results_hsc_star --diffsplice --read_length 101 --strandedness auto --max_cpus 10 --max_memory 22.GB -profile docker
```

rMATS needs no flag of its own: it runs whenever a genome aligner produces BAMs
and `--input` is a samplesheet. Keep `--strandedness auto` so RSeQC records the
inference per sample.

Nobody has checked whether the five published events move in this dataset. If
you run it, that is a genuinely open question and the answer is worth having.

## Reading the results

`OUTPUTS.md` sections 13, 14 and 15 describe `isoform_switch_output/`,
`dtu_output/` and `diffsplice_output/` file by file, including the direction of
each effect size and why the three tables do not join on `gene_id`. Section 18
covers the `reproduce/` folders, which redraw any published figure from the
results alone.
