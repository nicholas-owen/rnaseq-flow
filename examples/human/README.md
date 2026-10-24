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
nextflow run main.nf --build_indices --aligner salmon --genome_fasta refs/human/v116/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --gtf refs/human/v116/Homo_sapiens.GRCh38.116.gtf.gz --outdir idx/human_salmon --max_cpus 10 --max_memory 22.GB -profile docker
```

Two things come out of this: `idx/human_salmon/salmon_index`, and
`idx/human_salmon/transcripts.fa.gz`, the transcriptome `gffread` cut from the
genome using the GTF.

The transcriptome is derived rather than downloaded because Ensembl's own files
disagree with the GTF both ways. `cdna.all.fa.gz` carries 14,212 patch and
haplotype transcripts the primary-assembly GTF omits, whose copies take reads
from their primary counterparts. It also holds only 465,769 transcripts against
the GTF's 646,577, because it is cDNA only: the roughly 198,000 missing are
almost all lncRNA, shipped separately in Ensembl's `ncrna.fa`. Most lncRNAs are
poly-adenylated and so are captured by this library. Deriving the transcriptome
makes the index and the annotation describe the same transcripts by
construction, and includes the lncRNAs. The genome is still needed here because
it also supplies Salmon's decoys.

> **This changes results against an index built from `cdna.all.fa.gz`**, which
> quantified roughly 40% fewer features. The numbers in section 5 come from the
> derived transcriptome. Against the old one the counts were lower and ACTN1,
> one of the five signature genes, did not reach significance at all. Do not
> compare transcript-level results either side of this change.

## 4. Run

```bash
nextflow run main.nf --input hsc_data/samplesheet_hsc.csv --aligner salmon --salmon_index idx/human_salmon/salmon_index --gtf refs/human/v116/Homo_sapiens.GRCh38.116.gtf.gz --transcript_fasta idx/human_salmon/transcripts.fa.gz --gmt refs/human/gmt/c5_go_bp.gmt --outdir results_hsc_salmon --dtu --diffsplice --isoform_switch --max_cpus 10 --max_memory 22.GB -profile docker
```

No `--strandedness` flag. Salmon is invoked with `--libType A` and detects the
library type per sample, so the parameter never reaches it. The yeast example
passes `--strandedness unstranded` because its genome-aligner routes need it.

`--transcript_fasta` points at the derived `transcripts.fa.gz` from step 3, not
at the downloaded cDNA. `--isoform_switch` is skipped without it, and it has to
be the file the index was built from or the mismatch comes straight back.

## 5. What you should see

Counts move a little with Salmon's version and the annotation release. The
named genes are the durable part.

| Analysis | Result |
|---|---|
| DTU | 2287 transcripts at `padj < 0.05`, from 112,752 transcripts across 12,892 multi-transcript genes |
| Isoform switching | 538 genes, 773 switches across 643 isoforms; DEXSeq puts genes with DTU at 494 to 824 |
| diffSplice | 1027 genes at `FDR < 0.05` |

The ten highest-ranked switching genes by q-value:

```
AHCYL1  MYLK  TXN  NR1H2  FLNB  ZBTB38  SRI  PARP1  C3  GPSM3
```

MYLK at rank 2 is the one worth recognising: myosin light chain kinase, central
to the contractile apparatus that defines an activated stellate cell, and not
something this analysis was looking for.

## 6. The published expectation

Alli-Oke et al. 2025 profiled TGF-beta-induced myofibroblast differentiation in
primary skin fibroblasts, validated a five-event splicing signature by ddPCR
and AS-PCR, and recovered it in independent public RNA-seq of lung and skin
fibroblasts. The five are ACTN1 19A/19B, COL5A1 64A/64B, COL6A3 exon 4, FLNA
exon 30 and TPM1 6a/6b.

Four of the five carry significant transcript usage changes here, in hepatic
stellate cells rather than fibroblasts:

| Gene | Isoforms at `padj < 0.05` | Best `padj` |
|---|---|---|
| ACTN1 | 1 | 1.5e-16 |
| TPM1 | 4 | 2.6e-06 |
| FLNA | 2 | 3.7e-05 |
| COL6A3 | 3 | 3.3e-04 |

Check them with:

```bash
grep -iE '"(ACTN1|COL5A1|COL6A3|FLNA|TPM1)"' results_hsc_salmon/dtu_output/dtu_transcript_results.csv | awk -F, '$NF!="NA" && $NF+0<0.05' | cut -d, -f2,4,9
```

ACTN1 is the strongest of the four, and it is also the clearest demonstration
of why the transcriptome is derived from the GTF rather than downloaded. The
same transcript, `ENST00001044660`, was tested against an index built from
Ensembl's `cdna.all.fa.gz` and reached only `padj` 0.25. Isoform usage is a
proportion over a gene's transcripts, so when some of ACTN1's were missing from
the index, reads from the absent ones were forced onto the present ones and
blurred the proportions. With the full transcript set the 19A/19B usage shift
is unambiguous.

COL5A1 does not reach significance. It was tested rather than filtered out, so
that is a real negative rather than a gap in the data.

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
