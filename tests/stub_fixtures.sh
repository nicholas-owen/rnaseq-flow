#!/usr/bin/env bash
#
# Build the smallest input set a -stub-run needs.
#
# A stub run executes each process's stub: block instead of its script, so no
# tool, no container and no real data are involved. Nextflow still stages the
# declared input files though, so they have to exist on disk. Everything here is
# empty: nothing reads a byte of it.
#
# Usage: bash tests/stub_fixtures.sh <dir>
#
# Writes <dir>/samplesheet.csv, <dir>/genome.gtf.gz, <dir>/transcripts.fa.gz,
# <dir>/star_index/, <dir>/hisat2_index/, <dir>/salmon_index/ and eight empty
# FASTQ files under <dir>/reads/.
#
# Conditions are REF and TREAT. Both are valid R names, which the samplesheet
# validator now requires: DESeq2 builds its coefficient names from
# make.names()-sanitised levels, so a condition containing a hyphen or leading
# digit fails at results(name=) after everything expensive has run.
set -euo pipefail

DEST="${1:?usage: stub_fixtures.sh <dir>}"
mkdir -p "$DEST"
# Absolute, because the samplesheet records read paths and Nextflow resolves
# them against its own launch directory, not against this script's.
DEST="$(cd "$DEST" && pwd)"

mkdir -p "$DEST/reads" "$DEST/star_index" "$DEST/salmon_index" "$DEST/hisat2_index"

# Index directories are passed as paths and never opened by a stub, but an
# empty directory does not survive being staged, so each gets one file.
: > "$DEST/star_index/SA"
: > "$DEST/salmon_index/info.json"
# HISAT2, because it is the only genome aligner that takes no GTF at alignment
# time. That makes it the route where a run with no --gtf is legitimate, and so
# the one that can test differential expression skipping cleanly. STAR hands the
# annotation to --sjdbGTFfile and now refuses up front without it.
: > "$DEST/hisat2_index/genome.1.ht2"
: > "$DEST/genome.gtf.gz"
# For --isoform_switch, which needs the transcript sequences. Present so a test
# can reach the analysis-level gate rather than stopping at the missing-file
# notice before it.
: > "$DEST/transcripts.fa.gz"

printf 'sample,R1,R2,condition\n' > "$DEST/samplesheet.csv"
for sample in ctrl_rep1 ctrl_rep2 treat_rep1 treat_rep2; do
    : > "$DEST/reads/${sample}_1.fastq.gz"
    : > "$DEST/reads/${sample}_2.fastq.gz"
    case "$sample" in
        treat_*) condition=TREAT ;;
        *)       condition=REF ;;
    esac
    printf '%s,%s/reads/%s_1.fastq.gz,%s/reads/%s_2.fastq.gz,%s\n' \
        "$sample" "$DEST" "$sample" "$DEST" "$sample" "$condition" \
        >> "$DEST/samplesheet.csv"
done

echo "stub fixtures in $DEST"
cat "$DEST/samplesheet.csv"
