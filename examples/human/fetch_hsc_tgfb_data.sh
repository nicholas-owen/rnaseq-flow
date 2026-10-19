#!/usr/bin/env bash
#
# Fetch the GSE151251 human RNA-seq set and write a samplesheet.
#
#   examples/human/fetch_hsc_tgfb_data.sh [OUTDIR] [READS_PER_SAMPLE]
#
# Defaults: OUTDIR=hsc_data, READS_PER_SAMPLE=0 (keep full depth)
#
# Dataset: BioProject PRJNA635294 / GSE151251 (ENA). Six runs of primary human
# hepatic stellate cells, three untreated and three TGF-beta treated.
#
#   control (-> condition REF)   SRR11854548  SRR11854549  SRR11854550
#   TGF-beta                     SRR11854551  SRR11854552  SRR11854553
#
# All six are paired-end, 22.9 to 30.9 M pairs each, about 10 GB gzipped in
# total.
#
# Unlike the yeast set this defaults to FULL depth, because it is meant for
# alternative splicing and transcript usage. Those read junctions and isoform
# proportions, and both lose signal to subsampling much faster than gene-level
# counts do: the yeast script's 1 M pairs would leave almost nothing testable.
# Pass a read count to subsample anyway, e.g. 10000000.
#
# Note the scale before running. The FASTQs are the small part. A human STAR
# index is roughly 30 GB of disk and wants about 32 GB of RAM to build; a
# Salmon index for the transcript-usage route is a couple of GB and far less
# demanding. Build only the one for the route you intend to run.
#
# FASTQ URLs are resolved from the ENA API at run time rather than hard-coded,
# so the script keeps working if ENA reorganises its paths.
set -euo pipefail

OUTDIR="${1:-hsc_data}"
READS="${2:-0}"
ACCESSION=PRJNA635294

# Controls first, so the samplesheet's REF condition is unambiguous.
CONTROL=(SRR11854548 SRR11854549 SRR11854550)
TREATED=(SRR11854551 SRR11854552 SRR11854553)

# OUTDIR may be relative, so say plainly where the data is going. Running from
# the wrong directory is the easy way to re-download 10 GB you already have.
mkdir -p "$OUTDIR/raw"
ABS=$(cd "$OUTDIR" && pwd)
echo "Writing to $ABS"

# --- subsampling tool ------------------------------------------------------
# seqtk gives a proper random sample; without it, fall back to taking the first
# N reads, which exercises every process but is not a random sample and should
# not be used for real analysis. Only consulted when READS > 0.
SUBSAMPLE=none
if [ "$READS" -gt 0 ]; then
    if command -v seqtk >/dev/null 2>&1; then
        SUBSAMPLE=seqtk
    else
        SUBSAMPLE=head
        echo "seqtk not found: taking the first $READS pairs rather than a random sample."
        echo "Install seqtk (apt install seqtk, or conda install -c bioconda seqtk) for a real subsample."
    fi
fi

# --- run table -------------------------------------------------------------
RUNINFO="$OUTDIR/runinfo.tsv"
if [ ! -s "$RUNINFO" ]; then
    echo "Fetching the run table from ENA"
    curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACCESSION}&result=read_run&fields=run_accession,sample_title,library_layout,fastq_ftp,fastq_md5,read_count&format=tsv" \
        > "$RUNINFO.part"
    mv "$RUNINFO.part" "$RUNINFO"
fi

field() {   # field <run> <column-number>
    awk -F'\t' -v r="$1" -v c="$2" 'NR>1 && $1==r {print $c}' "$RUNINFO"
}

# --- download --------------------------------------------------------------
# A file is given its final name only once it is complete and its checksum
# matches, so an interrupted run resumes cleanly and a truncated file is never
# mistaken for a finished one.
download() {
    local run="$1"
    local urls md5s i url md5 dest
    urls=$(field "$run" 4 | tr ';' ' ')
    md5s=$(field "$run" 5 | tr ';' ' ')
    [ -n "$urls" ] || { echo "No FASTQ URLs for $run in $RUNINFO" >&2; exit 1; }

    i=0
    for url in $urls; do
        i=$((i + 1))
        md5=$(echo "$md5s" | cut -d' ' -f"$i")
        dest="$OUTDIR/raw/$(basename "$url")"
        if [ -s "$dest" ]; then
            echo "  have $(basename "$dest")"
            continue
        fi
        echo "  get  $(basename "$dest")"
        curl -fsSL --retry 4 --retry-delay 5 -o "$dest.part" "https://${url}"
        if command -v md5sum >/dev/null 2>&1; then
            local got
            got=$(md5sum "$dest.part" | cut -d' ' -f1)
            if [ "$got" != "$md5" ]; then
                rm -f "$dest.part"
                echo "  checksum mismatch for $(basename "$dest")" >&2
                exit 1
            fi
        fi
        mv "$dest.part" "$dest"
    done
}

for run in "${CONTROL[@]}" "${TREATED[@]}"; do
    echo "$run  $(field "$run" 2)"
    download "$run"
done

# --- optional subsample ----------------------------------------------------
READDIR="$OUTDIR/raw"
if [ "$SUBSAMPLE" != none ]; then
    mkdir -p "$OUTDIR/subset"
    READDIR="$OUTDIR/subset"
    for run in "${CONTROL[@]}" "${TREATED[@]}"; do
        for mate in 1 2; do
            src="$OUTDIR/raw/${run}_${mate}.fastq.gz"
            dst="$OUTDIR/subset/${run}_${mate}.fastq.gz"
            [ -s "$dst" ] && continue
            echo "  subsample $(basename "$dst")"
            if [ "$SUBSAMPLE" = seqtk ]; then
                # One seed for both mates, or the pairs fall out of step.
                seqtk sample -s42 "$src" "$READS" | gzip > "$dst.part"
            else
                zcat "$src" | head -n $((READS * 4)) | gzip > "$dst.part"
            fi
            mv "$dst.part" "$dst"
        done
    done
fi

# --- samplesheet -----------------------------------------------------------
# Absolute paths, so the samplesheet works from any launch directory.
SHEET="$OUTDIR/samplesheet_hsc.csv"
{
    echo "sample,R1,R2,condition"
    n=0
    for run in "${CONTROL[@]}"; do
        n=$((n + 1))
        echo "ctrl_rep${n},${ABS}/$(basename "$READDIR")/${run}_1.fastq.gz,${ABS}/$(basename "$READDIR")/${run}_2.fastq.gz,REF"
    done
    n=0
    for run in "${TREATED[@]}"; do
        n=$((n + 1))
        echo "tgfb_rep${n},${ABS}/$(basename "$READDIR")/${run}_1.fastq.gz,${ABS}/$(basename "$READDIR")/${run}_2.fastq.gz,TGFb"
    done
} > "$SHEET"

echo
echo "samplesheet: $SHEET"
echo "reads:       $READDIR  ($(du -sh "$READDIR" 2>/dev/null | cut -f1))"
if [ "$SUBSAMPLE" = head ]; then
    echo
    echo "Reads were taken from the head of each file, not sampled at random."
    echo "Fine for checking the pipeline runs; do not draw conclusions from it."
fi
