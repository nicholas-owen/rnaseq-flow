#!/usr/bin/env bash
#
# Fetch a small S. cerevisiae RNA-seq test set and write a samplesheet.
#
#   assets/test/fetch_yeast_test_data.sh [OUTDIR] [READS_PER_SAMPLE]
#
# Defaults: OUTDIR=test_data, READS_PER_SAMPLE=1000000
#
# Dataset: BioProject PRJDB13901 (DDBJ/ENA). Six runs of BY4741 carrying the
# empty pTOWug2836 vector, three in YPD and three in NaCl -- a 2 x 3 design,
# which is the minimum the pipeline's samplesheet validation accepts (>= 2
# conditions, >= 2 replicates each).
#
#   YPD  (control -> condition REF)   DRR392077  DRR392078  DRR392079
#   NaCl (stress)                     DRR392092  DRR392093  DRR392094
#
# Yeast is used deliberately: a 12 Mb genome indexes in ~2 GB of RAM, so the
# whole pipeline runs on a laptop or a memory-limited WSL instance. Human and
# mouse need far more, mostly for STAR/HISAT2 index building.
#
# FASTQ URLs are resolved from the ENA API at run time rather than hard-coded,
# so the script keeps working if ENA reorganises its paths.
set -euo pipefail

OUTDIR=${1:-test_data}
NREADS=${2:-1000000}

RUNS_YPD=(DRR392077 DRR392078 DRR392079)
RUNS_NACL=(DRR392092 DRR392093 DRR392094)

RAW="$OUTDIR/raw"
SUB="$OUTDIR/subset"
mkdir -p "$RAW" "$SUB"

# OUTDIR may be relative, so say plainly where the data is going. Running from
# the wrong directory is the easy way to re-download 3 GB you already have: from
# inside test_data/ the default resolves to test_data/test_data.
echo "Output directory: $(cd "$OUTDIR" && pwd)"
echo "  already downloaded: $(ls -1 "$RAW"/*.fastq.gz 2>/dev/null | wc -l)/12 raw files"
echo "  already subsampled: $(ls -1 "$SUB"/*.fastq.gz 2>/dev/null | wc -l)/12 subsets"

# --- subsampling tool ------------------------------------------------------
# seqtk gives a proper random sample; without it, fall back to taking the first
# N reads, which is fine for a smoke test (it exercises every process) but is
# not a random sample and should not be used for real analysis.
if command -v seqtk >/dev/null 2>&1; then
    SUBSAMPLE=seqtk
    echo "Subsampling with seqtk (random, seed 42)."
else
    SUBSAMPLE=head
    echo "NOTE: seqtk not found - taking the FIRST ${NREADS} reads instead of a"
    echo "      random sample. Fine for a pipeline smoke test; install seqtk"
    echo "      (conda install -c bioconda seqtk) for a real subsample."
fi

fetch_run() {
    local run=$1
    local urls
    # ENA always prepends run_accession to the requested fields, so take the
    # LAST column rather than the first -- $1 is the accession, not the URLs.
    urls=$(curl -sS --fail --max-time 120 -G \
        "https://www.ebi.ac.uk/ena/portal/api/filereport" \
        --data-urlencode "accession=${run}" \
        --data-urlencode "result=read_run" \
        --data-urlencode "fields=fastq_ftp" \
        --data-urlencode "format=tsv" | awk 'NR==2 {print $NF}')

    if [ -z "$urls" ]; then
        echo "ERROR: no FASTQ URLs returned by ENA for ${run}" >&2
        return 1
    fi

    local i=0
    for u in ${urls//;/ }; do
        i=$((i + 1))
        local dest="$RAW/${run}_${i}.fastq.gz"
        if [ -s "$dest" ]; then
            echo "  have ${run}_${i}"
        else
            echo "  downloading ${run}_${i} ..."
            curl -sS --fail --max-time 3600 -o "$dest.part" "https://${u#ftp://}"
            mv "$dest.part" "$dest"
        fi
    done
}

subsample_run() {
    local run=$1
    for i in 1 2; do
        local src="$RAW/${run}_${i}.fastq.gz"
        local dst="$SUB/${run}_${i}.fastq.gz"
        [ -s "$src" ] || continue
        if [ -s "$dst" ]; then
            echo "  have subset ${run}_${i}"
            continue
        fi
        rm -f "$dst.part"
        # Same seed for both mates so the pairs stay in step.
        if [ "$SUBSAMPLE" = seqtk ]; then
            seqtk sample -s42 "$src" "$NREADS" | gzip > "$dst.part"
        else
            # `head` closes the pipe once it has its lines, which sends SIGPIPE
            # to zcat. Under `set -o pipefail` that makes the whole pipeline
            # non-zero and `set -e` aborts the script -- leaving a stray .part
            # file and no samplesheet. Lift pipefail for this one command.
            set +o pipefail
            zcat "$src" | head -n $((NREADS * 4)) | gzip > "$dst.part"
            set -o pipefail
        fi
        if [ ! -s "$dst.part" ]; then
            echo "ERROR: subsampling produced nothing for ${run}_${i}" >&2
            rm -f "$dst.part"
            return 1
        fi
        mv "$dst.part" "$dst"
    done
}

echo
echo "== downloading =="
for r in "${RUNS_YPD[@]}" "${RUNS_NACL[@]}"; do
    echo "${r}:"
    fetch_run "$r"
done

echo
echo "== subsampling to ${NREADS} reads/sample =="
for r in "${RUNS_YPD[@]}" "${RUNS_NACL[@]}"; do
    echo "${r}:"
    subsample_run "$r"
done

# --- samplesheet -----------------------------------------------------------
# Absolute paths, so the sheet works from any working directory.
ABS=$(cd "$SUB" && pwd)
SHEET="$OUTDIR/samplesheet_yeast.csv"
{
    echo "sample,R1,R2,condition"
    n=0
    for r in "${RUNS_YPD[@]}";  do n=$((n+1)); echo "ypd_rep${n},${ABS}/${r}_1.fastq.gz,${ABS}/${r}_2.fastq.gz,REF"; done
    n=0
    for r in "${RUNS_NACL[@]}"; do n=$((n+1)); echo "nacl_rep${n},${ABS}/${r}_1.fastq.gz,${ABS}/${r}_2.fastq.gz,NaCl"; done
} > "$SHEET"

echo
echo "== done =="
du -sh "$SUB" 2>/dev/null | sed 's/^/  subset size: /'
echo "  samplesheet: $SHEET"
echo
cat "$SHEET" | sed 's/^/    /'
