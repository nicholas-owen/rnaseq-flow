#!/usr/bin/env python3

import sys
import os
import csv
import itertools
import subprocess

# Usage: run_rmats.py <samplesheet.csv> <gtf> <read_type> <read_length> <output_dir> <sample_bam.tsv>
#
# <sample_bam.tsv> is a two-column, tab-separated sample->BAM mapping written by
# the RMATS process from Nextflow's own meta.id. BAMs used to be matched to
# samples here by filename substring, and sample 'ctrl1' matched 'ctrl10.bam':
# wrong BAM, wrong splicing calls, no error (C3). The mapping is authoritative,
# so every disagreement below is an internal consistency failure -- none is
# reachable from an ordinary samplesheet -- and each exits hard rather than
# warning, because a warning here means rMATS runs on the wrong grouping.

def main():
    if len(sys.argv) != 7:
        print("Usage: run_rmats.py <samplesheet.csv> <gtf> <read_type> <read_length> <output_dir> <sample_bam.tsv>")
        sys.exit(1)

    samplesheet = sys.argv[1]
    gtf = sys.argv[2]
    read_type = sys.argv[3] # "paired" or "single"
    read_length = sys.argv[4]
    base_output_dir = sys.argv[5]
    bam_map_path = sys.argv[6]

    # Parse Samplesheet
    samples = {} # sample_id -> condition
    conditions = set()
    with open(samplesheet, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            samples[row['sample']] = row['condition']
            conditions.add(row['condition'])

    # Read the sample->BAM mapping, collecting every problem before failing so
    # one run reports them all (the same style as the samplesheet validation).
    sample_bams = {}
    problems = []
    with open(bam_map_path, 'r') as f:
        for lineno, line in enumerate(f, 1):
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) != 2:
                problems.append(f"line {lineno} of {bam_map_path}: expected 'sample<TAB>bam', got: {line!r}")
                continue
            sample, bam = parts
            if sample not in samples:
                problems.append(f"BAM '{bam}' is mapped to '{sample}', which is not in the samplesheet")
            elif sample in sample_bams:
                problems.append(f"sample '{sample}' is mapped to two BAMs: '{sample_bams[sample]}' and '{bam}'")
            elif not os.path.exists(bam):
                problems.append(f"BAM for sample '{sample}' not found on disk: '{bam}'")
            else:
                sample_bams[sample] = bam

    missing = sorted(set(samples) - set(sample_bams))
    if missing and not problems:
        problems.append("no BAM mapped for sample(s): " + ", ".join(missing))

    if problems:
        print("ERROR: the sample->BAM mapping disagrees with the samplesheet. "
              "This is a pipeline bug, not a samplesheet mistake:")
        for p in problems:
            print(f"  - {p}")
        sys.exit(1)

    # Sort conditions
    sorted_conds = sorted(list(conditions))
    # Move REF to start
    if 'REF' in sorted_conds:
        sorted_conds.remove('REF')
        sorted_conds.insert(0, 'REF')
    
    pairs = list(itertools.combinations(sorted_conds, 2))
    
    if not os.path.exists(base_output_dir):
        os.makedirs(base_output_dir)

    for c1, c2 in pairs:
        print(f"Running rMATS for {c1} vs {c2}...")
        
        # Every sample is guaranteed a BAM by the validation above, so an empty
        # group here can only mean that guarantee was broken. Fail rather than
        # skip: the old silent 'continue' is how a wrong grouping got analysed.
        b1 = [sample_bams[s] for s in samples if samples[s] == c1]
        b2 = [sample_bams[s] for s in samples if samples[s] == c2]

        if not b1 or not b2:
            print(f"ERROR: contrast {c1} vs {c2} has an empty group after validation; "
                  "this should be unreachable.")
            sys.exit(1)
            
        # Write temporary b1.txt and b2.txt
        out_subdir = os.path.join(base_output_dir, f"{c1}_vs_{c2}")
        if not os.path.exists(out_subdir):
            os.makedirs(out_subdir)
            
        b1_path = os.path.join(out_subdir, "b1.txt")
        b2_path = os.path.join(out_subdir, "b2.txt")
        
        with open(b1_path, 'w') as f: f.write(",".join(b1))
        with open(b2_path, 'w') as f: f.write(",".join(b2))
        
        # Construct rMATS command.
        #
        # --variable-read-length is essential here: --readLength is a single
        # number, and by default rMATS DISCARDS every read that is not exactly
        # that length. Reads reaching this step have been adapter/quality
        # trimmed by fastp, so their lengths vary -- without this flag most of
        # the data is silently thrown away and the splicing results are
        # meaningless rather than obviously wrong.
        cmd = [
            "rmats.py",
            "--b1", b1_path,
            "--b2", b2_path,
            "--gtf", gtf,
            "-t", read_type,
            "--readLength", str(read_length),
            "--variable-read-length",
            "--nthread", "4",
            "--od", out_subdir,
            "--tmp", os.path.join(out_subdir, "tmp")
        ]
        
        # Run
        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
