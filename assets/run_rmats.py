#!/usr/bin/env python3

import sys
import os
import csv
import itertools
import subprocess

# Usage: run_rmats.py <samplesheet.csv> <gtf> <read_type> <read_length> <output_dir>
#                     <sample_ids_csv> <bam_files...>
#
# <sample_ids_csv> is a comma-separated list of sample ids in the SAME ORDER as
# the BAM arguments that follow. It comes from Nextflow's own meta.id, so the
# pairing is authoritative rather than reconstructed. BAMs used to be matched to
# samples here by filename substring, and sample 'ctrl1' matched 'ctrl10.bam':
# wrong BAM, wrong splicing calls, no error (C3).
#
# Because the mapping is authoritative, every disagreement below is an internal
# consistency failure -- none is reachable from an ordinary samplesheet -- so
# each exits hard rather than warning. A warning here would mean rMATS running
# on the wrong grouping.

def main():
    if len(sys.argv) < 8:
        print("Usage: run_rmats.py <samplesheet.csv> <gtf> <read_type> <read_length> "
              "<output_dir> <sample_ids_csv> <bam_files...>")
        sys.exit(1)

    samplesheet = sys.argv[1]
    gtf = sys.argv[2]
    read_type = sys.argv[3] # "paired" or "single"
    read_length = sys.argv[4]
    base_output_dir = sys.argv[5]
    sample_ids = [s for s in sys.argv[6].split(',') if s]
    bam_files = sys.argv[7:]

    # Parse Samplesheet
    samples = {} # sample_id -> condition
    conditions = set()
    with open(samplesheet, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            samples[row['sample']] = row['condition']
            conditions.add(row['condition'])

    # Pair ids with BAMs by position, collecting every problem before failing so
    # one run reports them all (the same style as the samplesheet validation).
    problems = []
    if len(sample_ids) != len(bam_files):
        problems.append(f"{len(sample_ids)} sample id(s) but {len(bam_files)} BAM(s); "
                        "the workflow passed mismatched lists")

    sample_bams = {}
    for sample, bam in zip(sample_ids, bam_files):
        if sample not in samples:
            problems.append(f"BAM '{bam}' is paired with '{sample}', which is not in the samplesheet")
        elif sample in sample_bams:
            problems.append(f"sample '{sample}' appears twice: '{sample_bams[sample]}' and '{bam}'")
        elif not os.path.exists(bam):
            problems.append(f"BAM for sample '{sample}' not found on disk: '{bam}'")
        else:
            sample_bams[sample] = bam

    missing = sorted(set(samples) - set(sample_bams))
    if missing and not problems:
        problems.append("no BAM for sample(s): " + ", ".join(missing))

    if problems:
        print("ERROR: the sample-to-BAM pairing disagrees with the samplesheet. "
              "This is a pipeline bug, not a samplesheet mistake:")
        for p in problems:
            print(f"  - {p}")
        sys.exit(1)

    # Order conditions so the baseline is LAST, which puts it second in every pair and
    # therefore makes it the denominator -- matching how DESeq2 and edgeR orient
    # every contrast in this pipeline (see run_contrast() in assets/deseq2.R).
    #
    # This is load-bearing for interpretation, not cosmetics. itertools.
    # combinations takes the pair in list order, c1 becomes both the directory
    # prefix and --b1, and rMATS defines
    #     IncLevelDifference = mean(IncLevel1) - mean(IncLevel2)
    # where sample 1 is --b1. With REF first the directory read REF_vs_NaCl and
    # a positive value meant "higher inclusion in REF" -- the opposite of a
    # positive log2FoldChange in the DE tables sitting beside it. Anyone
    # comparing the two concluded they agreed when they disagreed (M19).
    #
    # REF last gives NaCl_vs_REF and PSI(NaCl) - PSI(REF): positive means higher
    # inclusion in the treatment, exactly as positive log2FC does.
    # The baseline is --reference_level (default "REF"), delivered as an
    # environment variable; see the note in assets/deseq2.R for why it travels
    # that way rather than as an argument.
    ref_level = os.environ.get('RNASEQ_FLOW_REFERENCE_LEVEL', 'REF')

    sorted_conds = sorted(list(conditions))
    if ref_level in sorted_conds:
        sorted_conds.remove(ref_level)
        sorted_conds.append(ref_level)
    else:
        # deseq2.R warns in this case rather than failing; match that, because a
        # silent fallback to alphabetical order decides which way every fold
        # change and every PSI difference points. main.nf fails the run at launch
        # when a level named explicitly is absent, so reaching here means the
        # default was in force and no condition happened to be called REF.
        print(f"Warning: baseline condition '{ref_level}' not found in the "
              "samplesheet. Using alphabetical order, so the last condition "
              "alphabetically becomes the denominator of each contrast. Set "
              "--reference_level, or name a condition 'REF', to control this.")

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
