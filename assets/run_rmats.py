#!/usr/bin/env python3

import sys
import os
import csv
import itertools
import subprocess

# Usage: run_rmats.py <samplesheet.csv> <gtf> <read_type> <read_length> <output_dir> <bam_files...>

def main():
    if len(sys.argv) < 7:
        print("Usage: run_rmats.py <samplesheet.csv> <gtf> <read_type> <read_length> <output_dir> <bam_files...>")
        sys.exit(1)

    samplesheet = sys.argv[1]
    gtf = sys.argv[2]
    read_type = sys.argv[3] # "paired" or "single"
    read_length = sys.argv[4]
    base_output_dir = sys.argv[5]
    bam_files = sys.argv[6:]

    # Parse Samplesheet
    samples = {} # sample_id -> condition
    conditions = set()
    with open(samplesheet, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            samples[row['sample']] = row['condition']
            conditions.add(row['condition'])

    # Map BAMs to Samples
    # We assume filename contains sample id. 
    # Logic: iterate bams, find matching sample id.
    sample_bams = {}
    for bam in bam_files:
        basename = os.path.basename(bam)
        matched = False
        for s in samples:
            # Simple match: sample id is at start or enclosed in . or _
            # Or just check if sample id in string. 
            # Risk: sample "A" matches "SampleA".
            # Better: Nextflow usually passes files with predictable names if we set it up.
            # Assuming BAM is named "<sample>.bam" or "<sample>.*.bam"
            if s in basename: 
                 sample_bams[s] = bam
                 matched = True
                 break
        if not matched:
            print(f"Warning: Could not match BAM {bam} to any sample in samplesheet.")

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
        
        # Get BAMs for c1
        b1 = [sample_bams[s] for s in samples if samples[s] == c1 and s in sample_bams]
        # Get BAMs for c2
        b2 = [sample_bams[s] for s in samples if samples[s] == c2 and s in sample_bams]
        
        if not b1 or not b2:
            print(f"Skipping contrast {c1} vs {c2} due to missing BAMs.")
            continue
            
        # Write temporary b1.txt and b2.txt
        out_subdir = os.path.join(base_output_dir, f"{c1}_vs_{c2}")
        if not os.path.exists(out_subdir):
            os.makedirs(out_subdir)
            
        b1_path = os.path.join(out_subdir, "b1.txt")
        b2_path = os.path.join(out_subdir, "b2.txt")
        
        with open(b1_path, 'w') as f: f.write(",".join(b1))
        with open(b2_path, 'w') as f: f.write(",".join(b2))
        
        # Construct rMATS command
        # rmats.py --b1 b1.txt --b2 b2.txt --gtf gtf -t paired --readLength 100 --nthread 4 --od out_dir --tmp tmp_dir
        cmd = [
            "rmats.py",
            "--b1", b1_path,
            "--b2", b2_path,
            "--gtf", gtf,
            "-t", read_type,
            "--readLength", str(read_length),
            "--nthread", "4",
            "--od", out_subdir,
            "--tmp", os.path.join(out_subdir, "tmp")
        ]
        
        # Run
        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
