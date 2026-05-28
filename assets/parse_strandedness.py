#!/usr/bin/env python3
"""Reduce an RSeQC infer_experiment.py report to a single strandedness verdict.

Usage: parse_strandedness.py <infer_experiment.txt>

Prints one of: forward | reverse | unstranded

RSeQC reports two 'explained by' fractions, in a fixed order: the first is the
forward/sense fraction, the second the reverse/antisense fraction. If one
clearly dominates the library is stranded that way; otherwise it is unstranded.
"""
import sys

if len(sys.argv) < 2:
    sys.exit("Usage: parse_strandedness.py <infer_experiment.txt>")

fracs = []
with open(sys.argv[1]) as fh:
    for line in fh:
        if "explained by" in line:
            try:
                fracs.append(float(line.rsplit(":", 1)[-1]))
            except ValueError:
                pass

if len(fracs) >= 2 and (fracs[0] + fracs[1]) > 0:
    forward_fraction = fracs[0] / (fracs[0] + fracs[1])
    if forward_fraction >= 0.8:
        verdict = "forward"
    elif forward_fraction <= 0.2:
        verdict = "reverse"
    else:
        verdict = "unstranded"
else:
    # Inference inconclusive or report unreadable -> safest default.
    verdict = "unstranded"

print(verdict)
