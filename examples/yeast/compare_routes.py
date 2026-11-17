#!/usr/bin/env python3
"""Count and compare the four aligner routes on the yeast example.

Prints the comparison table in this folder's README, and the concordance
figures that follow it, as markdown ready to paste.

    python3 examples/yeast/compare_routes.py            # from the repo root
    python3 examples/yeast/compare_routes.py --dir .    # results_* live here

Why this exists. The table has been wrong twice: once when the routes ran with
the wrong strandedness, and once when the Salmon and Kallisto rows were left
behind by the 1.8.0 transcriptome change. Both times the numbers were prose
that nothing regenerated. Running this and pasting its output is the whole
update.

Why csv.DictReader and not awk. The yeast annotation contains the gene
`ADE5,7`, whose comma sits inside a quoted field. `awk -F,` splits on it, every
later column shifts by one, and `baseMean` gets read as `log2FoldChange`. That
turned one strongly down-regulated gene into an up-regulated one, in three rows
of 5579, which is exactly how the up/down split came out wrong by one per route
the first time this table was rebuilt.
"""
import argparse
import csv
import glob
import os
import sys

ROUTES = ("star", "hisat2", "salmon", "kallisto")


def read_calls(path, padj_col, lfc_col, padj_cut, lfc_cut):
    """Return (tested, significant gene ids, up count) for one results CSV."""
    tested, sig, up = 0, set(), 0
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            padj = row.get(padj_col, "")
            if padj in ("", "NA"):
                continue
            tested += 1
            lfc = float(row.get(lfc_col, "nan"))
            if float(padj) < padj_cut and abs(lfc) > lfc_cut:
                sig.add(row.get("gene_id") or row.get("GeneID") or row[list(row)[0]])
                if lfc > 0:
                    up += 1
    return tested, sig, up


def collect(base, tool, padj_col, lfc_col, padj_cut, lfc_cut):
    """One tool across every route that has results."""
    out = {}
    for route in ROUTES:
        pattern = os.path.join(base, f"results_{route}", f"{tool}_output",
                               f"{tool}_results_*_vs_REF.csv")
        hits = sorted(glob.glob(pattern))
        if not hits:
            continue
        if len(hits) > 1:
            print(f"note: {route} has {len(hits)} contrasts, using {os.path.basename(hits[0])}",
                  file=sys.stderr)
        out[route] = read_calls(hits[0], padj_col, lfc_col, padj_cut, lfc_cut)
    return out


def overlaps(sets):
    """Smallest pairwise overlap, as a percentage of the smaller set."""
    worst, names = 100.0, None
    items = list(sets.items())
    for i, (a, sa) in enumerate(items):
        for b, sb in items[i + 1:]:
            smaller = min(len(sa), len(sb))
            if not smaller:
                continue
            pct = 100.0 * len(sa & sb) / smaller
            if pct < worst:
                worst, names = pct, (a, b)
    return worst, names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".", help="directory holding results_<route>/")
    ap.add_argument("--padj", type=float, default=0.05)
    ap.add_argument("--lfc", type=float, default=1.0)
    args = ap.parse_args()

    deseq2 = collect(args.dir, "deseq2", "padj", "log2FoldChange", args.padj, args.lfc)
    edger = collect(args.dir, "edger", "FDR", "logFC", args.padj, args.lfc)

    if not deseq2:
        sys.exit(f"no results_<route>/deseq2_output/ found under {os.path.abspath(args.dir)}")
    missing = [r for r in ROUTES if r not in deseq2]
    if missing:
        print(f"note: no results for {', '.join(missing)}; "
              "the concordance figures below cover the rest\n", file=sys.stderr)

    print(f"Thresholds: padj < {args.padj}, |log2FC| > {args.lfc}\n")
    print("| route | genes tested | significant | up | down |")
    print("| --- | --- | --- | --- | --- |")
    for route in ROUTES:
        if route not in deseq2:
            continue
        tested, sig, up = deseq2[route]
        label = {"star": "STAR", "hisat2": "HISAT2",
                 "salmon": "Salmon", "kallisto": "Kallisto"}[route]
        print(f"| {label} | {tested} | {len(sig)} | {up} | {len(sig) - up} |")

    for tool, data in (("DESeq2", deseq2), ("edgeR", edger)):
        if len(data) < 2:
            continue
        sets = {r: v[1] for r, v in data.items()}
        common = set.intersection(*sets.values())
        worst, pair = overlaps(sets)
        counts = ", ".join(str(len(sets[r])) for r in ROUTES if r in sets)
        print(f"\n{tool}: {counts} significant for "
              f"{', '.join(r for r in ROUTES if r in sets)}; "
              f"{len(common)} called by all {len(sets)}; "
              f"worst pairwise overlap {worst:.1f}% ({pair[0]} vs {pair[1]})")


if __name__ == "__main__":
    main()
