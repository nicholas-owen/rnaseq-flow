#!/usr/bin/env python3
"""Extract a gene_id -> gene_name / gene_biotype table from a GTF annotation.

The GTF is streamed line by line (so a multi-GB annotation needs only constant
memory) and one row per gene_id is written to `gene_info.tsv`:

    gene_id <tab> gene_name <tab> gene_biotype

This small table is used to annotate the differential-expression result tables
with readable gene symbols, so the outputs do not need a separate ID-mapping
step. A gzipped GTF is handled transparently.

Usage:  gtf2geneinfo.py <annotation.gtf[.gz]>
"""
import sys
import re
import gzip


def open_gtf(path):
    """Open a plain or gzipped GTF for text reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def attr(name, text):
    """Return the value of a GTF attribute (`name "value"`), or '' if absent."""
    m = re.search(name + r' "([^"]*)"', text)
    return m.group(1) if m else ""


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: gtf2geneinfo.py <annotation.gtf[.gz]>")

    gtf = sys.argv[1]
    genes = {}   # gene_id -> (gene_name, gene_biotype); first occurrence wins

    with open_gtf(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = cols[8]
            gid = attr("gene_id", attrs)
            if not gid or gid in genes:
                continue
            # gene_name is optional; fall back to the ID so a name always exists.
            name = attr("gene_name", attrs) or gid
            # Ensembl GTFs use gene_biotype; GENCODE uses gene_type.
            biotype = attr("gene_biotype", attrs) or attr("gene_type", attrs) or "NA"
            genes[gid] = (name, biotype)

    with open("gene_info.tsv", "w") as out:
        out.write("gene_id\tgene_name\tgene_biotype\n")
        for gid, (name, biotype) in genes.items():
            out.write(f"{gid}\t{name}\t{biotype}\n")

    sys.stderr.write(f"gtf2geneinfo: wrote {len(genes)} genes to gene_info.tsv\n")


if __name__ == "__main__":
    main()
