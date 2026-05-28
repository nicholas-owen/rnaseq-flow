#!/usr/bin/env python3

import sys
import os
import re
import requests
import datetime

# Usage: download_refs.py <species> <source> <outdir> [release]
#
#   release : an Ensembl release number (e.g. 102) to pin for a reproducible
#             download, or 'current' (default) for the latest release.
#             A leading 'v' is accepted and stripped (v102 -> 102).

def main():
    if len(sys.argv) < 4:
        print("Usage: download_refs.py <species> <source> <outdir> [release]")
        sys.exit(1)

    species = sys.argv[1].lower().replace(" ", "_")
    source  = sys.argv[2].lower()
    outdir  = sys.argv[3]
    release = sys.argv[4] if len(sys.argv) > 4 else "current"
    release = re.sub(r'^[vV]', '', release.strip()) or "current"

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    log_file = os.path.join(outdir, "download_log.txt")
    with open(log_file, "a") as log:
        log.write(f"Download started at: {datetime.datetime.now()}\n")
        log.write(f"Species: {species}\n")
        log.write(f"Source: {source}\n")
        log.write(f"Ensembl release: {release}\n")

    if source == "ensembl":
        download_ensembl(species, outdir, log_file, release)
    elif source == "ncbi":
        # Robust NCBI downloads need external tooling (ncbi-datasets / datasets).
        print("NCBI download via this script is experimental. Suggest using Ensembl for standard RNA-seq.")
        with open(log_file, "a") as log:
            log.write("Source 'ncbi' requested. Implementation limited in this script version.\n")
        download_ncbi_experimental(species, outdir, log_file)
    else:
        print(f"Unknown source: {source}. Supported: ensembl")
        sys.exit(1)


def download_file(url, outfile, log_file):
    print(f"Downloading {url}...")
    with open(log_file, "a") as log:
        log.write(f"Downloading: {url}\n")
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(outfile, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    print(f"Saved to {outfile}")


def download_ensembl(species, outdir, log_file, release):
    # Choose between the rolling 'current' release and a pinned release-<N>
    # directory. Pinning a release makes the download reproducible.
    if release in ("", "current"):
        base_ftp_dna = f"http://ftp.ensembl.org/pub/current_fasta/{species}/dna/"
        base_ftp_gtf = f"http://ftp.ensembl.org/pub/current_gtf/{species}/"
        print("Using the current (latest) Ensembl release.")
    else:
        base_ftp_dna = f"http://ftp.ensembl.org/pub/release-{release}/fasta/{species}/dna/"
        base_ftp_gtf = f"http://ftp.ensembl.org/pub/release-{release}/gtf/{species}/"
        print(f"Using pinned Ensembl release-{release} (reproducible).")

    # Informational only: report the latest assembly name. This call is
    # best-effort and never blocks the download -- the file finder below works
    # by scraping the directory listing and is assembly-agnostic.
    try:
        r = requests.get(f"http://rest.ensembl.org/info/assembly/{species}?",
                         headers={"Content-Type": "application/json"}, timeout=30)
        if r.ok:
            print(f"Latest assembly for {species}: {r.json().get('default_coord_system_version')}")
        else:
            print(f"(Ensembl REST did not recognise species '{species}'; continuing via FTP listing.)")
    except Exception as e:
        print(f"(Could not query Ensembl REST for assembly info: {e})")

    genus_species = species.capitalize()

    # ---- genome FASTA ----
    print("Finding genome FASTA...")
    try:
        r_dna = requests.get(base_ftp_dna)
        candidates = re.findall(rf'href="({genus_species}.*?\.dna\.primary_assembly\.fa\.gz)"', r_dna.text)
        if not candidates:
            candidates = re.findall(rf'href="({genus_species}.*?\.dna\.toplevel\.fa\.gz)"', r_dna.text)
        if candidates:
            fn = candidates[0]
            download_file(base_ftp_dna + fn, os.path.join(outdir, fn), log_file)
        else:
            msg = f"Could not locate a genome FASTA at {base_ftp_dna}"
            print(msg)
            with open(log_file, "a") as log: log.write(f"Error: {msg}\n")
    except Exception as e:
        print(f"Error fetching FASTA: {e}")
        with open(log_file, "a") as log: log.write(f"Error fetching FASTA: {e}\n")

    # ---- GTF annotation ----
    print("Finding GTF annotation...")
    try:
        r_gtf = requests.get(base_ftp_gtf)
        candidates_gtf = re.findall(rf'href="({genus_species}.*?\.gtf\.gz)"', r_gtf.text)
        valid_gtfs = [c for c in candidates_gtf if "abinitio" not in c and "chr.gtf" not in c]
        if valid_gtfs:
            fn = valid_gtfs[0]
            download_file(base_ftp_gtf + fn, os.path.join(outdir, fn), log_file)
        else:
            msg = f"Could not locate a GTF at {base_ftp_gtf}"
            print(msg)
            with open(log_file, "a") as log: log.write(f"Error: {msg}\n")
    except Exception as e:
        print(f"Error fetching GTF: {e}")
        with open(log_file, "a") as log: log.write(f"Error fetching GTF: {e}\n")


def download_ncbi_experimental(species, outdir, log_file):
    print("NCBI implementation requires external tools for reliability. Skipped.")


if __name__ == "__main__":
    main()
