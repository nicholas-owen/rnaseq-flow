#!/usr/bin/env python3

import sys
import gzip

def _open(path):
    """Open a GTF transparently, whether or not it is gzip-compressed."""
    with open(path, 'rb') as probe:
        is_gz = probe.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if is_gz else open(path, 'r')

def parse_gtf(gtf_file):
    transcripts = {}
    with _open(gtf_file) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            if parts[2] != 'exon': continue
            
            chrom = parts[0]
            start = int(parts[3]) - 1 # 0-based
            end = int(parts[4])
            strand = parts[6]
            
            attributes = parts[8]
            attr_dict = {}
            for attr in attributes.split(';'):
                if not attr.strip(): continue
                key_val = attr.strip().split(' ')
                if len(key_val) >= 2:
                    key = key_val[0]
                    val = key_val[1].strip('"')
                    attr_dict[key] = val
            
            tid = attr_dict.get('transcript_id')
            gid = attr_dict.get('gene_id')
            
            if tid:
                if tid not in transcripts:
                    transcripts[tid] = {
                        'chrom': chrom,
                        'strand': strand,
                        'exons': [],
                        'test_gene_id': gid
                    }
                transcripts[tid]['exons'].append((start, end))

    return transcripts

def write_bed12(transcripts):
    for tid, info in transcripts.items():
        chrom = info['chrom']
        strand = info['strand']
        exons = sorted(info['exons'])
        
        if not exons: continue
        
        tx_start = exons[0][0]
        tx_end = exons[-1][1]
        
        block_count = len(exons)
        block_sizes = []
        block_starts = []
        
        for start, end in exons:
            block_sizes.append(str(end - start))
            block_starts.append(str(start - tx_start))
            
        print(f"{chrom}\t{tx_start}\t{tx_end}\t{tid}\t0\t{strand}\t{tx_start}\t{tx_end}\t0\t{block_count}\t{','.join(block_sizes)},\t{','.join(block_starts)},")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: gtf2bed12.py <gtf_file>")
        sys.exit(1)
        
    data = parse_gtf(sys.argv[1])
    write_bed12(data)
