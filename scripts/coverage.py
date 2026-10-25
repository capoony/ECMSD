#!/usr/bin/env python3
"""Per-reference coverage statistics from a (gzipped) PAF file."""
import gzip
import math
import sys
from collections import defaultdict


def coverage(paf_file, out_file, min_mapq):
    diff_arrays = defaultdict(lambda: defaultdict(int))
    ref_lengths = {}
    # Dedupe on the full interval, not just (read, ref): mates share a read name
    # but are distinct alignments and must both count towards breadth.
    seen = set()

    opener = gzip.open if paf_file.endswith(".gz") else open
    with opener(paf_file, "rt") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            read_name, ref_name = parts[0], parts[5]
            ref_len, ref_start, ref_end = int(parts[6]), int(parts[7]), int(parts[8])
            ref_lengths[ref_name] = ref_len
            if int(parts[11]) < min_mapq:
                continue
            key = (read_name, ref_name, ref_start, ref_end)
            if key in seen:
                continue
            seen.add(key)
            diff_arrays[ref_name][ref_start] += 1
            diff_arrays[ref_name][ref_end] -= 1

    with open(out_file, "w") as out:
        out.write("reference\tref_length\tmean_coverage\tstd_coverage\tpct_covered\n")
        for ref_name, ref_len in sorted(ref_lengths.items()):
            diff = diff_arrays[ref_name]
            depth = covered = sum_d = sum_sq_d = 0
            for pos in range(ref_len):
                depth += diff.get(pos, 0)
                sum_d += depth
                sum_sq_d += depth * depth
                if depth > 0:
                    covered += 1
            mean_cov = sum_d / ref_len if ref_len > 0 else 0
            variance = (sum_sq_d / ref_len - mean_cov ** 2) if ref_len > 0 else 0
            std_cov = math.sqrt(max(variance, 0.0))
            pct_covered = (covered / ref_len * 100) if ref_len > 0 else 0
            out.write(f"{ref_name}\t{ref_len}\t{mean_cov:.2f}\t{std_cov:.2f}\t{pct_covered:.2f}\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: coverage.py <PAF[.gz]> <OUT> <MIN_MAPQ>")
    coverage(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    print(f"Coverage statistics written to: {sys.argv[2]}")
