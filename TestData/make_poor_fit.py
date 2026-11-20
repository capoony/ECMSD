#!/usr/bin/env python3
"""Make TestData/merged_poor_fit.fastq.gz from merged.fastq.gz.

Takes the first N reads and substitutes RATE of the bases (seeded), so the reads
fit their true reference poorly, like a haplotype missing from the database.
Usage: make_poor_fit.py merged.fastq.gz merged_poor_fit.fastq.gz [N=200000] [RATE=0.06]
"""
import gzip, random, sys

src, dst = sys.argv[1], sys.argv[2]
n = int(sys.argv[3]) if len(sys.argv) > 3 else 200000
rate = float(sys.argv[4]) if len(sys.argv) > 4 else 0.06
rng = random.Random(1)
with gzip.open(src, "rt") as f, gzip.open(dst, "wt") as o:
    for _ in range(n):
        h, s, p, q = (f.readline() for _ in range(4))
        if not q:
            break
        s = "".join(rng.choice([b for b in "ACGT" if b != c]) if c in "ACGT" and rng.random() < rate else c
                    for c in s.rstrip("\n"))
        o.write(f"{h}{s}\n{p}{q}")
