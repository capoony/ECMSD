#!/usr/bin/env python3
"""Extract Stage 1 candidate reads from the original FASTQ file(s).

Stage 1 of ECMSD maps all reads liberally against the mitochondrial database
and writes the names of every read with at least one reported alignment (any
MAPQ, including 0). This script streams the original FASTQ file(s) and writes
the records of those reads, unchanged (header, sequence and qualities), to new
FASTQ file(s) for the Stage 2 mapping.

Single-end / merged reads: one input, one output.
Paired-end reads: two inputs, two outputs. A pair is kept, with both mates and
in the original order, when either mate had a Stage 1 hit, so that Stage 2 can
still be run in paired mode.

Memory: only the set of hit names is held in memory; FASTQ records are streamed.
Gzip input is detected from the file content, gzip output from the file name.
"""
import argparse
import gzip
import re
import sys

# minimap2 appends /1 and /2 to mate names in paired mode (even if the FASTQ
# name already ends in /1 or /2), so strip any trailing mate suffixes on both
# the hit names and the FASTQ names before comparing.
_MATE_SUFFIX = re.compile(rb"(?:/[12])+$")


def normalise(name):
    return _MATE_SUFFIX.sub(b"", name)


def open_binary(path, mode="rb"):
    """Open a plain or gzipped file; gzip input is recognised by its magic bytes."""
    if mode == "rb":
        with open(path, "rb") as fh:
            magic = fh.read(2)
        return gzip.open(path, "rb") if magic == b"\x1f\x8b" else open(path, "rb")
    return gzip.open(path, "wb", compresslevel=1) if path.endswith(".gz") else open(path, "wb")


def load_hit_names(path):
    """Set of normalised read names, one name per line (tab-separated extra columns ignored)."""
    names = set()
    with open_binary(path) as fh:
        for line in fh:
            name = line.rstrip(b"\r\n").split(b"\t", 1)[0]
            if name:
                names.add(normalise(name))
    return names


def read_name(header):
    """Read name from a FASTQ header line: first word, without the leading '@'."""
    return header[1:].split(None, 1)[0] if len(header) > 1 else b""


def fastq_records(fh, path="<fastq>"):
    """Yield (name, record_bytes) for each 4-line FASTQ record."""
    while True:
        header = fh.readline()
        if not header:
            return
        if not header.strip():
            continue
        rest = [fh.readline() for _ in range(3)]
        if not header.startswith(b"@") or not rest[1].startswith(b"+") or not rest[2]:
            raise ValueError(f"{path}: malformed FASTQ record starting with {header[:50]!r}")
        if not header.endswith(b"\n"):
            header += b"\n"
        rest[2] = rest[2] if rest[2].endswith(b"\n") else rest[2] + b"\n"
        yield normalise(read_name(header.rstrip(b"\r\n"))), header + rest[0] + rest[1] + rest[2]


def extract_single(hits, fastq_in, fastq_out):
    n_in = n_out = 0
    with open_binary(fastq_in) as fin, open_binary(fastq_out, "wb") as fout:
        for name, record in fastq_records(fin, fastq_in):
            n_in += 1
            if name in hits:
                fout.write(record)
                n_out += 1
    return n_in, n_out


def extract_paired(hits, in1, in2, out1, out2):
    n_in = n_out = 0
    with open_binary(in1) as f1, open_binary(in2) as f2, \
            open_binary(out1, "wb") as o1, open_binary(out2, "wb") as o2:
        r1, r2 = fastq_records(f1, in1), fastq_records(f2, in2)
        for rec1 in r1:
            rec2 = next(r2, None)
            if rec2 is None:
                raise ValueError(f"{in2} has fewer reads than {in1}")
            if rec1[0] != rec2[0]:
                raise ValueError(
                    f"mate names differ at pair {n_in + 1}: {rec1[0].decode()} vs {rec2[0].decode()}; "
                    "paired FASTQ files must be in the same order")
            n_in += 1
            if rec1[0] in hits:
                o1.write(rec1[1])
                o2.write(rec2[1])
                n_out += 1
        if next(r2, None) is not None:
            raise ValueError(f"{in2} has more reads than {in1}")
    return n_in, n_out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--hits", required=True,
                   help="file with the names of reads that had a Stage 1 alignment (one per line; .gz allowed)")
    p.add_argument("--in", dest="fastq_in", nargs="+", required=True, metavar="FASTQ",
                   help="original FASTQ file (single-end/merged) or R1 R2 (paired-end)")
    p.add_argument("--out", dest="fastq_out", nargs="+", required=True, metavar="FASTQ",
                   help="output FASTQ file(s), same number as --in; .gz extension gives gzip output")
    args = p.parse_args(argv)

    if len(args.fastq_in) not in (1, 2) or len(args.fastq_in) != len(args.fastq_out):
        p.error("--in and --out need the same number of files (1 or 2)")

    hits = load_hit_names(args.hits)
    try:
        if len(args.fastq_in) == 1:
            n_in, n_out = extract_single(hits, args.fastq_in[0], args.fastq_out[0])
            unit = "reads"
        else:
            n_in, n_out = extract_paired(hits, *args.fastq_in, *args.fastq_out)
            unit = "read pairs"
    except ValueError as err:
        sys.exit(f"extract_candidates.py: {err}")
    print(f"Stage 1 hit names: {len(hits)}; candidate {unit} extracted: {n_out} of {n_in}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
