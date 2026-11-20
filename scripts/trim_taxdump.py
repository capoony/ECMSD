#!/usr/bin/env python3
"""Trim nodes.dmp and names.dmp to the lineages of the database taxids.

Reads the full taxonomy (plain or gzipped nodes.dmp/names.dmp) and writes only the taxids
of the reference FASTA files plus every ancestor up to root 1, so LinkTaxonomy.py lookups
still resolve. FASTA headers are <taxid> or <taxid>|<accession>. Outputs are written as
.tmp and moved into place, so a running analysis never reads a half-written file.
Taxids missing from the taxonomy are printed to stdout, one per line.
"""
import argparse
import gzip
import os
import sys


def opener(path):
    """Open plain or gzipped text; gzip is detected by content (temp names end in .tmp)."""
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", encoding="latin-1") if gz else open(path, encoding="latin-1")


def fasta_taxids(paths):
    ids = set()
    for p in paths:
        with opener(p) as fh:
            for line in fh:
                if line.startswith(">"):
                    ids.add(line[1:].split()[0].split("|")[0])
    return ids


def lineages(nodes_fh, taxids):
    """(taxids with all ancestors, taxids not in the taxonomy)."""
    parent = {}
    for line in nodes_fh:
        f = line.split("\t|\t")
        parent[f[0].strip()] = f[1].strip()
    keep, todo = set(), list(taxids)
    while todo:
        t = todo.pop()
        if t in keep or t not in parent:
            continue
        keep.add(t)
        todo.append(parent[t])
    return keep, sorted(t for t in taxids if t not in parent)


def write_subset(src, dest, keep):
    with opener(src) as fh, open(dest + ".tmp", "w", encoding="latin-1") as out:
        for line in fh:
            if line.split("\t|\t")[0].strip() in keep:
                out.write(line)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--nodes", required=True, help="full nodes.dmp(.gz)")
    p.add_argument("--names", required=True, help="full names.dmp(.gz)")
    p.add_argument("--out-nodes", required=True)
    p.add_argument("--out-names", required=True)
    p.add_argument("fasta", nargs="+", help="reference FASTA files (plain or gzipped)")
    a = p.parse_args(argv)
    taxids = fasta_taxids(a.fasta)
    if not taxids:
        sys.exit("Error: no taxids found in " + ", ".join(a.fasta))
    with opener(a.nodes) as fh:
        keep, missing = lineages(fh, taxids)
    write_subset(a.nodes, a.out_nodes, keep)
    write_subset(a.names, a.out_names, keep)
    os.replace(a.out_nodes + ".tmp", a.out_nodes)
    os.replace(a.out_names + ".tmp", a.out_names)
    for t in missing:
        print(t)
    print(f"Taxonomy trimmed to {len(keep)} taxa for {len(taxids)} database taxids "
          f"({len(missing)} not in the taxonomy).", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
