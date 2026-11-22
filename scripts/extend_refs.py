#!/usr/bin/env python3
"""Prepare extra mitogenomes for the ECMSD database (ECMSD --extend-db).

Input is either a FASTA file with headers <taxid>|<id> or a list of NCBI nucleotide
accessions (one per line), whose sequences and taxids are fetched with E-utilities.
Every record is checked against the full NCBI taxonomy of the database:
  - taxids merged into another taxid (merged.dmp) are replaced by the current one;
  - taxids not in the taxonomy are rejected;
  - ids already in the existing extra references or repeated in the input are skipped;
  - taxids that already have references (RefSeq or earlier extensions) are kept with a note;
  - records shorter than 5 kb and RefSeq NC_ accessions (probably already in the
    database) are kept with a warning.
Writes the accepted records with headers <taxid>|<id> (unmasked, for bbmask) and a
TSV log (id, taxid, length, note).
"""
import argparse
import gzip
import json
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
HEADER = re.compile(r"^(\d+)\|(\S+)$")
MIN_LEN = 5000
BATCH = 200


def opener(path):
    """Open plain or gzipped text; gzip is detected by content (temp names end in .tmp)."""
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", encoding="latin-1") if gz else open(path, encoding="latin-1")


def read_fasta(fh):
    """Yield (header first token, sequence)."""
    name, seq = None, []
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0] if len(line) > 1 else "", []
        elif line:
            seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def load_taxonomy(nodes, merged):
    """(set of taxids, old taxid -> current taxid)."""
    with opener(nodes) as fh:
        known = {line.split("\t|\t", 1)[0].strip() for line in fh}
    remap = {}
    if merged:
        with opener(merged) as fh:
            for line in fh:
                f = line.split("\t|\t")
                remap[f[0].strip()] = f[1].strip().rstrip("\t|")
    return known, remap


def existing_ids(path):
    if not path:
        return set()
    with opener(path) as fh:
        return {line[1:].split()[0].split("|", 1)[-1] for line in fh if line.startswith(">")}


def taxid_counts(*paths):
    """taxid -> number of references with that taxid in the given <taxid>|<id> FASTA files."""
    counts = {}
    for path in filter(None, paths):
        with opener(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    t = line[1:].split("|", 1)[0]
                    counts[t] = counts.get(t, 0) + 1
    return counts


def post(tool, **params):
    data = urllib.parse.urlencode(dict(db="nuccore", tool="ECMSD", **params)).encode()
    for attempt in range(3):
        try:
            with urllib.request.urlopen(EUTILS + tool, data=data, timeout=120) as r:
                return r.read().decode()
        except OSError as e:
            if attempt == 2:
                raise SystemExit(f"Error: NCBI {tool} failed: {e}")
            time.sleep(2)
    return ""


def fetch_accessions(accs):
    """Yield ('<taxid>|<accession.version>', sequence) for NCBI accessions."""
    for i in range(0, len(accs), BATCH):
        chunk = accs[i:i + BATCH]
        result = json.loads(post("esummary.fcgi", retmode="json", id=",".join(chunk)))["result"]
        taxid = {}
        for uid in result.get("uids", []):
            r = result[uid]
            taxid[r["accessionversion"]] = str(r["taxid"])
        found = {a.split(".")[0] for a in taxid}
        for a in chunk:
            if a.split(".")[0] not in found:
                print(f"Warning: {a} not found at NCBI, skipped.", file=sys.stderr)
        time.sleep(0.4)  # NCBI allows 3 requests per second without an API key
        fasta = post("efetch.fcgi", rettype="fasta", retmode="text", id=",".join(taxid))
        for name, seq in read_fasta(fasta.splitlines()):
            if name in taxid:
                yield f"{taxid[name]}|{name}", seq
        time.sleep(0.4)


def check(records, known, remap, skip_ids, have=None):
    """Yield (taxid, id, seq, note) for accepted records; print why others are dropped."""
    seen = set(skip_ids)
    for header, seq in records:
        m = HEADER.match(header)
        if not m:
            print(f"Rejected {header}: header must be <taxid>|<id>.", file=sys.stderr)
            continue
        taxid, rid = m.groups()
        notes = []
        if taxid in remap:
            notes.append(f"taxid {taxid} merged into {remap[taxid]}")
            taxid = remap[taxid]
        if taxid not in known:
            print(f"Rejected {rid}: taxid {taxid} is not in the database taxonomy "
                  "(newer than the stored taxdump? run ECMSD --update-db).", file=sys.stderr)
            continue
        if rid in seen:
            print(f"Skipped {rid}: already in the database or listed twice.", file=sys.stderr)
            continue
        seen.add(rid)
        if not seq:
            print(f"Rejected {rid}: empty sequence.", file=sys.stderr)
            continue
        if have and have.get(taxid):
            notes.append(f"taxid {taxid} already has {have[taxid]} reference(s) in the database")
        if len(seq) < MIN_LEN:
            notes.append(f"only {len(seq)} bp, probably not a complete mitogenome")
        if rid.startswith("NC_"):
            notes.append("RefSeq NC_ accession, probably already in the database")
        for n in notes:
            print(f"Warning {rid}: {n}.", file=sys.stderr)
        yield taxid, rid, seq, "; ".join(notes) or "-"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--fasta", help="FASTA with headers <taxid>|<id>")
    src.add_argument("--accessions", help="file with NCBI nucleotide accessions, one per line")
    p.add_argument("--nodes", required=True, help="full nodes.dmp(.gz)")
    p.add_argument("--merged", help="merged.dmp(.gz)")
    p.add_argument("--existing", help="existing extra references (skip their ids)")
    p.add_argument("--refs", help="main reference FASTA (to note taxids that are already in)")
    p.add_argument("--output", required=True, help="accepted records (FASTA)")
    p.add_argument("--log", required=True, help="accepted records (TSV: id, taxid, length, note)")
    a = p.parse_args(argv)
    known, remap = load_taxonomy(a.nodes, a.merged)
    if a.fasta:
        with opener(a.fasta) as fh:
            records = list(read_fasta(fh))
    else:
        with open(a.accessions) as fh:
            accs = list(dict.fromkeys(l.split()[0] for l in fh if l.strip() and not l.startswith("#")))
        records = list(fetch_accessions(accs))
    n = 0
    with open(a.output, "w") as out, open(a.log, "w") as log:
        for taxid, rid, seq, note in check(records, known, remap, existing_ids(a.existing),
                                                      taxid_counts(a.refs, a.existing)):
            out.write(f">{taxid}|{rid}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")
            log.write(f"{rid}\t{taxid}\t{len(seq)}\t{note}\n")
            n += 1
    print(f"{n} of {len(records)} records accepted.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
