#!/usr/bin/env python3
"""Turn the reported references into calls with a status (Mito_calls.txt).

Inputs: the evidence table of assess_references.py (columns call_status, call_group,
call_reason) and the read assignment of LinkTaxonomy.py (<prefix>_Mito_summary.txt, column
Weight). One row per call, sorted by proportion:
  confident     the reads fit this reference as expected for their own species and separate it
                from the other credible references
  ambiguous     references that share most of their compatible reads with each other (one
                row per group, e.g. 'Drosophila group: D. mauritiana / D. simulans'); the reads
                cannot tell which of them is present
  poor_fit      the reads fit the reference worse than reads fit their own species (or only one
                region of it fits): a relative of the source, not the source. A group whose
                members all fit poorly is one poor_fit row.
  low_evidence  credible only through reads at the minimap2 -N limit
When every call is poor_fit and all are in one genus, the sample line says so: the reads come
from an unsampled member of that genus (or a haplotype not in the database), and the named
species are only its closest relatives in the database.
"""
import argparse
import csv
import sys
from collections import defaultdict

HEADER = ["call", "status", "proportion", "weighted_reads", "taxids", "members", "reason"]


def load_calls(path):
    """reference -> (name, status, group, reason) for the credible references."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("credible") == "yes":
                out[r["reference"]] = (r["taxon_name"].replace("_", " "), r["call_status"],
                                       r["call_group"], r["call_reason"])
    return out


def load_weights(path):
    """reference -> (summed Weight, genus) from the Mito_summary.txt."""
    w, genus = defaultdict(float), {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            w[r["TaxID"]] += float(r.get("Weight") or 1)
            genus[r["TaxID"]] = r.get("genus", "NA")
    return w, genus


def abbrev(name):
    p = name.split()
    return f"{p[0][0]}. {' '.join(p[1:])}" if len(p) > 1 else name


def make_calls(calls, weights, genus):
    """List of call rows and the sample verdict line."""
    groups = defaultdict(list)
    for ref in weights:
        name, status, group, reason = calls.get(ref, (ref, "confident", "-", "-"))
        groups[group if group != "-" else "ref:" + ref].append((ref, name, status, reason))
    total = sum(weights.values())
    rows = []
    for key, ms in groups.items():
        ms.sort(key=lambda m: -weights[m[0]])
        w = sum(weights[m[0]] for m in ms)
        if len(ms) == 1:
            ref, name, status, reason = ms[0]
            call = name
        else:
            gen = {genus.get(m[0]) for m in ms}
            head = gen.pop() + " group" if len(gen) == 1 else "group"
            call = f"{head} ({key}): " + " / ".join(abbrev(m[1]) for m in ms)
            status = "poor_fit" if all(m[2] == "poor_fit" for m in ms) else "ambiguous"
            reason = "; ".join(f"{abbrev(m[1])}: {m[2]} {m[3]}" for m in ms)
        rows.append([call, status, w / total if total else 0.0, w, ",".join(m[0] for m in ms),
                     "|".join(m[1] for m in ms), reason])
    rows.sort(key=lambda r: -r[2])
    gens = {genus.get(ref) for ref in weights}
    if not rows:
        verdict = "no species reported"
    elif all(r[1] == "poor_fit" for r in rows) and len(gens) == 1 and None not in gens and "NA" not in gens:
        verdict = (f"reads fit no reference well: an unsampled {gens.pop()} species or haplotype; "
                   "the calls are its closest relatives in the database")
    elif all(r[1] == "confident" for r in rows):
        verdict = "all calls confident"
    else:
        verdict = "some calls are not confident, see status"
    return rows, verdict


def db_line(path):
    """One-line version of the database from its db_info.txt."""
    info = dict(l.rstrip("\n").split("\t", 1) for l in open(path) if "\t" in l)
    line = (f"RefSeq release {info.get('refseq_release', 'unknown')}, built {info.get('built', 'unknown')}, "
            f"taxonomy {info.get('taxdump_downloaded', 'unknown')}")
    if info.get("n_extra_references", "0") != "0":
        line += f", {info['n_extra_references']} extra references (extended {info.get('extended', 'unknown')})"
    return line


def report(rows, verdict):
    """Terminal summary: the sample verdict and one line per call."""
    lines = [f"Sample: {verdict}"]
    for r in rows:
        lines.append(f"  {r[1]:<13} {r[2]:6.3f}  {r[0]}")
    if any(r[1] != "confident" for r in rows):
        lines.append("  Only confident rows are species calls; see the reason column of Mito_calls.txt.")
    return "\n".join(lines)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--Evidence", required=True, help="evidence table of assess_references.py")
    p.add_argument("--Summary", required=True, help="<prefix>_Mito_summary.txt of LinkTaxonomy.py")
    p.add_argument("--output", required=True, help="calls table (TSV)")
    p.add_argument("--DbInfo", help="db_info.txt of the database, summarised in the header")
    a = p.parse_args(argv)
    weights, genus = load_weights(a.Summary)
    rows, verdict = make_calls(load_calls(a.Evidence), weights, genus)
    with open(a.output, "w") as out:
        out.write("# sample: " + verdict + "\n")
        if a.DbInfo:
            out.write("# db: " + db_line(a.DbInfo) + "\n")
        out.write("\t".join(HEADER) + "\n")
        for r in rows:
            out.write("\t".join(r[:2] + ["%.4f" % r[2], "%.1f" % r[3]] + r[4:]) + "\n")
    print(report(rows, verdict))
    print(f"Calls written to: {a.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
