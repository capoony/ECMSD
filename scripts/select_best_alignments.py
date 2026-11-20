#!/usr/bin/env python3
"""Select the best-scoring ("compatible") alignments of each read from the ECMSD PAF.

ECMSD maps the reads with `minimap2 -c --secondary=yes -N <n>`, so a
read can have several alignments to different references. They are compared by
alignment score (AS), never by MAPQ: closely related mitogenomes give a correct
alignment a low MAPQ (even 0), and MAPQ is not used to rank anything.

For every read (paired mates are separate reads, minimap2 names them NAME/1 and
NAME/2) the script keeps, for each reference, its best alignment (higher AS, then
primary over secondary, then lower NM, then first reported) and writes every
alignment whose AS equals the top AS of the read: the references the read is
compatible with. Ties are kept so that coverage can use them and so that the
assignment step can see which reads cannot tell references apart. Two tags are
appended to each written line:
  nc:i:N   number of references the read is compatible with (1 = unique best)
  dm:i:D   top AS minus the best AS on any other reference, only present when a
           lower-scoring competing reference was reported
  dt:i:D   like dm, but only over references of another taxid than all compatible ones
           (reference names "<taxid>|<accession>"; a bare name is its own taxid). Only
           written when it differs from dm, i.e. when several references of one taxid
           (haplotypes of a species) are in the database; -1 = no lower-scoring
           reference of another taxid was reported
  na:i:A   number of alignments minimap2 reported for the read (all references). When
           this reaches the -N limit the list of competitors is truncated, so a unique
           best alignment is only provisional
With --one-per-read only the single best alignment is written, without the tags.

All alignments of a read must be adjacent in the input, which is how minimap2 writes
them. Input and output may be '-' (stdin / stdout), so minimap2 can be piped in
directly and the (very large) all-alignments PAF never has to be stored.

Optional read filters, judged on the best alignment of the read, remove reads that
have no credible alignment, such as short chance matches: --min-mapq (highest MAPQ
of the compatible alignments), --min-aln-len, --min-aln-frac (aligned fraction of
the read) and --min-identity (1 - NM / alignment length).

Tags: tp, AS and NM are required for a reliable choice. Alignments missing AS are
ranked below all alignments that have it, a read without any AS falls back to its
primary (or first) alignment, and every such case is counted and reported. If no
alignment in the whole file carries AS (e.g. minimap2 was run without -c) the
script stops with an error.

Optionally every alignment is written to a TSV table (--table) with its tags, rank
and score difference to the best alignment of the read. This is slower and the
table is large, so it is only written on request.
"""
import argparse
import gzip
import sys

TABLE_HEADER = ("read\tref\tref_len\tref_start\tref_end\tquery_len\tquery_start\tquery_end\t"
                "strand\taligned_len\tmapq\ttp\tAS\tNM\tbest\tAS_delta\tn_alignments\n")
INF = float("inf")


def opener(path, mode):
    if path == "-":
        return sys.stdin if mode.startswith("r") else sys.stdout
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def parse_tags(tail):
    """Return (tp, AS, NM) from the optional PAF fields (text after the 12th tab).
    A tag that is absent is None."""
    tail = "\t" + tail
    tp = AS = NM = None
    i = tail.find("\ttp:A:")
    if i >= 0:
        tp = tail[i + 6:i + 7]
    i = tail.find("\tAS:i:")
    if i >= 0:
        j = tail.find("\t", i + 1)
        AS = int(tail[i + 6:] if j < 0 else tail[i + 6:j])
    i = tail.find("\tNM:i:")
    if i >= 0:
        j = tail.find("\t", i + 1)
        NM = int(tail[i + 6:] if j < 0 else tail[i + 6:j])
    return tp, AS, NM


def rank_key(tp, AS, NM):
    """Sort key: smaller is better. Ties are resolved by input order (first wins)."""
    return (0 if AS is not None else 1,
            -AS if AS is not None else 0,
            0 if tp == "P" else 1,
            NM if NM is not None else INF)


def select_best(paf_in, paf_out, table_out=None, one_per_read=False, min_mapq=0,
                min_aln_len=0, min_aln_frac=0.0, min_identity=0.0):
    stats = dict(lines=0, reads=0, with_AS=0, missing_AS=0, missing_NM=0, missing_tp=0,
                 reads_no_AS=0, reads_multi=0, reads_AS_tie=0, best_secondary=0, bad_lines=0,
                 reads_out=0, lines_out=0, reads_low_mapq=0, reads_short=0,
                 reads_low_fraction=0, reads_low_identity=0)
    tbl = opener(table_out, "wt") if table_out else None
    if tbl:
        tbl.write(TABLE_HEADER)
    check_quality = bool(min_aln_len or min_aln_frac or min_identity)

    # State of the read being read. by_ref maps each reference to its best alignment
    # [rank_key, line, tp, AS, NM, index, mapq, taxid]; the whole group is only kept for the table.
    name = None
    by_ref = {}
    n = 0
    group = []

    def flush(out):
        if n == 0:
            return
        stats["reads"] += 1
        entries = list(by_ref.values())
        scored = [e for e in entries if e[3] is not None]
        if scored:
            top = max(e[3] for e in scored)
            compat = [e for e in entries if e[3] == top]
            lower = [e[3] for e in scored if e[3] < top]
            margin = top - max(lower) if lower else None
            top_taxa = {e[7] for e in compat}
            lower_other = [e[3] for e in scored if e[3] < top and e[7] not in top_taxa]
            margin_taxid = top - max(lower_other) if lower_other else None
        else:
            stats["reads_no_AS"] += 1
            top = margin = margin_taxid = None
            compat = [min(entries, key=lambda e: (e[0], e[5]))]
        best = min(compat, key=lambda e: (e[0], e[5]))
        if n > 1:
            stats["reads_multi"] += 1
        if len(compat) > 1:
            stats["reads_AS_tie"] += 1
        if best[2] == "S":
            stats["best_secondary"] += 1

        keep = True
        if max(e[6] for e in compat) < min_mapq:
            stats["reads_low_mapq"] += 1
            keep = False
        elif check_quality:
            f = best[1].split("\t", 12)
            alen, qlen = int(f[10]), int(f[1])
            frac = (int(f[3]) - int(f[2])) / qlen if qlen > 0 else 0.0
            ident = None if best[4] is None or alen <= 0 else 1.0 - best[4] / alen
            if alen < min_aln_len:
                stats["reads_short"] += 1
                keep = False
            elif frac < min_aln_frac:
                stats["reads_low_fraction"] += 1
                keep = False
            elif ident is not None and ident < min_identity:
                stats["reads_low_identity"] += 1
                keep = False

        written = {best[5]} if one_per_read else {e[5] for e in compat}
        if keep:
            stats["reads_out"] += 1
            if one_per_read:
                out.write(best[1])
                stats["lines_out"] += 1
            else:
                tags = ("\tnc:i:%d" % len(compat) + ("" if margin is None else "\tdm:i:%d" % margin)
                        + ("" if margin_taxid == margin else
                           "\tdt:i:%d" % (-1 if margin_taxid is None else margin_taxid))
                        + "\tna:i:%d" % n)
                for e in compat:
                    out.write(e[1][:-1] + tags + "\n")
                stats["lines_out"] += len(compat)
        if tbl:
            for idx, (line, tp, AS, NM) in enumerate(group):
                f = line.rstrip("\r\n").split("\t")
                delta = "NA" if (AS is None or top is None) else top - AS
                tbl.write("\t".join([
                    f[0], f[5], f[6], f[7], f[8], f[1], f[2], f[3], f[4], f[10], f[11],
                    tp or "NA", "NA" if AS is None else str(AS), "NA" if NM is None else str(NM),
                    "1" if (keep and idx in written) else "0", str(delta), str(n)]) + "\n")

    with opener(paf_in, "rt") as fin, opener(paf_out, "wt") as out:
        for line in fin:
            f = line.split("\t", 12)
            if len(f) < 12:
                continue  # stray log text or blank line
            try:
                int(f[7]), int(f[8])
                mapq = int(f[11])
            except ValueError:
                stats["bad_lines"] += 1
                continue
            if not line.endswith("\n"):
                line += "\n"
            stats["lines"] += 1
            tp, AS, NM = parse_tags(f[12]) if len(f) > 12 else (None, None, None)
            if AS is None:
                stats["missing_AS"] += 1
            else:
                stats["with_AS"] += 1
            if NM is None:
                stats["missing_NM"] += 1
            if tp is None:
                stats["missing_tp"] += 1
            if f[0] != name:
                flush(out)
                name = f[0]
                by_ref = {}
                n = 0
                group = []
            key = rank_key(tp, AS, NM)
            cur = by_ref.get(f[5])
            if cur is None or key < cur[0]:
                by_ref[f[5]] = [key, line, tp, AS, NM, n, mapq, f[5].split("|", 1)[0]]
            if tbl:
                group.append((line, tp, AS, NM))
            n += 1
        flush(out)
    if tbl:
        tbl.close()
    return stats


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("paf_in", help="ECMSD PAF (plain or .gz, or '-' for stdin), made with minimap2 -c")
    p.add_argument("paf_out", help="output PAF (.gz gives gzip, '-' stdout)")
    p.add_argument("--table", help="optional TSV (.gz allowed) listing every alignment with its tags")
    p.add_argument("--one-per-read", action="store_true",
                   help="write only the single best alignment of each read, without nc/dm tags")
    p.add_argument("--min-mapq", type=int, default=0,
                   help="drop reads whose compatible alignments all have a MAPQ below this (default: 0 = off)")
    p.add_argument("--min-aln-len", type=int, default=0,
                   help="drop reads whose best alignment is shorter than this many bases (default: 0 = off)")
    p.add_argument("--min-aln-frac", type=float, default=0.0,
                   help="drop reads of which less than this fraction is aligned (default: 0 = off)")
    p.add_argument("--min-identity", type=float, default=0.0,
                   help="drop reads whose best alignment has a lower identity, 1 - NM/length (default: 0 = off)")
    args = p.parse_args(argv)

    s = select_best(args.paf_in, args.paf_out, args.table, args.one_per_read, args.min_mapq,
                    args.min_aln_len, args.min_aln_frac, args.min_identity)
    dropped = (s["reads_low_mapq"] + s["reads_short"] + s["reads_low_fraction"]
               + s["reads_low_identity"])
    summary = (f"Alignments: {s['lines']} for {s['reads']} reads; "
               f"{s['reads_multi']} reads with several alignments ({s['reads_AS_tie']} tied on AS "
               f"between references); {s['best_secondary']} reads whose best-AS alignment is a secondary one. "
               f"Written: {s['lines_out']} alignments of {s['reads_out']} reads")
    if dropped:
        summary += (f"; dropped {dropped} reads ({s['reads_low_mapq']} low MAPQ, {s['reads_short']} short "
                    f"alignment, {s['reads_low_fraction']} low aligned fraction, "
                    f"{s['reads_low_identity']} low identity)")
    print(summary, file=sys.stderr if args.paf_out == "-" else sys.stdout)
    if s["bad_lines"]:
        print(f"WARNING: {s['bad_lines']} malformed PAF lines skipped", file=sys.stderr)
    if s["lines"] and s["with_AS"] == 0:
        sys.exit("select_best_alignments.py: no alignment has an AS tag; the PAF must be "
                 "produced with minimap2 -c (base-level alignment)")
    for tag in ("AS", "NM", "tp"):
        if s[f"missing_{tag}"]:
            print(f"WARNING: {s['missing_' + tag]} of {s['lines']} alignments lack the {tag} tag",
                  file=sys.stderr)
    if s["reads_no_AS"]:
        print(f"WARNING: {s['reads_no_AS']} reads have no AS on any alignment; "
              "their primary (or first) alignment was used", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
