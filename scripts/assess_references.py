#!/usr/bin/env python3
"""Per-reference evidence for the taxonomic assignment of the Stage 2 ECMSD PAF.

Input is the PAF written by select_best_alignments.py: for every read the alignments
to all references that tie for its best alignment score (AS), tagged with nc (number
of compatible references) and dm (AS margin to the best lower-scoring competitor).
A PAF without these tags is also accepted; they are then derived from the alignments
of each read, but a competitor that was not reported cannot be seen.

For every reference two kinds of evidence are calculated:
  compatible    all reads for which the reference has the best AS, including reads that
                tie with other references. Evidence that the reference is compatible with
                the data, e.g. because it is almost identical to a species that is present.
  discriminatory reads for which this reference is the only one with the best AS, by a
                margin of at least --min-margin AS points over every competitor that was
                reported, and whose alignment is of adequate quality (length, aligned
                fraction of the read and identity). Independent evidence for presence.
and for each of them the breadth (% of positions covered), the mean depth, the
evenness (normalised Shannon entropy of the read positions over 20 equal windows of
the reference: 1 = uniform, 0 = all reads in one window). For the discriminatory reads
the table also gives the pooled identity (1 - summed NM / summed alignment length), the
fraction of reads without a mismatch, and the HOPS edit distance decline (-delta%,
Huebler et al. 2019): over the read counts with 0, 1, 2 and 3 mismatches, the summed
decreases between neighbouring counts divided by all summed changes, 1 for a strictly
declining distribution as expected for reads of the reference species. Shared reads are left out of the identity: they come from
conserved regions and match distant references as well as the right one. The identity
is only reported, never used to decide; compare references within a sample, because
sequencing errors, damage and haplotype differences shift all references of a sample.

PCR duplicates: a read whose compatible alignments have the same references, start,
end and strand as an earlier read is counted once (disable with --no-dedup). This is
exact for merged reads; for unmerged mates of equal length unrelated reads can share
a position at high depth, so duplicate_fraction is an upper bound there.

Step A, credible references: a reference is credible when its discriminatory reads
reach --min-disc-reads, cover at least --min-disc-breadth percent of the reference and
are spread evenly enough over it (--min-disc-evenness), and reach at least --min-window-ratio
times the number of 1 kb windows that as many reads placed at random would hit. Under random
placement this ratio stays above 0.6 for 20 or more reads, so 0.5 only removes references whose
evidence piles up in a few regions, e.g. a conserved or contaminated stretch of the reference. No threshold relates to the
abundance of other references, so a rare species is not removed for being rare.

Step B, which uses all compatible alignments to calculate the coverage of the credible
references, and the assignment of reads, is done by LinkTaxonomy.py with --Evidence.

Several references per species: reference names are "<taxid>" or "<taxid>|<accession>".
With --group-by-taxid all references of one taxid (e.g. haplotypes of a species) form one
unit: a read is discriminatory when its best alignments all lie on references of that one
taxid, by a margin over the best reference of another taxid (dt tag of
select_best_alignments.py, else dm). Breadth and evenness of a unit pool the positions of
its references (one alignment per read and unit); the table then has one row per taxid,
whose length is that of its longest reference.

A unique best alignment is only provisional evidence: minimap2 reports a limited number
of secondary alignments chosen before base-level alignment, so a better or equal
reference may not have been reported for a read. Discriminatory reads for which the
number of reported alignments exceeds --max-alns (the minimap2 -N of Stage 2) are counted
in the column discriminatory_at_limit.
"""
import argparse
import gzip
import math
import sys

NBINS = 20
SUBWIN = 100   # bp; the densest-window share slides in steps of SUBWIN
WINDOW = 1000  # bp; window of discriminatory_windows and max_1kb_share
HEADER = ["reference", "taxon_name", "ref_length", "compatible_reads", "discriminatory_reads",
          "weak_unique_reads", "shared_reads", "discriminatory_at_limit", "compatible_breadth", "discriminatory_breadth",
          "compatible_depth", "discriminatory_depth", "compatible_evenness",
          "discriminatory_evenness", "breadth_ratio", "discriminatory_windows", "window_ratio", "max_1kb_share", "discriminatory_identity", "discriminatory_exact", "edit_distance_decline", "duplicate_fraction", "credible", "reason",
          "credible_overlap", "call_status", "call_group", "call_reason"]


def opener(path, encoding=None):
    if path == "-":
        return sys.stdin
    return gzip.open(path, "rt", encoding=encoding) if path.endswith(".gz") else open(path, "rt", encoding=encoding)


def load_names(path):
    """taxid -> scientific name from an NCBI names.dmp."""
    names = {}
    with opener(path, "latin-1") as fh:
        for line in fh:
            if "scientific name" in line:
                its = [x.replace("\t", "") for x in line.rstrip().split("|")]
                names[its[0]] = its[1]
    return names


def iter_reads(fh):
    """Yield the PAF fields of all lines of one read; they are adjacent in the file."""
    name, group = None, []
    for line in fh:
        f = line.rstrip("\r\n").split("\t")
        if len(f) < 12:
            continue  # stray log text or blank line
        if f[0] != name:
            if group:
                yield group
            name, group = f[0], []
        group.append(f)
    if group:
        yield group


def parse_alignment(f):
    """(ref, ref_len, start, end, qlen, qspan, alen, AS, NM, dm, na, strand, dt) of one PAF line.
    dt is None without the tag and -1 when no reference of another taxid was reported."""
    AS = NM = dm = na = dt = None
    for tag in f[12:]:
        key = tag[:5]
        if key == "AS:i:":
            AS = int(tag[5:])
        elif key == "NM:i:":
            NM = int(tag[5:])
        elif key == "dm:i:":
            dm = int(tag[5:])
        elif key == "na:i:":
            na = int(tag[5:])
        elif key == "dt:i:":
            dt = int(tag[5:])
    return (f[5], int(f[6]), int(f[7]), int(f[8]), int(f[1]), int(f[3]) - int(f[2]), int(f[10]), AS, NM, dm, na, f[4], dt)


def ref_taxid(ref):
    """Taxid of a reference named "<taxid>" or "<taxid>|<accession>"."""
    return ref.split("|", 1)[0]


def identity(a):
    """1 - NM / alignment length, or None without NM or alignment."""
    alen, nm = a[6], a[8]
    return None if nm is None or alen <= 0 else 1.0 - nm / alen


def good_quality(a, min_aln_len, min_aln_frac, min_identity):
    """Is the alignment long, complete and accurate enough to count as evidence?"""
    alen, qlen, qspan, ident = a[6], a[4], a[5], identity(a)
    if alen < min_aln_len:
        return False
    if qlen > 0 and qspan / qlen < min_aln_frac:
        return False
    if ident is not None and ident < min_identity:
        return False
    return True


class RefEvidence:
    __slots__ = ("length", "compat", "disc", "weak", "limit", "cdiff", "ddiff", "cbins", "dbins", "dnm", "dalen", "dn", "dexact", "dhist", "dups", "dsub", "co")

    def __init__(self, length):
        self.length = length
        self.compat = self.disc = self.weak = self.limit = 0
        self.cdiff, self.ddiff = {}, {}      # sparse difference arrays: position -> delta depth
        self.cbins, self.dbins = [0] * NBINS, [0] * NBINS
        self.dnm = self.dalen = self.dn = self.dexact = 0   # summed NM, alignment length, reads, reads with NM 0
        self.dhist = [0] * 4                                 # reads with NM 0, 1, 2, 3
        self.dups = 0                                        # compatible reads skipped as PCR duplicates
        self.dsub = {}                                       # SUBWIN bp bin of the read midpoint -> discriminatory reads
        self.co = {}                                         # other reference -> reads compatible with both

    def add_identity(self, a):
        """Count a discriminatory read in the identity (reads without NM are skipped)."""
        if a[8] is None or a[6] <= 0:
            return
        self.dnm += a[8]
        self.dalen += a[6]
        self.dn += 1
        self.dexact += a[8] == 0
        if a[8] < len(self.dhist):
            self.dhist[a[8]] += 1

    def add(self, a, discriminatory):
        start, end = a[2], a[3]
        # bin on the alignment's own reference: references of a taxid unit differ in length
        mid = min(NBINS - 1, int((start + end) / 2 * NBINS / a[1])) if a[1] > 0 else 0
        for diff, bins in ((self.cdiff, self.cbins),) + (((self.ddiff, self.dbins),) if discriminatory else ()):
            diff[start] = diff.get(start, 0) + 1
            diff[end] = diff.get(end, 0) - 1
            bins[mid] += 1
        if discriminatory:
            b = (start + end) // 2 // SUBWIN
            self.dsub[b] = self.dsub.get(b, 0) + 1


def collect(fh, min_margin=10, min_aln_len=20, min_aln_frac=0.5, min_identity=0.9, max_alns=20, dedup=True,
            group_by_taxid=False):
    """Read the PAF once; return ({ref: RefEvidence}, read statistics). With group_by_taxid the
    keys are taxids and all references of a taxid count as one (see module doc)."""
    unit = ref_taxid if group_by_taxid else (lambda r: r)
    refs = {}
    stats = dict(reads=0, unique=0, shared=0, discriminatory=0, weak=0, duplicates=0)
    seen = set()   # molecule keys of the reads counted so far
    for group in iter_reads(fh):
        alns = [parse_alignment(f) for f in group]
        scores = [a[7] for a in alns if a[7] is not None]
        top = max(scores) if scores else None
        compat = [a for a in alns if top is None or a[7] == top]
        if group_by_taxid:
            # one alignment per unit (the first, as select_best_alignments writes the best first),
            # renamed to the unit so that everything below counts units
            first = {}
            for a in compat:
                first.setdefault(unit(a[0]), (unit(a[0]),) + a[1:])
            compat = list(first.values())
            alns = compat + [(unit(a[0]),) + a[1:] for a in alns if a[7] != top]
        compat_refs = {a[0] for a in compat}
        if dedup:
            key = frozenset((a[0], a[2], a[3], a[11]) for a in compat)
            if key in seen:
                stats["duplicates"] += 1
                for ref in compat_refs:
                    ev = refs.get(ref)
                    if ev is not None:  # always set: the first copy created it
                        ev.dups += 1
                continue
            seen.add(key)
        lower = [a[7] for a in alns if top is not None and a[7] is not None and a[7] < top
                 and a[0] not in compat_refs]
        dm = next((a[9] for a in compat if a[9] is not None), None)
        if group_by_taxid:
            dt = next((a[12] for a in compat if a[12] is not None), None)
            if dt is not None:
                dm = None if dt < 0 else dt
        if dm is None and lower:
            dm = top - max(lower)
        unique = len(compat_refs) == 1
        discriminatory = (unique and (dm is None or dm >= min_margin)
                          and good_quality(compat[0], min_aln_len, min_aln_frac, min_identity))
        at_limit = any(a[10] is not None and a[10] > max_alns for a in compat)
        stats["reads"] += 1
        stats["unique" if unique else "shared"] += 1
        stats["discriminatory" if discriminatory else "weak"] += 1 if unique else 0
        counted = set()
        for a in compat:
            ev = refs.get(a[0])
            if ev is None:
                ev = refs[a[0]] = RefEvidence(a[1])
            ev.length = max(ev.length, a[1])
            ev.add(a, discriminatory)
            if a[0] not in counted:
                counted.add(a[0])
                ev.compat += 1
                for other in compat_refs:
                    if other != a[0]:
                        ev.co[other] = ev.co.get(other, 0) + 1
                if unique:
                    if discriminatory:
                        ev.disc += 1
                        ev.add_identity(a)
                        ev.limit += at_limit
                    else:
                        ev.weak += 1
    return refs, stats


def sweep(diff):
    """(covered positions, summed depth) of a sparse difference array."""
    covered = total = depth = prev = 0
    for pos in sorted(diff):
        span = pos - prev
        if depth > 0:
            covered += span
            total += depth * span
        depth += diff[pos]
        prev = pos
    return covered, total


def evenness(bins):
    """Normalised Shannon entropy of read counts per window (1 uniform, 0 concentrated)."""
    n = sum(bins)
    if n <= 1:
        return 0.0
    return max(0.0, -sum(c / n * math.log(c / n) for c in bins if c) / math.log(len(bins)))


def breadth_ratio(covered_frac, depth):
    """Observed breadth / breadth expected for reads placed at random, 1 - exp(-depth) (Lander-Waterman)."""
    return covered_frac / (1.0 - math.exp(-depth)) if depth > 0 else None


def spatial(sub, length, n):
    """(fraction of WINDOW-bp windows with a discriminatory read, share of reads in the densest
    circular WINDOW-bp window) from read midpoints binned per SUBWIN bp."""
    if n == 0 or length <= 0:
        return None, None
    nsub = -(-length // SUBWIN)
    k = min(WINDOW // SUBWIN, nsub)
    counts = [sub.get(i, 0) for i in range(nsub)]
    run = sum(counts[:k])
    best = run
    for i in range(1, nsub):  # circular genome: the window wraps around the end
        run += counts[(i + k - 1) % nsub] - counts[i - 1]
        best = max(best, run)
    nwin = -(-length // WINDOW)
    hit = len({i * SUBWIN // WINDOW for i, c in enumerate(counts) if c})
    return hit / nwin, best / n


def edit_distance_decline(hist):
    """HOPS -delta%: summed decreases / summed changes between neighbouring counts; None without change."""
    d = [b - a for a, b in zip(hist, hist[1:])]
    total = sum(abs(x) for x in d)
    return sum(-x for x in d if x < 0) / total if total else None


def window_ratio(hit_frac, length, n):
    """Fraction of WINDOW-bp windows hit / fraction expected for n reads placed at random, 1 - (1 - 1/W)^n."""
    if hit_frac is None:
        return None
    w = -(-length // WINDOW)
    return hit_frac / (1.0 - (1.0 - 1.0 / w) ** n)


def judge(ev, min_disc_reads, min_disc_breadth, min_disc_evenness, min_window_ratio=0.0):
    """Step A. Return (row metrics, credible, reasons)."""
    L = ev.length
    cc, ct = sweep(ev.cdiff)
    dc, dt = sweep(ev.ddiff)
    m = dict(compat_breadth=100.0 * cc / L if L else 0.0, disc_breadth=100.0 * dc / L if L else 0.0,
             compat_depth=ct / L if L else 0.0, disc_depth=dt / L if L else 0.0,
             compat_even=evenness(ev.cbins), disc_even=evenness(ev.dbins),
             disc_ident=1.0 - ev.dnm / ev.dalen if ev.dalen else None,
             disc_exact=ev.dexact / ev.dn if ev.dn else None,
             decline=edit_distance_decline(ev.dhist),
             dup_frac=ev.dups / (ev.compat + ev.dups) if ev.compat + ev.dups else None)
    m["breadth_ratio"] = breadth_ratio(dc / L, dt / L) if L else None
    m["disc_windows"], m["max_1kb"] = spatial(ev.dsub, L, ev.disc)
    m["window_ratio"] = window_ratio(m["disc_windows"], L, ev.disc)
    reasons = []
    if ev.disc < min_disc_reads:
        reasons.append("few_discriminatory_reads")
    if m["disc_breadth"] < min_disc_breadth:
        reasons.append("low_discriminatory_breadth")
    if m["disc_even"] < min_disc_evenness:
        reasons.append("uneven_discriminatory_coverage")
    if (m["window_ratio"] or 0.0) < min_window_ratio:
        reasons.append("clustered_discriminatory_reads")
    return m, not reasons, reasons


def fmt(x, digits):
    return "NA" if x is None else "%.*f" % (digits, x)


def jaccard(refs, a, b):
    """Share of the reads compatible with a or b that are compatible with both."""
    both = refs[a].co.get(b, 0)
    union = refs[a].compat + refs[b].compat - both
    return both / union if union else 0.0


def call(refs, metrics, credible, min_disc_reads, min_decline=0.9, min_breadth_ratio=0.5, max_overlap=0.5):
    """Label the credible references: {ref: (overlap, status, group, reason)}. With --group-by-taxid
    the references are taxid units, so overlap, decline and breadth ratio are per unit.

    poor_fit      edit distance decline < min_decline, or the discriminatory reads cover less than
                  min_breadth_ratio of the breadth expected from their depth (reads from one region)
    ambiguous     shares more than max_overlap (Jaccard) of its compatible reads with another credible
                  reference; such references are linked into one group (single linkage)
    low_evidence  credible only through discriminatory reads at the -N limit
    confident     none of the above. poor_fit outranks ambiguous, which outranks low_evidence.
    """
    cred = sorted(credible, key=lambda r: (-refs[r].disc, r))
    parent = {r: r for r in cred}

    def root(r):
        while parent[r] != r:
            r = parent[r]
        return r

    overlap = {r: [] for r in cred}
    for i, a in enumerate(cred):
        for b in cred[i + 1:]:
            j = jaccard(refs, a, b)
            if j >= 0.05:
                overlap[a].append((b, j))
                overlap[b].append((a, j))
            if j > max_overlap:
                parent[root(b)] = root(a)
    members = {}
    for r in cred:
        members.setdefault(root(r), []).append(r)
    gid = {}
    for n, (rt, ms) in enumerate((x for x in members.items() if len(x[1]) > 1), 1):
        for r in ms:
            gid[r] = "G%d" % n
    out = {}
    for r in cred:
        ev, m = refs[r], metrics[r]
        ratio = m["breadth_ratio"]
        poor, why = [], []
        if m["decline"] is not None and m["decline"] < min_decline:
            poor.append("edit_distance_decline<%g" % min_decline)
        if ratio is not None and ratio < min_breadth_ratio:
            poor.append("breadth_ratio<%g" % min_breadth_ratio)
        low = ev.disc - ev.limit < min_disc_reads
        if r in gid:
            why.append("shares>%g_with:%s" % (max_overlap, ",".join(x for x in members[root(r)] if x != r)))
        if low:
            why.append("credible_only_with_reads_at_N_limit")
        status = "poor_fit" if poor else "ambiguous" if r in gid else "low_evidence" if low else "confident"
        out[r] = (",".join("%s:%.2f" % (b, j) for b, j in sorted(overlap[r], key=lambda x: -x[1])) or "-",
                  status, gid.get(r, "-"), ";".join(poor + why) or "-")
    return out


def write_table(refs, out, names, min_disc_reads, min_disc_breadth, min_disc_evenness, min_window_ratio=0.0,
                **call_args):
    n_credible = 0
    rows, metrics, credible = [], {}, []
    for ref, ev in refs.items():
        m, ok, reasons = judge(ev, min_disc_reads, min_disc_breadth, min_disc_evenness, min_window_ratio)
        n_credible += ok
        metrics[ref] = m
        if ok:
            credible.append(ref)
        rows.append((ev.disc, ev.compat, ref, ev, m, ok, reasons))
    calls = call(refs, metrics, credible, min_disc_reads, **call_args)
    rows.sort(key=lambda r: (-r[0], -r[1], r[2]))
    out.write("\t".join(HEADER) + "\n")
    for _, _, ref, ev, m, ok, reasons in rows:
        ovl, status, group, why = calls.get(ref, ("-", "not_credible", "-", "-"))
        out.write("\t".join([
            ref, names.get(ref_taxid(ref), "NA").replace(" ", "_"), str(ev.length), str(ev.compat), str(ev.disc),
            str(ev.weak), str(ev.compat - ev.disc - ev.weak), str(ev.limit),
            "%.2f" % m["compat_breadth"], "%.2f" % m["disc_breadth"],
            "%.2f" % m["compat_depth"], "%.2f" % m["disc_depth"],
            "%.3f" % m["compat_even"], "%.3f" % m["disc_even"],
            fmt(m["breadth_ratio"], 3), fmt(m["disc_windows"], 3), fmt(m["window_ratio"], 3), fmt(m["max_1kb"], 3),
            fmt(m["disc_ident"], 4), fmt(m["disc_exact"], 3),
            fmt(m["decline"], 3), fmt(m["dup_frac"], 3),
            "yes" if ok else "no", "pass" if ok else ";".join(reasons),
            ovl, status, group, why]) + "\n")
    return n_credible


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--PAF", required=True, help="Stage 2 PAF from select_best_alignments.py (plain or .gz)")
    p.add_argument("--output", required=True, help="evidence table (TSV)")
    p.add_argument("--Names", help="NCBI names.dmp, only used to add the scientific name to the table")
    p.add_argument("--min-margin", type=int, default=10,
                   help="minimum AS margin over the best competing reference for a discriminatory read (default: 10)")
    p.add_argument("--min-aln-len", type=int, default=20,
                   help="minimum alignment length of a discriminatory read (default: 20)")
    p.add_argument("--min-aln-frac", type=float, default=0.5,
                   help="minimum aligned fraction of a discriminatory read (default: 0.5)")
    p.add_argument("--min-identity", type=float, default=0.9,
                   help="minimum identity, 1 - NM/alignment length, of a discriminatory read (default: 0.9)")
    p.add_argument("--max-alns", type=int, default=20,
                   help="the minimap2 -N used in Stage 2; discriminatory reads with more reported alignments are counted as at the limit (default: 20)")
    p.add_argument("--no-dedup", action="store_true",
                   help="count PCR duplicates (reads with the same references, start, end and strand) as separate reads")
    p.add_argument("--min-disc-reads", type=int, default=20,
                   help="credible: minimum number of discriminatory reads (default: 20)")
    p.add_argument("--min-disc-breadth", type=float, default=5.0,
                   help="credible: minimum %% of the reference covered by discriminatory reads (default: 5)")
    p.add_argument("--min-disc-evenness", type=float, default=0.5,
                   help="credible: minimum evenness (0-1) of the discriminatory reads over the reference (default: 0.5)")
    p.add_argument("--group-by-taxid", action="store_true",
                   help="count all references of one taxid (names <taxid>|<accession>) as one unit")
    p.add_argument("--min-window-ratio", type=float, default=0.5,
                   help="credible: minimum ratio of the fraction of 1 kb windows hit by discriminatory reads to the fraction expected for as many reads placed at random; 0 disables (default: 0.5)")
    p.add_argument("--min-decline", type=float, default=0.9,
                   help="call: poor_fit below this edit distance decline of the discriminatory reads (default: 0.9)")
    p.add_argument("--min-breadth-ratio", type=float, default=0.5,
                   help="call: poor_fit when the discriminatory breadth is below this fraction of the breadth expected from their depth (default: 0.5)")
    p.add_argument("--max-overlap", type=float, default=0.5,
                   help="call: credible references sharing more than this fraction (Jaccard) of their compatible reads form an ambiguous group (default: 0.5)")
    args = p.parse_args(argv)

    with opener(args.PAF) as fh:
        refs, stats = collect(fh, args.min_margin, args.min_aln_len, args.min_aln_frac, args.min_identity,
                              args.max_alns, dedup=not args.no_dedup, group_by_taxid=args.group_by_taxid)
    names = load_names(args.Names) if args.Names else {}
    with open(args.output, "w") as out:
        n_credible = write_table(refs, out, names, args.min_disc_reads, args.min_disc_breadth,
                                 args.min_disc_evenness, args.min_window_ratio, min_decline=args.min_decline,
                                 min_breadth_ratio=args.min_breadth_ratio, max_overlap=args.max_overlap)
    print(f"Reference evidence: {stats['reads']} reads ({stats['duplicates']} PCR duplicates skipped), {stats['shared']} shared between references, "
          f"{stats['discriminatory']} discriminatory, {stats['weak']} with a unique but weak best alignment; "
          f"{n_credible} of {len(refs)} references credible. Written to: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
