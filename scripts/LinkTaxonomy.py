import functools
from collections import defaultdict as d
import sys
from optparse import OptionParser

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)

#########################################################   CODE   #########################################################################

parser.add_option("--Nodes", dest="Tax", help="NCBI node dmp file")
parser.add_option("--Names", dest="Names", help="NCBI names dmp file")
parser.add_option("--Coverage", dest="Coverage",
                  help="Coverage file produced by the coverage calculation step (TSV with columns: reference, ref_length, mean_coverage, std_coverage, pct_covered)")
parser.add_option("--CovThreshold", dest="CovThreshold",
                  help="Minimum percentage of reference covered by reads (0-100) to retain a reference (default: 50)", default=50)
parser.add_option("--PAF", dest="PAF",
                  help="PAF output file with SeqID in column1 and TaxID in column2")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--TaxonomicHierarchy", dest="TaxonomicHierarchy",
                  help="Taxonomic rank to use as the label in the ref_summary (e.g. subspecies, species, genus, family). Default: species. Where a reference is not resolved to the requested rank, the label walks up to the next coarser rank that is populated",
                  default="species")
parser.add_option("--MapQuality", dest="MapQuality",
                  help="Minimum mapping quality for primary alignments (default: 20)",
                  default=20, type="int")

parser.add_option("--Evidence", dest="Evidence",
                  help="Reference evidence table written by assess_references.py. When given, only references marked credible are reported, and each read is assigned to the credible reference it is compatible with (any equally best-scoring alignment in the PAF); a read compatible with several credible references is split between them in the proportion of their discriminatory reads (column Weight)")

(options, args) = parser.parse_args()

parents = {}
rank_dict = {}
names = {}

RANKS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies"]
RANK_INDEX = {rank: i for i, rank in enumerate(RANKS)}


def load_data(x):
    ''' import data either from a gzipped or uncompressed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


@functools.lru_cache(maxsize=None)
def taxon_trace(node):
    """Trace the taxonomic path from a node to the root."""
    rank = []
    name_path = []
    while True:
        rank.append(rank_dict[node])
        name_path.append(names.get(node, ""))

        if node == '1':
            break

        if node in parents:
            node = parents[node]
        else:
            sys.exit(f"{node}\tSomething may be wrong!")

    return "|".join(reversed(rank)), "|".join(reversed(name_path))


def load_coverage(coverage_file, threshold):
    """
    Loads the coverage file and returns a set of taxIDs whose percentage of
    covered reference positions meets or exceeds the given threshold.

    The coverage file is expected to have columns:
        reference  ref_length  mean_coverage  std_coverage  pct_covered
    The reference name is the taxID itself.
    Threshold is a percentage value (0-100).
    """
    passing = set()
    with open(coverage_file, "r") as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            ref_name = parts[0]
            try:
                pct_covered = float(parts[4])
            except ValueError:
                continue
            if pct_covered >= threshold:
                passing.add(ref_name)
    return passing


def fallback_label(name_path, chosen_rank):
    """
    Pick the label for `chosen_rank` out of a 9-element name_path (domain..subspecies,
    "NA" where unresolved). Most references are only resolved to species level, so
    where the requested rank is "NA" walk up to the next coarser rank that is
    populated rather than labelling them all "Unknown".
    """
    idx = RANK_INDEX[chosen_rank]
    label = name_path[idx]
    while label == "NA" and idx > 0:
        idx -= 1
        label = name_path[idx]
    if label == "NA":
        label = "Unknown"
    return label


def main():

    if not options.Tax or not options.Names or not options.PAF or not options.OUT or not options.Coverage:
        parser.print_help()
        sys.exit(1)

    print("Loading taxonomy data from nodes.dmp and names.dmp ...")
    with load_data(options.Tax) as node:
        for line in node:
            its = [x.replace("\t", "") for x in line.rstrip().split("|")]
            parents[its[0]] = its[1]
            rank_dict[its[0]] = its[2]

    print("Loading names data from names.dmp ...")
    with load_data(options.Names) as name:
        for line in name:
            its = [x.replace("\t", "") for x in line.rstrip().split("|")]
            if "scientific name" in line:
                names[its[0]] = its[1]


def load_best_hits(paf_path, min_mapq):
    """Single-pass: keep the first alignment per read that passes MQ.
    Returns a dict: read_name -> PAF fields of the first passing alignment.
    """
    winners = {}
    with load_data(paf_path) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            read_name = fields[0]
            mapq = int(fields[11])
            if mapq < min_mapq:
                continue
            winners.setdefault(read_name, fields)

    print(f"PAF loaded: {len(winners)} reads passed MQ filter.")
    return winners


def load_evidence(path):
    """References marked credible in the table written by assess_references.py.
    Returns a dict: reference -> number of discriminatory reads."""
    credible = {}
    with open(path, "r") as fh:
        header = next(fh).rstrip("\n").split("\t")
        try:
            i_ref, i_ok = header.index("reference"), header.index("credible")
            i_disc = header.index("discriminatory_reads")
        except ValueError:
            sys.exit(f"{path} is not a reference evidence table written by assess_references.py")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > max(i_ref, i_ok, i_disc) and parts[i_ok] == "yes":
                credible[parts[i_ref]] = int(parts[i_disc])
    return credible


def split_weights(members, credible):
    """Share of one read per reference: 1 for a single reference; for several credible
    references, their discriminatory reads / the sum over all of them (equal shares if
    the sum is 0, or without an evidence table). Returns a dict: reference -> weight."""
    if len(members) == 1:
        return {members[0]: 1.0}
    disc = [credible[r] if credible else 0 for r in members]
    total = sum(disc)
    if total == 0:
        return {r: 1.0 / len(members) for r in members}
    return {r: n / total for r, n in zip(members, disc)}


def load_compatible_hits(paf_path, min_mapq):
    """Per read: the PAF fields of its first alignment and the distinct references of all
    its compatible (equally best-scoring) alignments, as written by select_best_alignments.py.
    A read is kept when the best MAPQ of its alignments reaches min_mapq.
    Returns a dict: read_name -> (fields, [reference, ...]).
    """
    reads = {}
    with load_data(paf_path) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            mapq = int(fields[11])
            entry = reads.get(fields[0])
            if entry is None:
                reads[fields[0]] = [fields, [fields[5]], mapq]
            else:
                if fields[5] not in entry[1]:
                    entry[1].append(fields[5])
                entry[2] = max(entry[2], mapq)
    return {read: (e[0], e[1]) for read, e in reads.items() if e[2] >= min_mapq}


def run():
    if options.TaxonomicHierarchy not in RANK_INDEX:
        sys.exit(f"Unknown rank '{options.TaxonomicHierarchy}'. Choose from: {', '.join(RANKS)}")
    # Load taxonomy first (populates parents, rank_dict, names module-level dicts)
    main()

    passing_taxids = load_coverage(options.Coverage, float(options.CovThreshold))
    print(
        f"References passing coverage threshold (>= {options.CovThreshold}% of reference covered): {len(passing_taxids)}")

    credible = None
    if options.Evidence:
        credible = load_evidence(options.Evidence)
        passing_taxids &= set(credible)
        print(f"References with credible discriminatory evidence: {len(credible)}; "
              f"passing the coverage threshold as well: {len(passing_taxids)}")

    print("Loading PAF alignments...")
    if credible is None:
        best_hits = load_best_hits(options.PAF, options.MapQuality)
        hits = ((read, fields, [fields[5]]) for read, fields in best_hits.items())
    else:
        compatible = load_compatible_hits(options.PAF, options.MapQuality)
        print(f"PAF loaded: {len(compatible)} reads passed MQ filter.")
        hits = ((read, fields, refs) for read, (fields, refs) in compatible.items())

    # Count reads per reference (ref_name, taxId) for ranking
    ref_read_counts = d(int)
    ref_to_taxid = {}
    ref_to_name = {}
    ref_to_species = {}
    ref_to_subspecies = {}
    path_cache = {}
    n_assigned = n_shared = n_unassigned = 0

    def name_path_of(taxId):
        if taxId not in path_cache:
            node_path, name_path_str = taxon_trace(taxId)
            rank_name = dict(zip(node_path.split("|"), name_path_str.split("|")))
            path_cache[taxId] = [rank_name.get(rank, "NA") for rank in RANKS]
        return path_cache[taxId]

    with open(options.OUT+".txt", 'wt') as export:
        export.write(
            "SeqID\tTaxID\tLength\tMappingQuality\tdomain\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tsubspecies\tWeight\n")
        for seqId, lines, refs in hits:
            members = [r for r in refs if r in passing_taxids and r in names]
            if not members:
                n_unassigned += 1
                continue
            # A read compatible with several credible references is split between them in the
            # proportion of their discriminatory reads; only credible references are ever
            # reported, never a higher taxonomic level
            if len(members) > 1:
                n_shared += 1
            n_assigned += 1
            for taxId, weight in split_weights(members, credible).items():
                name_path = name_path_of(taxId)
                Tax = "\t".join(["_".join(x.split()) for x in name_path])
                export.write(seqId + "\t" + taxId + "\t" + lines[1] + "\t" + lines[11] + "\t" +
                             Tax + "\t" + f"{weight:.6g}" + "\n")

            for ref_name in members:
                ref_read_counts[ref_name] += 1
                if ref_name not in ref_to_taxid:
                    ref_path = name_path_of(ref_name)
                    ref_to_taxid[ref_name] = ref_name
                    label = fallback_label(ref_path, options.TaxonomicHierarchy.lower())
                    ref_to_name[ref_name] = "_".join(label.split())
                    species_label = ref_path[7] if ref_path[7] != "NA" else "Unknown sp"
                    ref_to_species[ref_name] = "_".join(species_label.split())
                    subspecies_label = ref_path[8] if ref_path[8] != "NA" else "NA"
                    ref_to_subspecies[ref_name] = "_".join(subspecies_label.split())

    if credible is not None:
        print(f"Reads assigned: {n_assigned}, of which {n_shared} are compatible with several credible "
              f"references and were split between them in proportion to their discriminatory reads; "
              f"reads without a credible reference: {n_unassigned}")

    # Write ranked reference summary for plotting
    with open(options.OUT+".ref_summary.txt", 'wt') as ref_out:
        ref_out.write(
            "ref_name\ttaxid\ttaxon_name\tspecies_name\tsubspecies_name\tread_count\n")
        for ref_name, count in sorted(ref_read_counts.items(), key=lambda x: -x[1]):
            ref_out.write(
                f"{ref_name}\t{ref_to_taxid[ref_name]}\t{ref_to_name[ref_name]}\t{ref_to_species[ref_name]}\t{ref_to_subspecies[ref_name]}\t{count}\n")

    print(f"Reference summary written to: {options.OUT}.ref_summary.txt")


if __name__ == "__main__":
    run()
