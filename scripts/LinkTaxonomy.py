import functools
from collections import defaultdict as d
import sys
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

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
                  help="Taxonomic rank to use as the label in the ref_summary (e.g. subspecies, species, genus, family). Default: species. If subspecies is chosen but a reference has no subspecies rank, the species name is used instead",
                  default="species")
parser.add_option("--MapQuality", dest="MapQuality",
                  help="Minimum mapping quality for primary alignments (default: 20)",
                  default=20, type="int")

(options, args) = parser.parse_args()
parser.add_option_group(group)

parents = {}
rank_dict = {}
names = {}

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

# Load taxonomy first (populates parents, rank_dict, names module-level dicts)
main()

passing_taxids = load_coverage(options.Coverage, float(options.CovThreshold))
print(
    f"References passing coverage threshold (>= {options.CovThreshold}% of reference covered): {len(passing_taxids)}")


def get_alignment_score(fields):
    """Return the best available alignment score from PAF optional fields."""
    for f in fields[12:]:
        if f.startswith("ms:i:"):
            return int(f[5:])
    return 0


def load_best_hits(paf_path, min_mapq):
    """Single-pass: keep every primary alignment that passes MQ.
    Returns a dict: read_name -> list of PAF fields for that read.
    """
    winners = d(list)
    with load_data(paf_path) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            read_name = fields[0]
            mapq = int(fields[11])
            if mapq < min_mapq:
                continue
            winners[read_name].append(fields)

    print(f"PAF loaded: {len(winners)} reads passed MQ filter.")
    return winners


print("Loading PAF alignments...")
best_hits = load_best_hits(options.PAF, options.MapQuality)

# Count reads per reference (ref_name, taxId) for ranking
ref_read_counts = d(int)
ref_to_taxid = {}
ref_to_name = {}
ref_to_species = {}
ref_to_subspecies = {}

with open(options.OUT+".txt", 'wt') as export:
    export.write(
        "SeqID\tTaxID\tLength\tMappingQuality\tdomain\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tsubspecies\n")
    for seqId, hit_list in best_hits.items():
        lines = hit_list[0]
        ref_name = lines[5]
        taxId = ref_name

        if taxId not in passing_taxids:
            continue

        if taxId not in names:
            continue

        node_path, name_path = taxon_trace(taxId)
        rank_name = dict(zip(node_path.split("|"), name_path.split("|")))
        name_path = []
        for rank in ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies"]:
            if rank in rank_name:
                name_path.append(rank_name[rank])
            else:
                name_path.append("NA")

        Tax = "\t".join(["_".join(x.split()) for x in name_path])
        export.write(seqId + "\t" + taxId + "\t" +
                     lines[1] + "\t" + lines[11] + "\t" + Tax + "\n")

        ref_read_counts[ref_name] += 1
        ref_to_taxid[ref_name] = taxId
        rank_index = {"domain": 0, "kingdom": 1, "phylum": 2, "class": 3,
                      "order": 4, "family": 5, "genus": 6, "species": 7,
                      "subspecies": 8}
        chosen_rank = options.TaxonomicHierarchy.lower()
        idx = rank_index.get(chosen_rank, 7)
        label = name_path[idx]
        # most references are only resolved to species level; fall back rather
        # than labelling them all "Unknown"
        if label == "NA" and chosen_rank == "subspecies":
            label = name_path[7]
        if label == "NA":
            label = "Unknown"
        ref_to_name[ref_name] = "_".join(label.split())
        species_label = name_path[7] if name_path[7] != "NA" else "Unknown sp"
        ref_to_species[ref_name] = "_".join(species_label.split())
        subspecies_label = name_path[8] if name_path[8] != "NA" else "NA"
        ref_to_subspecies[ref_name] = "_".join(subspecies_label.split())

# Write ranked reference summary for plotting
with open(options.OUT+".ref_summary.txt", 'wt') as ref_out:
    ref_out.write(
        "ref_name\ttaxid\ttaxon_name\tspecies_name\tsubspecies_name\tread_count\n")
    for ref_name, count in sorted(ref_read_counts.items(), key=lambda x: -x[1]):
        ref_out.write(
            f"{ref_name}\t{ref_to_taxid[ref_name]}\t{ref_to_name[ref_name]}\t{ref_to_species[ref_name]}\t{ref_to_subspecies[ref_name]}\t{count}\n")

print(f"Reference summary written to: {options.OUT}.ref_summary.txt")
