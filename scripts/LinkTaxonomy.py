from collections import defaultdict as d
import sys
import functools
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
                  help="Taxonomic rank to use as the label in the ref_summary (e.g. species, genus, family). Default: species",
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
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
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
    The taxID is extracted from the reference name as reference.split('|')[1].
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
                taxid = ref_name.split("|")[1]
                passing.add(taxid)
    return passing


def main():

    #check if all required options are provided
    if not options.Tax or not options.Names or not options.PAF or not options.OUT:
        parser.print_help()
        sys.exit(1) 

    print("Loading taxonomy data from nodes.dmp and names.dmp ...")
    with load_data(options.Tax) as node:

        for line in node:
            # Split the line by tab and remove leading/trailing whitespace
            its = [x.replace("\t", "") for x in line.rstrip().split("|")]
            parents[its[0]] = its[1]
            rank_dict[its[0]] = its[2]

    print("Loading names data from names.dmp ...")
    with load_data(options.Names) as name:
        for line in name:
            its = [x.replace("\t", "") for x in line.rstrip().split("|")]
            if "scientific name" in line:
                names[its[0]] = its[1]

passing_taxids = load_coverage(options.Coverage, float(options.CovThreshold))
print(
    f"References passing coverage threshold (>= {options.CovThreshold}% of reference covered): {len(passing_taxids)}")


def get_alignment_score(fields):
    """Return the best available alignment score from PAF optional fields.
    Primary alignments carry ms:i (mate score / alignment score).
    Secondary alignments carry s1:i (chaining score) instead.
    Falls back to 0 if neither tag is present."""
    for f in fields[12:]:
        if f.startswith("ms:i:"):
            return int(f[5:])
        if f.startswith("s1:i:"):
            return int(f[5:])
    return 0


def load_best_hits(paf_path, min_mapq):
    """Two-pass winner-takes-all disambiguation.

    For every read, collect all alignments (primary + secondary) that pass
    basic quality filters, then keep only the alignment(s) that map to the
    single best-scoring reference.  This prevents reads from a species A
    contributing coverage to a closely-related species B via secondary hits.

    Returns a dict: read_name -> list of raw PAF line strings (winners only).
    """
    # First pass: determine the best score per read
    best_score = {}          # read_name -> highest alignment score seen
    best_ref = {}            # read_name -> ref_name of that best alignment

    with load_data(paf_path) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            read_name = fields[0]
            mapq = int(fields[11])
            is_secondary = any(f == "tp:A:S" for f in fields[12:])
            if not is_secondary and mapq < min_mapq:
                continue
            score = get_alignment_score(fields)
            if score > best_score.get(read_name, -1):
                best_score[read_name] = score
                best_ref[read_name] = fields[5]

    # Second pass: keep only alignments to the winner reference
    winners = d(list)
    with load_data(paf_path) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            read_name = fields[0]
            ref_name = fields[5]
            mapq = int(fields[11])
            is_secondary = any(f == "tp:A:S" for f in fields[12:])
            if not is_secondary and mapq < min_mapq:
                continue
            if best_ref.get(read_name) == ref_name:
                winners[read_name].append(fields)

    n_disambiguated = sum(
        1 for rd, hits in winners.items() if len(hits) > 1
    )
    print(f"Winner-takes-all disambiguation: {len(winners)} reads retained, "
          f"{n_disambiguated} reads had multiple refs collapsed to best hit.")
    return winners


print("Loading PAF and disambiguating secondary alignments (winner-takes-all)...")
best_hits = load_best_hits(options.PAF, options.MapQuality)

# Count reads per reference (ref_name, taxId) for ranking
ref_read_counts = d(int)
ref_to_taxid = {}
ref_to_name = {}
ref_to_species = {}

with open(options.OUT+".txt", 'wt') as export:
    # write header to output file
    export.write(
        "SeqID\tTaxID\tLength\tMappingQuality\tdomain\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
    for seqId, hit_list in best_hits.items():
        # Each read now maps to exactly one reference; take the first (and only) entry
        lines = hit_list[0]
        ref_name = lines[5]
        taxId = ref_name.split("|")[1]

        if taxId not in passing_taxids:
            continue

            if tax_id not in names:
                continue

            node_path, name_path = taxon_trace(tax_id)
            # loop through the following ranks "no rank|cellular root|domain|kingdom|phylum|class|order|family|genus|species" and make new list from name path. replace with NA if not available
            rank_name = dict(zip(node_path.split("|"), name_path.split("|")))
            name_path = []
            for rank in ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]:
                if rank in rank_name:
                    name_path.append(rank_name[rank])
                else:
                    name_path.append("NA")

            Tax = "\t".join(["_".join(x.split()) for x in name_path])
            export.write(sequence_id+"\t"+tax_id+"\t" +
                        lines[1]+"\t"+lines[11] + "\t"+Tax+"\n")

        # Track read counts per ref_name for ranking
        ref_read_counts[ref_name] += 1
        ref_to_taxid[ref_name] = taxId
        # Build a clean taxon label based on the chosen taxonomic hierarchy
        rank_index = {"domain": 0, "kingdom": 1, "phylum": 2, "class": 3,
                      "order": 4, "family": 5, "genus": 6, "species": 7}
        chosen_rank = options.TaxonomicHierarchy.lower()
        idx = rank_index.get(chosen_rank, 7)  # default to species
        label = name_path[idx] if name_path[idx] != "NA" else "Unknown"
        # For ranks above species the NCBI name is already a single word;
        # for species it is a full binomial ("Genus epithet") — use it as-is.
        ref_to_name[ref_name] = "_".join(label.split())
        # Always store the species binomial independently (used in plot titles)
        species_label = name_path[7] if name_path[7] != "NA" else "Unknown sp"
        ref_to_species[ref_name] = "_".join(species_label.split())

# Write ranked reference summary for pafr plotting
with open(options.OUT+".ref_summary.txt", 'wt') as ref_out:
    ref_out.write("ref_name\ttaxid\ttaxon_name\tspecies_name\tread_count\n")
    for ref_name, count in sorted(ref_read_counts.items(), key=lambda x: -x[1]):
        ref_out.write(
            f"{ref_name}\t{ref_to_taxid[ref_name]}\t{ref_to_name[ref_name]}\t{ref_to_species[ref_name]}\t{count}\n")

print(f"Reference summary written to: {options.OUT}.ref_summary.txt")


    sys.exit()

if __name__ == "__main__":
    main()
