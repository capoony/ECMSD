import numpy as np
from collections import defaultdict as d
import math
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
parser.add_option("--Bins", dest="bins",
                  help="bin size for RMUS calculation (default: 1000)", default=1000)
parser.add_option("--RMUS", dest="RMUS",
                  help="minimum Read Mapping Uniformity Score (default: 0.5)", default=0.5)
parser.add_option("--PAF", dest="PAF",
                  help="PAF output file with SeqID in column1 and TaxID in column2")
parser.add_option("--output", dest="OUT", help="Output file")

(options, args) = parser.parse_args()
parser.add_option_group(group)


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


def calculate_rmus_from_paf(paf_file, bin_size, OUT):
    """
    Calculates the Read Mapping Uniformity Score (RMUS) from a PAF file.

    Parameters:
        paf_file (str): Path to the PAF file.
        ref_length (int): Length of the reference sequence.
        bin_size (int): Size of each bin (default: 1000 bp).
        ref_name (str): (Optional) Reference name to filter for a specific target.

    Returns:
        float: RMUS value between 0 and 1.
    """

    refDict = d(int)
    BinDict = d(lambda: d(int))

    with load_data(paf_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue  # Not a valid PAF line
            target = fields[5].split("|")[1]  # Extract the target sequence ID
            target_start = int(fields[7])
            target_end = int(fields[8])
            ref_length = int(fields[6])
            refDict[target] = int(fields[6])
            # make partially over
            num_bins = math.ceil(ref_length / bin_size)
            start_bin = target_start // bin_size
            end_bin = min(num_bins, (target_end // bin_size) + 1)
            for b in range(start_bin, end_bin):
                BinDict[target][b] += 1

    RMUS = d(float)
    OUT.write("Name\tTaxID\tTotalReads\tRMUS\tBins\n")
    printlist = d(list)
    for k in refDict.keys():
        ref_length = refDict[k]
        total_reads = sum(BinDict[k].values())

        # Convert to probabilities
        probs = [x/total_reads for x in BinDict[k].values()]

        # Calculate entropy
        entropy = sum([x * math.log2(x) for x in probs if x > 0])

        # Normalize by maximum possible entropy
        max_entropy = math.log2(len(probs))
        rmus = entropy / max_entropy if max_entropy > 0 else -0.0
        RMUS[k] = rmus*-1
        if k in names:
            NAME = "_".join(names[k].split())
        else:
            NAME = "Unknown"
        printlist[total_reads].append([NAME, k, total_reads, rmus*-1, ' | '.join(
            [str(x) for x in [str(bin_size*y)+':'+str(z) for y, z in BinDict[k].items()]])])

    # Sort the printlist by total_reads in descending order
    sorted_reads = sorted(printlist.keys(), reverse=True)
    for reads in sorted_reads:
        for entry in printlist[reads]:
            OUT.write("\t".join([str(x) for x in entry]) + "\n")
    OUT.close()
    return RMUS


parents = {}
rank_dict = {}
names = {}

# Load the taxonomy data from nodes.dmp and names.dmp files
with load_data(options.Tax) as node:

    for line in node:
        # Split the line by tab and remove leading/trailing whitespace

        its = [x.replace("\t", "") for x in line.rstrip().split("|")]
        parents[its[0]] = its[1]
        rank_dict[its[0]] = its[2]

with load_data(options.Names) as name:
    for line in name:
        its = [x.replace("\t", "") for x in line.rstrip().split("|")]
        if "scientific name" in line:
            names[its[0]] = its[1]

RMUS = calculate_rmus_from_paf(options.PAF, int(
    options.bins), open(options.OUT+".RMUS.txt", 'wt'))

with load_data(options.PAF) as PAF, open(options.OUT+".txt", 'wt') as export:
    # write header to output file
    export.write(
        "SeqID\tTaxID\tLength\tMappingQuality\tdomain\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
    for line in PAF:
        lines = line.rstrip().split("\t")
        if len(lines) < 2:
            continue
        seqId = lines[0]
        taxId = lines[5].split("|")[1]
        refLen = lines[6]
        refStart = lines[7]
        refEnd = lines[8]

        if RMUS[taxId] < float(options.RMUS):
            continue

        if taxId not in names:
            continue

        node_path, name_path = taxon_trace(taxId)
        # loop through the following ranks "no rank|cellular root|domain|kingdom|phylum|class|order|family|genus|species" and make new list from name path. replace with NA if not available
        RankName = dict(zip(node_path.split("|"), name_path.split("|")))
        name_path = []
        for rank in ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]:
            if rank in RankName:
                name_path.append(RankName[rank])
            else:
                name_path.append("NA")

        Tax = "\t".join(["_".join(x.split()) for x in name_path])
        export.write(seqId+"\t"+taxId+"\t" +
                     lines[1]+"\t"+lines[11] + "\t"+Tax+"\n")


sys.exit()
