import sys
import os
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--Taxid", dest="TX", help="")
parser.add_option("--output", dest="OUT", help="Output file")

(options, args) = parser.parse_args()
parser.add_option_group(group)


def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    # if no file is given, read from STDIN
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


def load_taxid_dict(fh):
    """Parse a tab-separated name<tab>taxid stream into a name->taxid dict."""
    taxid_dict = d(str)
    for l in fh:
        name, taxid = l.rstrip().split("\t")
        taxid_dict[name] = taxid
    return taxid_dict


def rename_fasta(in_fh, out_fh, taxid_dict):
    """
    Rewrite a FASTA stream, replacing each header with its taxid (drops the
    sequence description). Keeps only the first sequence seen per taxid;
    later sequences mapping to an already-written taxid (or to a name not in
    taxid_dict) are dropped entirely.
    """
    seen_taxids = set()
    skip = False
    for line in in_fh:
        if line.startswith(">"):
            name = line[1:].rstrip().split(" ")[0]
            if name in taxid_dict and taxid_dict[name] not in seen_taxids:
                skip = False
                seen_taxids.add(taxid_dict[name])
                out_fh.write(f">{taxid_dict[name]}\n")
                continue
            else:
                skip = True
        if not skip:
            out_fh.write(line)


def run():
    if not options.IN or not options.OUT or not options.TX:
        parser.print_help()
        sys.exit(1)

    if not os.path.isfile(options.IN):
        print(f"Error: Input file {options.IN} does not exist.")
        sys.exit(1)
    if not os.path.isfile(options.TX):
        print(f"Error: Taxid file {options.TX} does not exist.")
        sys.exit(1)

    with load_data(options.TX) as tx_fh:
        taxid_dict = load_taxid_dict(tx_fh)

    with load_data(options.IN) as f, open(options.OUT, "wt") as o:
        rename_fasta(f, o, taxid_dict)


if __name__ == "__main__":
    run()
