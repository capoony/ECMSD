import sys
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
parser.add_option("--logical", dest="log",
                  help="logical parameter", action="store_true")
parser.add_option("--param", dest="param",
                  help="numerical parameter", default=1)

(options, args) = parser.parse_args()
parser.add_option_group(group)

#


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


# Create a dictionary to map names to taxids
# The dictionary will be used to rename the FASTA headers
# to include the taxid in the format: >kraken:taxid|<taxid>
# The taxid will be used to filter the sequences based on the taxid
TaxidDict = d(str)
for l in load_data(options.TX):
    # Split the line into name and taxid
    # The input file is expected to be a tab-separated file with two columns:
    # name<tab>taxid
    name, taxid = l.rstrip().split("\t")
    TaxidDict[name] = taxid

# Create a set to keep track of taxids that have already been written
# This will prevent duplicate entries in the output file
# The set will be used to filter the sequences based on the taxid
# If a taxid is already in the set, it will be skipped
# If a taxid is not in the set, it will be added to the set and written to the output file
IDlist = d(str)
# Initialize a flag to skip lines that do not match the taxid
SKIP = False
# Read the input FASTA file and rename the headers
with load_data(options.IN) as f, open(options.OUT, "wt") as o:
    # Iterate through each line in the input file
    for line in f:
        # Check if the line starts with '>' indicating a FASTA header
        if line.startswith(">"):
            # Extract the name from the header line
            name = line[1:].rstrip().split(" ")[0]
            # Check if the name is in the TaxidDict and if the taxid is not already in IDlist
            if name in TaxidDict and TaxidDict[name] not in IDlist:
                # If the taxid is not in IDlist, add it to IDlist and write the new header
                SKIP = False
                IDlist[TaxidDict[name]]
                o.write(f">kraken:taxid|{TaxidDict[name]}\n")
                continue
            # If the name is not in TaxidDict or the taxid is already in IDlist, skip the line
            else:
                SKIP = True
        # If the line does not start with '>', it is a sequence line
        # Write the sequence line only if SKIP is False
        # This ensures that only sequences with valid taxids are written to the output file
        # If SKIP is True, the sequence line will be skipped
        # This prevents writing sequences that do not match the taxid criteria
        if SKIP == False:
            o.write(line)
