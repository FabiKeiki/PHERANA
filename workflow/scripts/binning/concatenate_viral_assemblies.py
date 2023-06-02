
import sys
import os
import argparse
from Bio import SeqIO


# This script parses the viral asemblies for VAMB
# It creaes a concat file with contigs renamed in the VAMB fashion and a metadata table of all contigs and wether they are retained or no
# 
# ​

# This script parses the viral asemblies for VAMB
# It creaes a concat file with contigs renamed in the VAMB fashion and a metadata table of all contigs and wether they are retained or no
# 
# ​
parser = argparse.ArgumentParser(
    description="""Creates the input FASTA file for Vamb.
Input should be one or more FASTA files, each from a sample-specific assembly.
If keepnames is False, resulting FASTA can be binsplit with separator 'C'.""",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=False,
)

parser.add_argument("outpath", help="Path to output FASTA file")
parser.add_argument("inpaths", help="Paths to input FASTA file(s)", nargs="+")
args = parser.parse_args()

concat= open(args.outpath, "w+")

for file in args.inpaths:
    sample_file = os.path.basename(file)
    sam = sample_file.split('_')[0]
    print("\nParsing " + sam)
    for entry in SeqIO.parse(file, "fasta"):
        length=len(entry.seq)
        print("-------------sequence: " + str(entry.id) + " Length="+str(length))
        
        if length>=2000:
            old=entry.id
            new=sam + "C" + old
            concat.write(">" + new + "\n" + str(entry.seq) +"\n")
concat.close()