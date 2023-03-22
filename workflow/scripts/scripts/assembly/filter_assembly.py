###########################filter assembly############################
# This script parse the assembled contig files (for megahit and spades) 
#and applies a filter using the defined in parameters

import os
import sys
import re
from Bio import SeqIO

length_threshold = snakemake.params["length_t"]
coverage_threshold = snakemake.params["cov_t"]

contigs_filt=[]
contigs= SeqIO.parse(open(snakemake.input[0]),'fasta')

# looping through contigs in fasta file. The regexp is used to catch the coverage (multi in megahit, cov in spades)
for contig in contigs:
    coverage = float(re.findall(r"(?:cov_|multi=)([\d.]+)", contig.description)[0])
    length = len(contig.seq)
    if length>1000 and coverage>1:
        contigs_filt.append(contig)

SeqIO.write(contigs_filt, snakemake.output[0], "fasta")