################################# get_ins_fasta ##########################
#This script extract the bins sequesnces from the assemblies for MAGS selection
# we keep only bins with length > 500'000 bp
##########################################################################

from Bio import SeqIO
import pandas as pd
import ast
import os

directory = snakemake.output['outdir']
if not os.path.exists(directory):
    os.makedirs(directory)

# This function is used to convert the string representation of a list to a list
def convert_contig(value):
    return ast.literal_eval(value)

bins = pd.read_csv(snakemake.input['bins'],sep='\t', converters={'contig': convert_contig})

# loop over contig list of each row and add sample name + "C" to each contig in the list
for index, rows in bins.iterrows():
    contigs = rows['contig']
    for i in range(len(contigs)):
        contigs[i] = rows['sample'] + "C" + contigs[i]
    bins.at[index,'contig'] = contigs

# keep only bins with length > 500'000 bp
bins = bins[bins['length'] > 500000]
print(bins)

record_dict = SeqIO.to_dict(SeqIO.parse(snakemake.input['concat_assembly'],"fasta"))
bin_fasta = {}

for index, rows in bins.iterrows():
    bin_fasta = {}
    for contig in rows['contig']:
        bin_fasta[contig] = record_dict[contig]
    with open(snakemake.output['outdir']+'/'+ rows['binname']+".fasta", "w") as output:
        SeqIO.write(bin_fasta.values(), output, "fasta")

