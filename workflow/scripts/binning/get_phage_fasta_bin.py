################################# output bin phage fasta ##########################
#This script extract the phages sequences from the bins identidied as phages. and
# concatenate the sequences of the same bin into one sequence

from Bio import SeqIO
import pandas as pd
import ast


# This function is used to convert the string representation of a list to a list
def convert_contig(value):
    return ast.literal_eval(value)

bins = pd.read_csv(snakemake.input['phages_id'],sep='\t', converters={'contig': convert_contig})

# loop over contig list of each row and add sample name + "C" to each contig in the list
for index, rows in bins.iterrows():
    contigs = rows['contig']
    for i in range(len(contigs)):
        contigs[i] = rows['sample'] + "C" + contigs[i]
    bins.at[index,'contig'] = contigs

record_dict = SeqIO.to_dict(SeqIO.parse(snakemake.input['concat_assembly'],"fasta"))



bin_fasta = {}

for index, rows in bins.iterrows():
    for contig in rows['contig']:
        bin_fasta[contig] = record_dict[contig]
        #change contig name and id to bin name¨
        bin_fasta[contig].id = rows['binname']
        bin_fasta[contig].description = ""

# concatenate similar ids into one sequence
bin_fasta_unique = {}
for key, value in bin_fasta.items():
    if value.id not in bin_fasta_unique:
        bin_fasta_unique[value.id] = value
    else:
        bin_fasta_unique[value.id].seq = bin_fasta_unique[value.id].seq + value.seq

# write the new phages dict to a fasta file
SeqIO.write(bin_fasta_unique.values(), snakemake.output[0], "fasta")


