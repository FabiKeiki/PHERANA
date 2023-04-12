################################# output phage fasta ##########################
#This script extract the phages sequences from the assemblies. A filter based
# on the confidence is applied : only sequences with confidence_thres >= x,
# will be selected, (x ranging from 1 to 3)

from Bio import SeqIO
import pandas as pd

phages = pd.read_csv(snakemake.input['phages_id'],sep='\t')


# filter phages according to score 
phages = phages[phages['confidence'].astype(int) >= snakemake.params['confidence_thres']]

# Load the names of the contigs from DataFrame A
contigs = phages['contig'].to_list()
print(len(contigs))

phages_seq = {}

# Iterate through the records in the input FASTA file
for record in SeqIO.parse(snakemake.input['concat_assembly'], "fasta"):
    # Check whether the record ID is in the set of contig names
    if record.id in contigs:
        # Add the selected record to the dictionary, with the ID as the key
        phages_seq[record.id] = record

# Write the selected records to the output file
with open(snakemake.output[0], "w") as output:
    SeqIO.write(phages_seq.values(), output, "fasta")

