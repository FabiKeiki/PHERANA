####### output phage fasta#############

from Bio import SeqIO
import pandas as pd

phages = pd.read_csv('91P_alltools_score.tsv',sep='\t')

# Load the names of the contigs from DataFrame A
contigs = phages['contig'].to_list()
print(len(contigs))

phages_seq = {}

# Iterate through the records in the input FASTA file
for record in SeqIO.parse("91P_concat_assembly.fasta", "fasta"):
    # Check whether the record ID is in the set of contig names
    if record.id in contigs:
        # Add the selected record to the dictionary, with the ID as the key
        phages_seq[record.id] = record

# Write the selected records to the output file
with open("91P_phages.fasta", "w") as output_handle:
    SeqIO.write(phages_seq.values(), output_handle, "fasta")

