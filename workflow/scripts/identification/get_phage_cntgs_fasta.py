################################# output single contig phage fasta ##########################
#This script extract the phages sequences from the assemblies. A filter based
# on the confidence is applied : only sequences with confidence_thres >= x,
# will be selected, (x ranging from 1 to 3)

from Bio import SeqIO
import pandas as pd

phages = pd.read_csv(snakemake.input['phages_id'],sep='\t')


# filter phages according to score 
phages = phages[phages['confidence'].astype(int) >= 1]

# #keep only non-bineed contigs (binned = False)
phages = phages[phages['binned'] == False]

# add sample name + "C" to contig name
phages['contig'] = phages['sample'] + "C" + phages['contig']

# # Load the names of the contigs from DataFrame A
contigs = phages['contig'].to_list()
print(contigs)

phages_seq = {}

record_dict = SeqIO.to_dict(SeqIO.parse(snakemake.input['concat_assembly'],"fasta"))
    
# get all fasta sequences of contigs from the assemblies
for contig in contigs:
   phages_seq[contig] = record_dict[contig]

with open(snakemake.output[0], "w") as output:
    SeqIO.write(phages_seq.values(), output, "fasta")