############################## Trim provirus ##############################
# This script is used to trim the proviruses sequences from the phages bins 
#fasta file using checkv proviruses fasta file
# All proviruses sequences belonging to the same phage are merged into one 
# sequence
###########################################################################
from Bio import SeqIO

provir = SeqIO.to_dict(SeqIO.parse(snakemake.input['provir'],"fasta"))
phages = SeqIO.to_dict(SeqIO.parse(snakemake.input['phages'],"fasta"))

#loop over the proviruses dict and keep only what's before the last underscore in the ID
for key, value in provir.items():
    provir[key].id = value.id.rsplit('_',1)[0]

# duplicate IDs are to be merged into one sequence
# create a new dict with unique IDs and merged sequences
provir_unique = {}
for key, value in provir.items():
    if value.id not in provir_unique:
        provir_unique[value.id] = value
    else:
        provir_unique[value.id].seq = provir_unique[value.id].seq + value.seq

# replace phages sequences value with proviruses sequences when they have the same ID
for key, value in phages.items():
    if value.id in provir_unique:
        phages[value.id].seq = provir_unique[value.id].seq

# write the new phages dict to a fasta file
SeqIO.write(phages.values(), snakemake.output[0], "fasta")