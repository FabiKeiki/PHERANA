################################# filter_phages _bacteria ##########################
# This script parses the checkm output and filters the phages identified as bacteria
# based on the completeness score. Only phages with completeness score < 5 are kept.

import pandas as pd
import ast
from Bio import SeqIO


# Read the TSV file into a DataFrame
checkm = pd.read_csv(snakemake.input['checkm'], sep='\t', header=None)

# Parse checkm output
checkm = checkm.iloc[:, [0, -1]] 

# Rename the columns
checkm.columns = ['Phage', 'Properties']

# Convert string representation of dictionary to an actual dictionary
checkm['Properties'] = checkm['Properties'].apply(ast.literal_eval)

# Extract completeness and contamination from the Properties column
checkm['Completeness'] = checkm['Properties'].apply(lambda x: x.get('Completeness', None))
checkm['Contamination'] = checkm['Properties'].apply(lambda x: x.get('Contamination', None))
checkm['lineage'] = checkm['Properties'].apply(lambda x: x.get('marker lineage', None))

# Drop the original Properties column
checkm = checkm.drop('Properties', axis=1)

# export checkm dataframe to tsv
checkm.to_csv(snakemake.output['checkm_parsed'], sep='\t', index=False)

# kkep only entries >5 completeness
checkm_cont = checkm[checkm['Completeness'] >= 5]

record_dict = SeqIO.to_dict(SeqIO.parse(snakemake.input['phages'],"fasta"))

#remove phages entries in fasta file where completeness is >5
for index, rows in checkm_cont.iterrows():
    record_dict.pop(rows['Phage'], None)

#remove phages entries in fasta file where lineage is not root
for index, rows in checkm.iterrows():
    if rows['lineage'] != 'root':
        record_dict.pop(rows['Phage'], None)
        
# write the new phages dict to a fasta file
SeqIO.write(record_dict.values(), snakemake.output['phages_filtered'], "fasta")
