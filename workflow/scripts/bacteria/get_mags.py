################################# get_ins_fasta ##########################
#This script parse the checkm result to keep only the bins with completeness > 75%
##########################################################################

from Bio import SeqIO
import pandas as pd
import ast
import os

# if not existiong make the output directory
if not os.path.exists(snakemake.output[0]):
    os.makedirs(snakemake.output[0])


# Read the TSV file into a DataFrame
checkm = pd.read_csv(snakemake.input["checkm"], sep='\t', header=None)



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


# Keep only entries >75 completeness and < 10 contamination
checkm_cont = checkm[(checkm['Completeness'] >= 75) & (checkm['Contamination'] <= 10)]


# Move all fasta files matching the bins in the checkm_cont dataframe to a new folder
for index, rows in checkm_cont.iterrows():
    print("cp " + snakemake.input["bins"]+"/"+rows['Phage']+".fasta " + snakemake.output[0])
    os.system("cp " + snakemake.input["bins"]+"/"+rows['Phage']+".fasta " + snakemake.output[0])