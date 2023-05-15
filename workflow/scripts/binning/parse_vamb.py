############################################## parse_vamb   ##################
# This script generates cluster fasta files and a parsed version of the vamb output
# with the number of contigs per bin and the sample name
import os 
import pandas as pd
from Bio import SeqIO


# make directory to store the bins if it does not exist:
if not os.path.exists(snakemake.output['outdir']):
    os.makedirs(snakemake.output['outdir'])
print("making directory : " + snakemake.output['outdir'])
# modify the following line to add a header:

bins = pd.read_csv(snakemake.input[0],sep='\t', names=["bin","contig"],index_col=False)

# filter bins df to keep only rows with bins with more than one contig
bins = bins.groupby('bin').filter(lambda x: len(x) > 1)

# add a colum with sample name : sample name is in the contigs name (second column) before the "C"
bins['sample'] = bins.iloc[:,1].str.split("C").str[0]

# trim the sample name from contigs name
bins['contig'] = bins.iloc[:,1].str.split("C").str[1]

# contigs should be aggregated in a list for each bin and keep sample colum
bins = bins.groupby(['bin','sample'])['contig'].apply(list).reset_index(name='contigs')

# add a column with number of contigs
bins['n_contigs'] = bins['contigs'].str.len()

#export bins df to tsv
bins.to_csv(snakemake.output['vamb_parsed'],index=False,sep='\t')

# create loop through the df and for each bin loop through the list of contigs and print the bin name, sample name and contig name

for index, rows in bins.iterrows():
    bin_seq = {}
    record_dict = SeqIO.to_dict(SeqIO.parse(snakemake.params[0]+rows['sample']+"_concat_assembly.fasta", "fasta"))
    for contig in rows['contigs']:
        bin_seq[contig] = record_dict[contig]
    with open(snakemake.output['outdir']+'/'+ rows['bin']+".fasta", "w") as output:
        SeqIO.write(bin_seq.values(), output, "fasta")



