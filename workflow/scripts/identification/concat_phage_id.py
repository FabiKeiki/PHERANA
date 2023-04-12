################### concat phage id ##################
import pandas as pd
import re

phid_list = snakemake.input
phid_dfs= []

for phid in phid_list:
    df = pd.read_csv(phid, sep='\t')
    sample= re.search(r'\d+[B|P]', phid).group()
    df['sample'] = sample
    phid_dfs.append(df)
    

phid_all = pd.concat(phid_dfs)

phid_all.to_csv(snakemake.output[0],sep="\t",index=False)