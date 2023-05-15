#######Parse virsorter #############
# Keep max score only, lt2genes (detected genes in sequences ) nan values are replaced with a score of zero


import pandas as pd

virsorter = pd.read_csv(snakemake.input[0], usecols = ['seqname','max_score'],sep='\t')
  
virsorter.columns = ['contig','score']
virsorter['contig'] = virsorter['contig'].str.replace(r':?\|\|(full|lt2gene)', '', regex=True)
virsorter = virsorter.fillna(0)
virsorter['tool'] = "virsorter"
print(virsorter)


virsorter.to_csv(snakemake.output[0],sep='\t', index=False)