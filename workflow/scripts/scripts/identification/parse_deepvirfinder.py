#######Parse deepvirfinder#############
# Remove all p.values >0.05


import pandas as pd

deepvirfinder = pd.read_csv(snakemake.input[0], usecols = ['name','pvalue','score'],sep='\t')

deepvirfinder = deepvirfinder[deepvirfinder['pvalue'] < 0.05][['name', 'score']]
  
deepvirfinder.columns = ['contig','score']
deepvirfinder['tool'] = 'deepvirfinder'

deepvirfinder['contig'] = deepvirfinder['contig'].str.split(' ').str[0]
deepvirfinder.to_csv(snakemake.output[0],sep='\t', index=False)