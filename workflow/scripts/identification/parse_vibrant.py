#######Parse vibrant #############
# add a score of one to all identified viruses ;


import pandas as pd

vibrant = pd.read_csv(snakemake.input[0])

vibrant.columns = ['contig']
vibrant['score'] = 1
vibrant['tool'] = 'vibrant'


vibrant.to_csv(snakemake.output[0],sep='\t', index=False)