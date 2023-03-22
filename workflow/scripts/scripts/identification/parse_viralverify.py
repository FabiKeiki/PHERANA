#######Parse viral verify #############
# Assign scores to virus : 1


import pandas as pd

viralverify = pd.read_csv(snakemake.input[0], usecols = ['Contig name','Prediction','Score'],sep=',')

  
viralverify.columns = ['contig','prediction','score']
viralverify = viralverify[viralverify["prediction"] == "Virus"]


viralverify = viralverify[["contig", "score"]]
viralverify['tool'] ='viralverify'


viralverify.to_csv(snakemake.output[0],sep='\t', index=False)