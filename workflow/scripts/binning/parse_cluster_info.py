######################Parse cluster info #############
# this script parses the clustering of bins
# to generate a vamb like table ouptut with the clusters


import pandas as pd
print('importing input :')
file = snakemake.input[0]

print(file)
cluster = pd.read_csv(snakemake.input[0],sep='\t',header=None)

# keepn only rows with C in first column:
cluster = cluster[cluster[0].str.contains("C")]
print(cluster)
# define name of bin as the elements before "_clusters.tsv"
bin = str(file).split("_clusters_info.tsv")[0].split("/")[-1]
print('bin :' +bin)

# new df with bin name in first column and a list of clusters in second column : vae_8 [k126_1,k126]
cluster = pd.DataFrame([[bin,cluster[8].tolist()]],columns=['bin','clusters'])

# add third column with number of clusters
cluster['n_clusters'] = cluster['clusters'].str.len()

# export to tsv
cluster.to_csv(snakemake.output[0],index=False,sep='\t')

