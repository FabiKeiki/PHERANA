################################## get_phage_bin #####################################################
#This script has mutltple purposes:
# (1) to export a files containing stats on the bins (number of contigs, length of the bin, sample name)
# (2) to add phamb for single contigs bins to the contigs scoring file
# (3) to extract two files: one containing info for single contigs bins and one for multi contigs bins
#######################################################################################################
import pandas as pd

# import vamb,phamb and contigs_score 

vamb = pd.read_csv(snakemake.input['vamb'],sep='\t',names=['binname','contig'])
phamb = pd.read_csv(snakemake.input['phamb'],sep='\t')
contigs_score = pd.read_csv(snakemake.input['score'],sep='\t')


#change contigs score tools to numeric object
contigs_score.iloc[:,2:] = contigs_score.iloc[:,2:].apply(pd.to_numeric,errors='coerce')
# change phamb probability to numeric object
phamb['probability'] = phamb['probability'].apply(pd.to_numeric,errors='coerce')


print(contigs_score.dtypes)
print(phamb.dtypes)

#(1) Exporting vamb and phamb stats

## merge phamb and vamb on binname
merged = pd.merge(phamb,vamb,on='binname',how='left')

## add a column with sample name.sample name is in the contig column before the "C"
merged['sample'] = merged['contig'].str.split("C").str[0]

## Trim the sample name from contigs name
merged['contig'] = merged['contig'].str.split("C").str[1]

## create a new df where each row is a bin and add a columm with the list of contigs belonging to this bin
merged = merged.groupby(['binname','sample','label','probability'])['contig'].apply(list).reset_index(name='contig')

## loop trough the contig column and add a column with summed length of the contigs for each bin
merged['length'] = merged['contig'].apply(lambda x: sum([int(i.split("_")[3]) for i in x]))


## export merged without the probibality column to csv
merged_out = merged[['binname','sample','label','length','contig']]
merged_out.to_csv(snakemake.output['bins'],sep='\t',index=False)

# Getting phages_id : bins and single contigs

## new df with bins of label viral
phages= merged[merged['label'] == 'viral']

# Keep only bins with porbability > 0.8
phages = phages[phages['probability'] > 0.8]


#(2) add single contigs bins to the contigs scoring file

## create new phages df with bins only conating one contig and unlist contig column
single_contigs = phages[phages['contig'].str.len() == 1]
single_contigs['contig'] = single_contigs['contig'].str[0]


## rename probability column to phamb and only keep sample,score and contig columns
single_contigs = single_contigs.rename(columns={'probability':'phamb'})
single_contigs = single_contigs[['sample','contig','phamb']]


## merge single contigs with contigs_score on contig and sample
phage_cntg = pd.merge(single_contigs,contigs_score,on=['sample','contig'],how='right')

# # convert all tools columns to float
# phage_cntg.iloc[:,2:] = phage_cntg.iloc[:,2:].astype('float')


## replace nan with 0
phage_cntg = phage_cntg.fillna(0)


## set score
score_3 = 1
score_2 = 0.8
score_1 = 0.5
score_0 = 0.3


## for all but the two first columns, add change the value with 3 if score >score_2, 2 if score > score_1 and 1 if score > score_0, 0 otherwise
phage_cntg.iloc[:,2:] = phage_cntg.iloc[:,2:].applymap(lambda x: 3 if x > score_2 else (2 if x > score_1 else (1 if x > score_0 else 0)))

## discard conitgs with contig score = 0
phage_cntg = phage_cntg[phage_cntg.iloc[:,2:].sum(axis=1) > 0]

## add a colomn confidence with the mean of all tools
phage_cntg['confidence'] = phage_cntg[['phamb','viralverify','vibrant','virsorter','deepvirfinder']].mean(axis=1)


## round confidence to the closest integer
phage_cntg['confidence'] = phage_cntg['confidence'].round(0)


##  Extract multi contigs bins from their list and add them to a new list
multi_contigs = []
for i in phages['contig']:
    if len(i) > 1:
        for j in i:
            multi_contigs.append(j)


print(phage_cntg['contig'])

## for each contigs in phages_cntg, test if it is in multi_contigs, if yes True, if no False to a new column called binned
phage_cntg['binned'] = phage_cntg['contig'].isin(multi_contigs)


#(3) extract two files: one containing info for single contigs bins and one for multi contigs bins

phage_cntg.to_csv(snakemake.output['phages_cntgs'],sep='\t',index=False)


## create a df containing only multi contigs bins
phages_bins = phages[phages['contig'].str.len() > 1]

phages_bins.to_csv(snakemake.output['phages_bins'],sep='\t',index=False)
