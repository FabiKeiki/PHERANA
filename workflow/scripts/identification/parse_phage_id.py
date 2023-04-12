#################### parse_identification ########################
# This script assemble all the phages id from different tools into
# a summary table. A confidence is then caclulated based on the
# score of each tool. Three level of confidence are calculated :
# 1: low
# 2: medium
# 3: high
# Viralverify score is set to one because only
# all the sequences correspond to identified viruses by the tool 
# (uncertains excluded,(see andrade-martinezComputationalToolsAnalysis2022) )
##################################################################

import pandas as pd
import re


# getting thresholds from snakemake
hscore_thres = 0.8
hcount_thres = 3

lscore_thres = 0.5
lcount_thres = 3

#importing input file (phages id)
virsorter = pd.read_csv(snakemake.input['virsorter'],sep='\t')
viralverify = pd.read_csv(snakemake.input['viralverify'],sep='\t')
deepvirfinder = pd.read_csv(snakemake.input['deepvirfinder'],sep='\t')
vibrant = pd.read_csv(snakemake.input['vibrant'],sep='\t')

phage_id = pd.concat([virsorter,viralverify,deepvirfinder,vibrant])
phage_id = phage_id.pivot_table(index= 'contig',columns='tool', values='score').fillna(0)
phage_id= phage_id.reset_index()

# setting new score for viralverify
phage_id['viralverify'] = 1


# counting number of tools passing score thresholds
phage_id['high_count'] = (phage_id.drop(['contig'],axis=1) > hscore_thres).sum(axis=1)
phage_id['low_count'] = (phage_id.drop(['contig','high_count'],axis=1) > lscore_thres).sum(axis=1)

# categorizing using a threshold for the minumum number of tools in passing thresholds
phage_id['confidence'] = 1
phage_id.loc[phage_id['high_count'] >= hcount_thres , 'confidence'] = 3 
phage_id.loc[(phage_id['confidence'] != 3) & (phage_id['low_count'] >= lcount_thres), 'confidence'] = 2


phage_id.drop(['high_count','low_count'],axis=1).to_csv(snakemake.output[0],sep="\t",index=False)



