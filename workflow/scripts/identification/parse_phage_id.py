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


sample = snakemake.wildcards.sample

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

# add a coulmn with the sample name at first position
phage_id.insert(0,'sample',sample)

phage_id.to_csv(snakemake.output[0],sep="\t",index=False)



