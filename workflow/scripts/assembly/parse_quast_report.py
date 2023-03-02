################# parse_quast_report   ##################
# This script takes a few important variables out of the 
# quast report and assemble all samples in a single file#

import pandas as pd
import re

qrs_list = snakemake.input

#concat all reports together in one panda df and  keeping column of interest
qrs_all = pd.concat(pd.read_csv(qrs,sep='\t',\
    usecols = ['Assembly','# contigs (>= 0 bp)','# contigs (>= 1000 bp)','Total length','N50','N75','L50','L75'])\
     for qrs in qrs_list)

qrs_all.columns=['assembly','contigs','contigs_1000','length','N50','N75','L50','L75']


# seperating assembly field (wich contains both sample and assembler) into two fields


sample = [s.split("_")[0] if\
          len(s.split("_")) ==3 else\
          s.split("_")[0] + "_filt"  \
          for s in qrs_all['assembly'].to_list()]
assembly = [s.split("_")[1] for s in qrs_all['assembly'].to_list()]

qrs_all.insert(0,'sample',sample,True)
qrs_all['assembly'] = assembly

qrs_all.to_csv(snakemake.output[0],sep="\t",index=False)

