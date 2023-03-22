####### parse_identification#############




import pandas as pd

virsorter = pd.read_csv('91P_virsorter_score_parsed.tsv',sep='\t')
viralverify = pd.read_csv('91P_viralverify_score_parsed.tsv',sep='\t')
deepvirfinder = pd.read_csv('91P_deepvirfinder_score_parsed.tsv',sep='\t')

phage_id = pd.concat([virsorter,viralverify,deepvirfinder])
phage_id = phage_id.pivot_table(index= 'contig',columns='tool', values='score').fillna(0)
phage_id= phage_id.reset_index()

phage_id.to_csv('91P_alltools_score.tsv',sep="\t",index=False)

