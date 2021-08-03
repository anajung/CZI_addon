#!/usr/bin/env python

import pandas as pd
#import sys

def join_lineage(pangolin_lineage, nextclade_lineage):
    pangolin = pd.read_csv(pangolin_lineage)
    nextclade = pd.read_csv(nextclade_lineage, sep = '\t')
    pangolin.rename(columns={'taxon': 'Sample_Name','lineage':'Pangolin_lineage'},inplace=True)
    nextclade.rename(columns={'seqName': 'Sample_Name','clade':'Nextclade_lineage'},inplace=True)
    final_lineage = pd.merge(left=pangolin,right=nextclade,left_on='Sample_Name',right_on='Sample_Name')
    final_lineage.insert(2, 'Nextclade_lineage', final_lineage.pop("Nextclade_lineage"))
    filtered_lineage = final_lineage[final_lineage.Pangolin_lineage != 'None']
    final_lineage.to_csv('joined_lineage.csv',index=False)
    filtered_lineage.to_csv('filtered_joined_lineage.csv',index=False)

pangolin_lineage = "!{pangolin_lineage}"
nextclade_lineage = "!{nextclade_lineage}"
join_lineage(pangolin_lineage, nextclade_lineage)