# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 11:56:13 2025

@author: hayat
"""

import pandas as pd

pfam_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/enriched_PFAMs_all_widespread_significant_with_annotation.tsv", sep="\t", quoting=3)  # quoting=3 means csv.QUOTE_NONE
cog_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/enriched_COG_categories__all_widespread_significant.tsv", sep="\t", quoting=3)
kegg_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/enriched_KEGG_Pathways_all_widespread_significant.tsv", sep="\t", quoting=3)

pfam_df = pfam_df.rename(columns={'COG_category': 'COG', 'KEGG_pathways': 'KEGG'})
cog_df = cog_df.rename(columns={'COG_category': 'COG'})
kegg_df = kegg_df.rename(columns={'KEGG_Pathway': 'KEGG'})

pfam_cog_merged = pd.merge(pfam_df, cog_df, on='COG', how='outer', suffixes=('_PFAM', '_COG'))

# Split KEGG pathways into lists
pfam_cog_merged['KEGG_list'] = pfam_cog_merged['KEGG'].str.split(',')

# Explode so each KEGG pathway has its own row
pfam_cog_expanded = pfam_cog_merged.explode('KEGG_list')

# Strip whitespace
pfam_cog_expanded['KEGG_list'] = pfam_cog_expanded['KEGG_list'].str.strip()

kegg_df = kegg_df.rename(columns={
    'Foreground_Count': 'Foreground_Count_KEGG',
    'Background_Count': 'Background_Count_KEGG',
    'Odds_Ratio': 'Odds_Ratio_KEGG',
    'P_Value_Adjusted': 'P_Value_Adjusted_KEGG'
})

final_merged = pd.merge(
    pfam_cog_expanded,
    kegg_df,
    left_on='KEGG_list',
    right_on='KEGG',
    how='left',
    suffixes=('', '_KEGG')
)


agg_cols = {
    'Foreground_Count_PFAM': 'first',
    'Background_Count_PFAM': 'first',
    'Odds_Ratio_PFAM': 'first',
    'P_Value_Adjusted_PFAM': 'first',
    'COG': 'first',
    'Foreground_Count_COG': 'first',
    'Background_Count_COG': 'first',
    'Odds_Ratio_COG': 'first',
    'P_Value_Adjusted_COG': 'first',
    'KEGG': lambda x: ','.join(x.dropna().unique()),
    'Foreground_Count_KEGG': 'max',
    'Background_Count_KEGG': 'max',
    'Odds_Ratio_KEGG': 'max',
    'P_Value_Adjusted_KEGG': 'min'
}

final_agg = final_merged.groupby('PFAMs').agg(agg_cols).reset_index()
