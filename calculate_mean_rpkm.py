# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 07:09:26 2025
@author: hayat
"""

import pandas as pd

# --- File paths ---
plasmid_list_file = "C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_plasmid_names.txt"
rpkm_matrix_file = "C:/Users/hayat/Downloads/R_files/data/CoverM_MAPPING_rpkm_Plasmid_Contigs_ouput.tsv"
output_overview = "C:/Users/hayat/Downloads/R_files/data/top_widespread_plasmids_with_mean_total_rpkm.csv"
output_mean = "C:/Users/hayat/Downloads/R_files/data/widespread_mean_rpkm.csv"
output_total = "C:/Users/hayat/Downloads/R_files/data/widespread_total_rpkm.csv"
output_presence = "C:/Users/hayat/Downloads/R_files/data/widespread_presence_count.csv"

# --- Step 1: Read plasmid list (no header) ---
plasmid_list = pd.read_csv(plasmid_list_file, header=None, names=['Plasmid'])
plasmid_list['Plasmid'] = plasmid_list['Plasmid'].astype(str).str.strip()

# --- Step 2: Read RPKM matrix ---
rpkm_df = pd.read_csv(rpkm_matrix_file, sep="\t")
if 'Plasmid' not in rpkm_df.columns:
    rpkm_df.rename(columns={rpkm_df.columns[0]: 'Plasmid'}, inplace=True)
rpkm_df['Plasmid'] = rpkm_df['Plasmid'].astype(str).str.strip()

# --- Step 3: Filter and compute metrics ---
filtered_rpkm = rpkm_df[rpkm_df['Plasmid'].isin(plasmid_list['Plasmid'])].copy()
numeric_cols = filtered_rpkm.drop(columns=['Plasmid']).apply(pd.to_numeric, errors='coerce')

filtered_rpkm['Mean_RPKM'] = numeric_cols.mean(axis=1).round(2)
filtered_rpkm['Total_RPKM'] = numeric_cols.sum(axis=1).round(2)
filtered_rpkm['Presence_Count'] = (numeric_cols > 0).sum(axis=1)

# --- Step 4: Save full summary ---
output_df = filtered_rpkm[['Plasmid', 'Mean_RPKM', 'Total_RPKM', 'Presence_Count']]
output_df.to_csv(output_overview, index=False)

# --- Step 5: Save iTOL-compatible CSVs ---
output_df[['Plasmid', 'Mean_RPKM']].to_csv(output_mean, index=False)
output_df[['Plasmid', 'Total_RPKM']].to_csv(output_total, index=False)
output_df[['Plasmid', 'Presence_Count']].to_csv(output_presence, index=False)

print("✅ RPKM summary and individual iTOL files saved.")
