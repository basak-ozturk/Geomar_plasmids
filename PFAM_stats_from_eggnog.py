# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 08:33:13 2025

@author: hayat
"""

import os
import pandas as pd
from collections import Counter

cluster_folder = "C:/Users/hayat/Downloads/R_files/data/ani_cluster_annotations"
pfam_stats_output = os.path.join(cluster_folder, "cluster_pfam_stats.tsv")

stats = []

for fname in os.listdir(cluster_folder):
    if fname.startswith("cluster_") and fname.endswith(".tsv"):
        cluster_id = fname.replace("cluster_", "").replace(".tsv", "")
        path = os.path.join(cluster_folder, fname)
        
        df = pd.read_csv(path, sep="\t")
        
        if "PFAMs" not in df.columns and "PFAMs_list" not in df.columns:
            print(f"[Warning] No PFAM column in cluster {cluster_id}")
            continue
        
        # Use whichever column exists
        pfam_col = "PFAMs_list" if "PFAMs_list" in df.columns else "PFAMs"
        
        # Parse PFAMs: convert strings like "['PF00193', 'PF00270']" into lists
        pfam_lists = df[pfam_col].dropna().apply(eval)  # careful: assumes stringified list
        
        # Flatten all PFAMs
        all_pfam = [pfam for sublist in pfam_lists for pfam in sublist]
        pfam_counter = Counter(all_pfam)
        
        stats.append({
            "Cluster": cluster_id,
            "Proteins_with_Pfam": len(pfam_lists),
            "Total_Pfam_Count": len(all_pfam),
            "Unique_Pfam_Count": len(set(all_pfam)),
            "Top_10_Pfams": ", ".join([f"{k}({v})" for k, v in pfam_counter.most_common(10)])
        })

# Save summary table
stats_df = pd.DataFrame(stats)
stats_df.sort_values("Cluster", inplace=True)
stats_df.to_csv(pfam_stats_output, sep="\t", index=False)
