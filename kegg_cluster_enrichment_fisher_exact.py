# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 09:51:08 2025

@author: hayat
"""

import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import os
import re
from collections import defaultdict



cluster_files = {
    'Cluster1': 'C:/Users/hayat/Downloads/R_files/data/pathway_cluster_1.txt',
    'Cluster2': 'C:/Users/hayat/Downloads/R_files/data/pathway_cluster_2.txt',
    'Cluster3': 'C:/Users/hayat/Downloads/R_files/data/pathway_cluster_3.txt',
    'Cluster4': 'C:/Users/hayat/Downloads/R_files/data/pathway_cluster_4.txt',
    'Cluster5': 'C:/Users/hayat/Downloads/R_files/data/pathway_cluster_5.txt'
}

cluster_plasmids = {}
for cname, fname in cluster_files.items():
    with open(fname) as f:
        cluster_plasmids[cname] = set(line.strip() for line in f)

df_ann = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/KEGG_pathway_per_plasmid_with_names.txt", sep='\t', header=None, names=["Plasmid", "Pathways"])

def parse_kos(pathway_str):
    if pd.isna(pathway_str):
        return []

    entries = pathway_str.split(',')
    parsed = []

    for entry in entries:
        entry = entry.strip()
        # Match cases like "DNA replication (ko03030)" or just "ko03030"
        match = re.match(r'^(?:(.*?)\s*\()?(ko\d{5})\)?$', entry)
        if match:
            desc, ko = match.groups()
            desc = desc.strip() if desc else ""
            parsed.append((desc, ko))

    return parsed


# Create plasmid → list of (description, ko)
plasmid_to_kos = defaultdict(list)

for _, row in df_ann.iterrows():
    ko_list = parse_kos(row["Pathways"])
    for desc, ko in ko_list:
        plasmid_to_kos[row["Plasmid"]].append((desc, ko))
        
   # Get unique KO + description pairs
all_kos = set()
for kos in plasmid_to_kos.values():
    all_kos.update(kos)

# Group by KO, but keep one description (arbitrary or most common)
ko_to_desc = {}
for desc, ko in all_kos:
    if ko not in ko_to_desc:
        ko_to_desc[ko] = desc  # could improve by counting frequency

# Count all KO terms in full set
global_ko_counts = defaultdict(int)
for kos in plasmid_to_kos.values():
    for _, ko in kos:
        global_ko_counts[ko] += 1

results = []
all_plasmids = set(plasmid_to_kos.keys())

for cluster_name, cluster_set in cluster_plasmids.items():
    background = all_plasmids - cluster_set

    # Count in cluster
    cluster_ko_counts = defaultdict(int)
    for plasmid in cluster_set:
        for _, ko in plasmid_to_kos.get(plasmid, []):
            cluster_ko_counts[ko] += 1

    # Count in background
    background_ko_counts = defaultdict(int)
    for plasmid in background:
        for _, ko in plasmid_to_kos.get(plasmid, []):
            background_ko_counts[ko] += 1

    total_fg = sum(cluster_ko_counts.values())
    total_bg = sum(background_ko_counts.values())

    for ko in global_ko_counts:
        a = cluster_ko_counts[ko]
        b = total_fg - a
        c = background_ko_counts[ko]
        d = total_bg - c
    
        if a < 10 or (a + c) < 5:
            continue
    
        table = [[a, b], [c, d]]
        oddsratio, pval = stats.fisher_exact(table, alternative='greater')
    
        # Compute 95% CI (approximate, log2 scale)
        or_log2 = np.log2(oddsratio) if oddsratio > 0 else 0
        se = np.sqrt(1/(a + 0.5) + 1/(d + 0.5))  # add pseudocount to avoid division by zero
        ci_low = or_log2 - 1.96 * se
        ci_high = or_log2 + 1.96 * se
    
        results.append({
            'Cluster': cluster_name,
            'KO': ko,
            'Description': ko_to_desc.get(ko, ""),
            'a': a, 'b': b, 'c': c, 'd': d,
            'OddsRatio': oddsratio,
            'log2(OR)': or_log2,
            'CI_lower': ci_low,
            'CI_upper': ci_high,
            'P': pval
        })


results_df = pd.DataFrame(results)
results_df['AdjP'] = results_df.groupby('Cluster')['P'].transform(lambda p: multipletests(p, method='fdr_bh')[1])

def forest_plot(df, cluster_name):
    sub = df[(df['Cluster'] == cluster_name) & (df['AdjP'] < 0.05)].copy()
    sub = sub[sub['a'] >= 10]  # Apply foreground count cutoff

    # Sort by enrichment strength (lowest adj P first)
    sub = sub.sort_values('AdjP').head(5)  # Only top 5
    sub = sub.sort_values('OddsRatio') 
    if sub.empty:
        print(f"No significant KO pathways in {cluster_name} after filtering.")
        return

    # Clean labels to avoid double parentheses
    sub['Label'] = sub.apply(
        lambda x: f"{x['Description'].strip()} ({x['KO']})" if x['Description'] else x['KO'],
        axis=1
    )

    plt.figure(figsize=(8, len(sub) * 0.6 + 2))  # slightly taller if needed
    plt.errorbar(
        sub['log2(OR)'], sub['Label'],
        xerr=[sub['log2(OR)'] - sub['CI_lower'], sub['CI_upper'] - sub['log2(OR)']],
        fmt='o', color='black', ecolor='gray', capsize=3
    )
    plt.axvline(0, color='red', linestyle='--')
    plt.xlabel('log2(Odds Ratio)')
    plt.title(f'Top {len(sub)} Enriched KEGG Pathways in {cluster_name}')
    plt.tight_layout()
    plt.savefig(f"C:/Users/hayat/Downloads/R_files/graphs/{cluster_name}_ko_forestplot.png")
    plt.show()
    plt.close()



for cluster in cluster_plasmids:
    forest_plot(results_df, cluster)

    
results_df.to_csv("ko_enrichment_with_CI.tsv", sep='\t', index=False)

summary = (
    results_df[results_df['AdjP'] < 0.05]
    .groupby('Cluster')
    .apply(lambda df: df[['KO', 'Description', 'OddsRatio', 'AdjP']].sort_values('AdjP'))
)

summary.to_csv("C:/Users/hayat/Downloads/R_files/data/ko_enrichment_summary_significant.tsv", sep='\t')

sig_dir = "C:/Users/hayat/Downloads/R_files/data/significant_kos_per_cluster"
os.makedirs(sig_dir, exist_ok=True)

for cluster in results_df['Cluster'].unique():
    sig = results_df[(results_df['Cluster'] == cluster) & (results_df['AdjP'] < 0.05)]
    if not sig.empty:
        sig.to_csv(f"{sig_dir}/{cluster}_significant.tsv", sep='\t', index=False)


