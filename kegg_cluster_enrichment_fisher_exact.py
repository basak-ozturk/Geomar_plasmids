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
#import os
import re
from collections import defaultdict
#import seaborn as sns



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

def combined_forest_plots(df, clusters, save_path):
    # Set layout (2 rows, 3 columns for 5 clusters)
    n = len(clusters)
    ncols = 3
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 6, nrows * 4), sharex=False)
    axes = axes.flatten()

    for idx, cluster_name in enumerate(clusters):
        ax = axes[idx]
        sub = df[(df['Cluster'] == cluster_name) & (df['AdjP'] < 0.05)].copy()
        sub = sub[sub['a'] >= 10]
        sub = sub.sort_values('AdjP').head(5).sort_values('OddsRatio')

        if sub.empty:
            ax.axis('off')
            ax.set_title(f'{cluster_name} (no sig. pathways)')
            continue

        sub['Label'] = sub.apply(
            lambda x: f"{x['Description'].strip()} ({x['KO']})" if x['Description'] else x['KO'],
            axis=1
        )

        ax.errorbar(
            sub['log2(OR)'], sub['Label'],
            xerr=[sub['log2(OR)'] - sub['CI_lower'], sub['CI_upper'] - sub['log2(OR)']],
            fmt='o', color='black', ecolor='gray', capsize=3
        )
        ax.axvline(0, color='red', linestyle='--')
        ax.set_title(cluster_name)
        ax.set_xlabel("log2(Odds Ratio)")
        if idx % ncols == 0:
            ax.set_ylabel("")  # Labels already shown

    # Hide any unused subplots
    for j in range(idx + 1, len(axes)):
        axes[j].axis('off')

    fig.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.show()
    plt.close()



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
   # plt.savefig(f"C:/Users/hayat/Downloads/R_files/graphs/{cluster_name}_ko_forestplot.png")
    plt.show()
    plt.close()



def unified_forest_plot(df, save_path):
    top_combined = []

    # Extract numeric cluster identifiers (e.g., from 'Cluster5' → 5)
    df = df.copy()
    df['ClusterNum'] = df['Cluster'].str.extract(r'(\d+)').astype(int)

    for cluster in df['ClusterNum'].unique():
        sub = df[(df['ClusterNum'] == cluster) & (df['AdjP'] < 0.05) & (df['a'] >= 10)].copy()
        sub = sub.sort_values('AdjP').head(5)  # Top 5
        top_combined.append(sub)

    combined_df = pd.concat(top_combined)

    # Clean labels
    combined_df['Label'] = combined_df.apply(
        lambda x: f"{x['Description'].strip()} ({x['KO']})" if x['Description'] else x['KO'],
        axis=1
    )

    # Sort by log2 OR for nicer display
    combined_df = combined_df.sort_values('log2(OR)')

    # Plot
    plt.figure(figsize=(10, len(combined_df) * 0.5 + 2))
    ax = plt.gca()

    # Custom cluster color mapping
    cluster_to_color = {
        1: "#2596be",
        2: "#e89e14",
        3: "#25be63",
        4: "#e85b14",
        5: "#ec75d8"
    }

    for idx, row in combined_df.iterrows():
        ax.errorbar(
            row['log2(OR)'],
            row['Label'],
            xerr=[[row['log2(OR)'] - row['CI_lower']], [row['CI_upper'] - row['log2(OR)']]],
            fmt='o',
            color=cluster_to_color.get(row['ClusterNum'], 'black'),  # fallback color if cluster not in map
            ecolor='gray',
            capsize=3,
            label=row['Cluster']
        )

    ax.axvline(0, color='red', linestyle='--')
    ax.set_xlabel("log2(Odds Ratio)")
    ax.set_title("Top Enriched KEGG Pathways by Cluster")

    # Custom legend without duplicates
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    
    # Sort by numeric cluster (assuming labels are like 'Cluster1', 'Cluster2', etc.)
    sorted_items = sorted(by_label.items(), key=lambda x: int(str(x[0]).strip().replace("Cluster", "")))
    sorted_labels, sorted_handles = zip(*sorted_items)
    
    ax.legend(sorted_handles, sorted_labels, title="Cluster")


    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.show()
    plt.close()

# Define numeric columns to round
numeric_cols = ['OddsRatio', 'log2(OR)', 'CI_lower', 'CI_upper', 'P', 'AdjP']
results_df[numeric_cols] = results_df[numeric_cols].round(2)

# Filter for OddsRatio > 1
results_df_filtered = results_df[results_df['OddsRatio'] > 1]

# Columns to include (no Description)
columns_for_table = ['Cluster', 'KO', 'OddsRatio', 'log2(OR)', 'CI_lower', 'CI_upper', 'AdjP']

# Collect top 10 per cluster
top10_list = []
for cluster in sorted(results_df_filtered['Cluster'].unique()):
    sub = results_df_filtered[results_df_filtered['Cluster'] == cluster].copy()
    sub = sub.sort_values('OddsRatio', ascending=False).head(10)
    top10_list.append(sub[columns_for_table])

# Combine into a single DataFrame
top10_combined = pd.concat(top10_list)

# Generate a single LaTeX longtable
latex_table_str = top10_combined.to_latex(
    index=False,
    caption="Top 10 enriched KOs per cluster (Odds Ratio > 1)",
    label="tab:kegg_enrichment",
    float_format="%.2f",
    longtable=True
)

# Print the LaTeX code
print(latex_table_str)



# Prepare a DataFrame to collect top 10 per cluster
top10_per_cluster = pd.DataFrame()

# Iterate over clusters
for cluster in sorted(results_df_filtered['Cluster'].unique()):
    sub = results_df_filtered[results_df_filtered['Cluster'] == cluster].copy()
    
    # Sort by OddsRatio descending
    sub = sub.sort_values('OddsRatio', ascending=False)
    
    # Take top 10
    sub_top = sub.head(10)
    
    # Append to the combined DataFrame
    top10_per_cluster = pd.concat([top10_per_cluster, sub_top[columns_for_table]])

# Save to CSV
top10_per_cluster.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/ko_enrichment_top10_per_cluster.csv",
    sep='\t',
    index=False
)

print("CSV saved with top 10 KOs per cluster.")





summary = (
    results_df[results_df['AdjP'] < 0.05]
    .groupby('Cluster')
    .apply(lambda df: df[['KO', 'Description', 'OddsRatio', 'AdjP']].sort_values('AdjP'))
)

#summary.to_csv("C:/Users/hayat/Downloads/R_files/data/ko_enrichment_summary_significant.tsv", sep='\t')

#sig_dir = "C:/Users/hayat/Downloads/R_files/data/significant_kos_per_cluster"
#os.makedirs(sig_dir, exist_ok=True)

# for cluster in results_df['Cluster'].unique():
#     sig = results_df[(results_df['Cluster'] == cluster) & (results_df['AdjP'] < 0.05)]
#     if not sig.empty:
#         sig.to_csv(f"{sig_dir}/{cluster}_significant.tsv", sep='\t', index=False)

# combined_forest_plots(
#     results_df,
#     clusters=list(cluster_plasmids.keys()),
#     #save_path="C:/Users/hayat/Downloads/R_files/graphs/combined_ko_forestplot_panel.png"
# )

#unified_forest_plot(
 #   results_df,
  #  save_path="C:/Users/hayat/Downloads/R_files/graphs/unified_ko_forestplot.svg"
#)

