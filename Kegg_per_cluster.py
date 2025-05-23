import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

# Load cluster data
def load_clusters(cluster_file):
    clusters = {}
    with open(cluster_file, 'r') as f:
        current_cluster = None
        for line in f:
            line = line.strip()
            if line.startswith("Cluster"):
                current_cluster = line.split()[1].strip(':')
                clusters[current_cluster] = []
            elif current_cluster is not None:
                clusters[current_cluster].append(line)
    return clusters

# Load KEGG annotations
def load_kegg_annotations(kegg_file):
    kegg_annotations = {}
    with open(kegg_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                kegg_annotations[parts[0]] = parts[1].split(',')
    return kegg_annotations

# Count KEGG occurrences per cluster
def count_kegg_per_cluster(clusters, kegg_annotations):
    kegg_counts = {}
    kegg_percentages = {}
    
    for cluster, plasmids in clusters.items():
        all_keggs = []
        for plasmid in plasmids:
            if plasmid in kegg_annotations:
                all_keggs.extend(kegg_annotations[plasmid])
        
        total_plasmids = len(plasmids)
        kegg_freq = Counter(all_keggs)
        
        # Compute percentage occurrence
        kegg_percent = {kegg: (count / total_plasmids) * 100 for kegg, count in kegg_freq.items()}
        
        # Filter KEGGs that occur in at least 10% of plasmids
        frequent_keggs = {kegg: perc for kegg, perc in kegg_percent.items() if perc >= 10}
        
        kegg_counts[cluster] = frequent_keggs
        kegg_percentages[cluster] = kegg_percent
    
    return kegg_counts, kegg_percentages

# Save KEGG percentages to a file
def save_kegg_percentages(kegg_percentages, output_file):
    with open(output_file, 'w') as f:
        for cluster, percentages in kegg_percentages.items():
            for kegg, perc in percentages.items():
                f.write(f"{cluster}\t{kegg}\t{perc:.2f}\n")

# Plot histograms
def plot_histograms(kegg_counts):
    for cluster, kegg_freq in kegg_counts.items():
        plt.figure(figsize=(10, 5))
        plt.bar(kegg_freq.keys(), kegg_freq.values())
        plt.xticks(rotation=90)
        plt.xlabel("KEGG Number")
        plt.ylabel("Percentage of Plasmids")
        plt.title(f"KEGG Distribution in Cluster {cluster}")
        plt.show()

# File paths
cluster_file = "C:\\Users\\hayat\\Downloads\\R_files\\data\\plasmids_per_cluster_ko_std.txt"  
kegg_file = "C:\\Users\\hayat\\Downloads\\R_files\\data\\KEGG_ko_per_plasmid.txt"
output_file = "C:\\Users\\hayat\\Downloads\\R_files\\data\\kegg_percentages_per_cluster.txt"

# Run analysis
clusters = load_clusters(cluster_file)
kegg_annotations = load_kegg_annotations(kegg_file)
kegg_counts, kegg_percentages = count_kegg_per_cluster(clusters, kegg_annotations)
save_kegg_percentages(kegg_percentages, output_file)
plot_histograms(kegg_counts)


