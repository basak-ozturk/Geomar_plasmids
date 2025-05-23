#import pandas as pd
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

# Load COG annotations
def load_COG_annotations(COG_file):
    COG_annotations = {}
    with open(COG_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                COG_annotations[parts[0]] = parts[1].split(',')
    return COG_annotations

# Count COG occurrences per cluster
def count_COG_per_cluster(clusters, COG_annotations):
    COG_counts = {}
    COG_percentages = {}
    
    for cluster, plasmids in clusters.items():
        all_COGs = []
        for plasmid in plasmids:
            if plasmid in COG_annotations:
                all_COGs.extend(COG_annotations[plasmid])
        
        total_plasmids = len(plasmids)
        COG_freq = Counter(all_COGs)
        
        # Compute percentage occurrence
        COG_percent = {COG: (count / total_plasmids) * 100 for COG, count in COG_freq.items()}
        
        # Filter COGs that occur in at least 10% of plasmids
        frequent_COGs = {COG: perc for COG, perc in COG_percent.items() if perc >= 10}
        
        COG_counts[cluster] = frequent_COGs
        COG_percentages[cluster] = COG_percent
    
    return COG_counts, COG_percentages

# Save COG percentages to a file
def save_COG_percentages(COG_percentages, output_file):
    with open(output_file, 'w') as f:
        for cluster, percentages in COG_percentages.items():
            for COG, perc in percentages.items():
                f.write(f"{cluster}\t{COG}\t{perc:.2f}\n")

# Plot histograms
def plot_histograms(COG_counts):
    for cluster, COG_freq in COG_counts.items():
        plt.figure(figsize=(10, 5))
        plt.bar(COG_freq.keys(), COG_freq.values())
        plt.xticks(rotation=90)
        plt.xlabel("COG Number")
        plt.ylabel("Percentage of Plasmids")
        plt.title(f"COG Distribution in Cluster {cluster}")
        plt.show()

cluster_file = "C:\\Users\\hayat\\Downloads\\R_files\\data\\plasmids_per_cluster_COG_std.txt"  
COG_file = "C:\\Users\\hayat\\Downloads\\R_files\\data\\COGs_per_plasmid.txt"
output_file = "C:\\Users\\hayat\\Downloads\\R_files\\data\\COG_percentages_per_cluster.txt"

# Run analysis
clusters = load_clusters(cluster_file)
COG_annotations = load_COG_annotations(COG_file)
COG_counts, COG_percentages = count_COG_per_cluster(clusters, COG_annotations)
save_COG_percentages(COG_percentages, output_file)
plot_histograms(COG_counts)


