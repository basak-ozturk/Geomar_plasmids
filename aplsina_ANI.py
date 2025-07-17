import pandas as pd
#import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt

#from scipy.cluster.hierarchy import fcluster
#from scipy.cluster.hierarchy import to_tree
# Load FastANI output
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/Kristina_fastani_output.txt", sep="\t", header=None)
df.columns = ["query", "ref", "ani", "shared_frags", "total_frags"]

# Clean filenames (remove paths and extensions)
df["query"] = df["query"].apply(lambda x: x.split("/")[-1].replace(".fasta", ""))
df["ref"] = df["ref"].apply(lambda x: x.split("/")[-1].replace(".fasta", ""))

def clean_name(name):
    base = name.split("/")[-1]
    for ext in [".fasta", ".fna", "_fna", "_fasta"]:
        base = base.replace(ext, "")
    return base

df["query"] = df["query"].apply(clean_name)
df["ref"] = df["ref"].apply(clean_name)

# All unique plasmids
plasmids = sorted(set(df["query"]).union(df["ref"]))
dist_matrix = pd.DataFrame(100.0, index=plasmids, columns=plasmids)
dist_matrix.index = dist_matrix.index.str.replace(r".fna$|_fasta$", "", regex=True)
dist_matrix.columns = dist_matrix.columns.str.replace(r".fna$|_fasta$", "", regex=True)


# Fill distances (100 - ANI)
for _, row in df.iterrows():
    dist = 100 - row["ani"]
    dist_matrix.loc[row["query"], row["ref"]] = dist
    dist_matrix.loc[row["ref"], row["query"]] = dist

# Save matrix
#dist_matrix.to_csv("C:/Users/hayat/Downloads/R_files/data/ani_distance_matrix.tsv", sep="\t")


# Convert to condensed format for linkage
condensed_dist = squareform(dist_matrix.values)
linkage_matrix = linkage(condensed_dist, method='average')

sns.set(style="white")
plt.figure(figsize=(10, 8))
g = sns.clustermap(dist_matrix,
    row_linkage=linkage_matrix,
    col_linkage=linkage_matrix, cmap="magma", xticklabels=True, yticklabels=True)
plt.suptitle("Clustered ANI Distance Matrix of the Aplysina Plasmids (100 - ANI%)", y=1.02)
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/Kristina_clustered_ani_heatmap_no_labels.png", dpi=300, bbox_inches='tight')
plt.show()

# Plot dendrogram
plt.figure(figsize=(10, 6))
dendrogram(linkage_matrix, labels=dist_matrix.index.tolist(), leaf_rotation=90)
plt.title("Plasmid Dendrogram (ANI Distance, UPGMA)")
#plt.xticks([], [])  # Removes both ticks and labels
plt.tight_layout()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/Kristina_plasmid_tree_ani_dendrogram.png", dpi=300)
plt.show()
