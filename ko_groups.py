import pandas as pd
#import numpy as np
import umap
import hdbscan
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import DBSCAN
import numpy as np
import umap.plot
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# Load data
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/KEGG_ko_per_plasmid.txt", sep="\t", header=None, names=["Plasmid", "KO"])

# Convert KO to list
df["KO"] = df["KO"].apply(lambda x: x.split(","))

# Create binary presence/absence matrix
binary_matrix = df.explode("KO").assign(Value=1).pivot(index="Plasmid", columns="KO", values="Value").fillna(0)

# Filter KO that appear in at least 10 plasmids
filtered_columns = binary_matrix.columns[binary_matrix.sum(axis=0) >= 10]
binary_matrix_filtered = binary_matrix[filtered_columns]

# HDBSCAN clustering before UMAP
hdbscan_model = hdbscan.HDBSCAN(min_cluster_size=8, min_samples=5, cluster_selection_epsilon=1.0)
hdbscan_labels = hdbscan_model.fit_predict(binary_matrix_filtered)

# Add HDBSCAN cluster labels to the data
binary_matrix_filtered["HDBSCAN_Cluster"] = hdbscan_labels

# UMAP dimensionality reduction with densmap
umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, metric='euclidean', random_state=42, densmap=True)
umap_embedding = umap_model.fit_transform(binary_matrix_filtered.drop(columns=["HDBSCAN_Cluster"]))

# Convert UMAP results into a DataFrame
umap_df = pd.DataFrame(umap_embedding, columns=["UMAP1", "UMAP2"], index=binary_matrix_filtered.index).reset_index()
umap_df["HDBSCAN_Cluster"] = hdbscan_labels

# Plot UMAP with HDBSCAN clusters
plt.figure(figsize=(7, 5))
sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", hue=umap_df["HDBSCAN_Cluster"].astype(str), palette="Dark2", alpha=0.7)
plt.title("HDBSCAN Clustering of UMAP Results")
plt.legend(title="Cluster")
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/UMAP_HDBSCAN_plot_ko.png", dpi=300, bbox_inches='tight')
plt.show()

# DBSCAN clustering on UMAP
dbscan_model = DBSCAN(eps=1, min_samples=8)
umap_df["DBSCAN_Cluster"] = dbscan_model.fit_predict(umap_df[["UMAP1", "UMAP2"]])

# Plot UMAP with DBSCAN clusters
plt.figure(figsize=(7, 5))
sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", hue=umap_df["DBSCAN_Cluster"].astype(str), palette="Dark2", alpha=0.7)
plt.title("DBSCAN Clustering of UMAP Results")
plt.legend(title="Cluster")
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/UMAP_DBSCAN_plot_ko_py.png", dpi=300, bbox_inches='tight')
plt.show()


# KDE Density plot with cluster colors
plt.figure(figsize=(7, 5))
sns.kdeplot(
    data=umap_df, x="UMAP1", y="UMAP2", hue=umap_df["DBSCAN_Cluster"].astype(str),
    fill=True, alpha=0.5, palette="viridis"
)

#plt.legend(title="Cluster", fontsize=8, title_fontsize=10)

# Title and save
plt.title("Density-based Visualization of DBSCAN Clustering")
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/UMAP_DBSCAN_Density.png", dpi=300, bbox_inches='tight')
plt.show()



# Get list of plasmids per cluster
plasmid_list = umap_df.groupby("DBSCAN_Cluster")["Plasmid"].apply(list).to_dict()

# Save plasmid list to text file
with open("C:/Users/hayat/Downloads/R_files/data/plasmids_per_cluster_ko_py.txt", "w") as f:
    for cluster, plasmids in plasmid_list.items():
        f.write(f"Cluster {cluster}:\n")
        f.write("\n".join(plasmids) + "\n\n")
