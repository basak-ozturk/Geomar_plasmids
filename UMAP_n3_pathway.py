import pandas as pd
import umap
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN
import numpy as np
import matplotlib.colors as mcolors
# Load data
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/KEGG_pathway_per_plasmid.txt", sep="\t", header=None, names=["Plasmid", "Pathway"])

# Convert Pathway to list
df["Pathway"] = df["Pathway"].apply(lambda x: x.split(","))

# Create binary presence/absence matrix
binary_matrix = df.explode("Pathway").assign(Value=1).pivot(index="Plasmid", columns="Pathway", values="Value").fillna(0)

# Filter Pathways that appear in at least 10 plasmids
filtered_columns = binary_matrix.columns[binary_matrix.sum(axis=0) >= 10]
binary_matrix_filtered = binary_matrix[filtered_columns]

# Standardize the data
scaler = StandardScaler()
binary_matrix_std = scaler.fit_transform(binary_matrix_filtered)

# UMAP dimensionality reduction with 3 components
umap_model = umap.UMAP(n_components=3, n_neighbors=15, min_dist=0.1, metric='euclidean', random_state=42)
umap_embedding = umap_model.fit_transform(binary_matrix_std)



# Convert UMAP results into a DataFrame
umap_df = pd.DataFrame(umap_embedding, columns=["UMAP1", "UMAP2", "UMAP3"], index=binary_matrix_filtered.index)

dbscan_model = DBSCAN(eps=3, min_samples=10)
umap_df["DBSCAN_Cluster"] = dbscan_model.fit_predict(umap_df[["UMAP1", "UMAP2", "UMAP3"]])


# Get unique cluster labels
cluster_labels = sorted(umap_df["DBSCAN_Cluster"].unique())
n_clusters = len(cluster_labels)

# Create a color list: red for -1, dark blue for 0, cyan for 1, and a colormap for other clusters
colors = ['red', 'darkblue', 'cyan'] + [plt.cm.viridis(i) for i in np.linspace(0, 1, n_clusters-3)]
cmap = mcolors.ListedColormap(colors)

# Create color mapping
color_mapping = {label: idx for idx, label in enumerate(cluster_labels)}

# Map cluster labels to color indices
mapped_colors = umap_df["DBSCAN_Cluster"].map(color_mapping)

# Create a 3D scatter plot colored by DBSCAN clusters
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

scatter = ax.scatter(
    umap_df["UMAP1"],
    umap_df["UMAP2"],
    umap_df["UMAP3"],
    c=mapped_colors,
    cmap=cmap,
    s=20,
    alpha=0.7
)

# Add axis labels and title
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.set_zlabel("UMAP3")
plt.title("3D UMAP Visualization Colored by DBSCAN Clusters")

# Add color bar for clusters
cbar = plt.colorbar(scatter, ax=ax, shrink=0.5)
cbar.set_label("Cluster Label")
cbar.set_ticks(range(len(cluster_labels)))
cbar.set_ticklabels(cluster_labels)

# Save and show the plot
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/UMAP_3D_DBSCAN_plot.png", dpi=300, bbox_inches='tight')
plt.show()
