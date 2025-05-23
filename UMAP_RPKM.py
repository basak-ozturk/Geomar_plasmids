# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 09:51:09 2025

@author: hayat
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
#from sklearn.decomposition import PCA
import umap
#from sklearn.cluster import DBSCAN
import umap.plot

df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/CoverM_MAPPING_rpkm_Plasmid_Contigs_ouput.tsv",sep='\t', index_col=0)  
#cluster_list = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/plasmids_per_cluster_ko_std.txt", header=None)[0].tolist() 



with open("C:/Users/hayat/Downloads/R_files/data/plasmids_per_cluster_ko_std.txt", 'r') as f:
    sequence_names = [line.strip() for line in f if line.strip() and not line.startswith('Cluster')]

# Filter rows
filtered_df = df.loc[df.index.isin(sequence_names)]

# Filter columns
column_prefixes = [name.rsplit('_', 1)[0] for name in sequence_names]
filtered_df = filtered_df.loc[:, filtered_df.columns.str.split('_').str[0].isin(column_prefixes)]

# Step 5: Save the result to a new TSV file (optional)
filtered_df.to_csv('filtered_output.tsv', sep='\t')

# Print the result
print(filtered_df)


# Standardize the data
scaler = StandardScaler()
filtered_matrix_std = scaler.fit_transform(filtered_df)

# Apply UMAP
umap_model = umap.UMAP(n_neighbors=20, min_dist=0.1, metric='cosine', random_state=42, densmap=True)
umap_result = umap_model.fit_transform(filtered_matrix_std)
umap_df = pd.DataFrame(umap_result, columns=["UMAP1", "UMAP2"], index=filtered_df.index)

 #Apply DBSCAN clustering
# dbscan = DBSCAN(eps=2, min_samples=5)
# clusters = dbscan.fit_predict(filtered_matrix_std)
# umap_df["Cluster"] = clusters


# Function to calculate centroids
# def calculate_centroids(df):
#     return df.groupby("Cluster")[df.columns[:2]].mean()

# umap_centroids = calculate_centroids(umap_df)

# Plot UMAP
plt.figure(figsize=(8, 6))
sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", palette="tab10", legend=True)
# for cluster, (x, y) in umap_centroids.iterrows():
#      plt.text(x, y, str(cluster), fontsize=12, fontweight='bold', ha='center', va='center', bbox=dict(facecolor='white', edgecolor='black', alpha=0.6))
plt.title("UMAP of RPKM DATA")
#plt.legend(title="Cluster")
plt.show()

umap.plot.points(umap_model, values=np.arange(5890), theme='viridis')
