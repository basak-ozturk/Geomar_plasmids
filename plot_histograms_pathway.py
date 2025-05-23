import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_histograms(tsv_file, output_dir):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Load the TSV file
    df = pd.read_csv(tsv_file, sep="\t", header=None, names=["Cluster", "KO", "Percentage", "Pathway"])

    # Convert data types
    df["Pathway"] = df["Pathway"].astype(str)  # Ensure pathway names are strings
    df["Percentage"] = pd.to_numeric(df["Percentage"], errors="coerce")  # Convert percentages to floats
    df = df.dropna(subset=["Percentage"])  # Drop rows with NaN percentages

    # Group by cluster and generate histograms
    for cluster, cluster_df in df.groupby("Cluster"):
        cluster_df = cluster_df.sort_values(by="Percentage", ascending=True)  # Ensure correct order

        plt.figure(figsize=(10, 6))
        plt.barh(cluster_df["Pathway"], cluster_df["Percentage"], color="skyblue")  # Use Pathway names

        plt.xlabel("Percentage")
        plt.ylabel("Pathway")
        plt.title(f"Pathway Distribution for Cluster {cluster}")

        plt.xticks(rotation=45, ha="right")  # Rotate x-axis labels for readability
        plt.tight_layout()

        # Save the figure
        output_path = os.path.join(output_dir, f"cluster_{cluster}.png")
        plt.savefig(output_path)
        plt.close()

# Example usage
plot_histograms("C:\\Users\\hayat\\Downloads\\R_files\\data\\eggnog_all_plasmids_pathway_per_plasmid.tsv", 
                "C:\\Users\\hayat\\Downloads\\R_files\\graphs\\pathways_names_per_cluster")


