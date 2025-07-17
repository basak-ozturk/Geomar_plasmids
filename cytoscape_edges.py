import pandas as pd

# Load FastANI output
fastani_file = "C:/Users/hayat/Downloads/R_files/data/fastani_output.txt"
df = pd.read_csv(fastani_file, sep="\t", header=None,
                 names=["Query", "Reference", "ANI", "Fragments", "Total"])

def clean_name(path):
    return path.split("/")[-1].replace(".fasta", "").replace(".fa", "").replace(".fna", "")

df["Query"] = df["Query"].apply(clean_name)
df["Reference"] = df["Reference"].apply(clean_name)

df = df[df["Query"] != df["Reference"]]

df["Pair"] = df.apply(lambda x: tuple(sorted([x["Query"], x["Reference"]])), axis=1)
df = df.drop_duplicates("Pair")

df_edges = df[["Query", "Reference", "ANI"]].rename(columns={
    "Query": "Source", "Reference": "Target"
})
df_edges.to_csv("C:/Users/hayat/Downloads/R_files/data/cytoscape_edges.csv", index=False)

print("Saved: cytoscape_edges.csv")

