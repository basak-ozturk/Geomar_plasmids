

import pandas as pd

# === Load plasmid lists ===
geodia_plasmids = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/plasmids_in_Geodia.tsv", sep="\t")["Plasmid"]
agelas_plasmids = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/plasmids_in_Agelas.tsv", sep="\t")["Plasmid"]

geodia_set = set(geodia_plasmids)
agelas_set = set(agelas_plasmids)

# === Read eggNOG annotations ===
eggnog_file =  "C:/Users/hayat/Downloads/R_files/data/eggnog_output.emapper.annotations"  

# Step 1: Skip metadata lines, but keep the real header
with open(eggnog_file) as f:
    lines = f.readlines()

# Keep lines that do NOT start with '##'
data_lines = [line for line in lines if not line.startswith("##")]

# Step 2: Now read it into pandas
from io import StringIO
annotations = pd.read_csv(StringIO("".join(data_lines)), sep="\t")

# Step 3: Extract plasmid name from protein ID
annotations["Plasmid"] = annotations["#query"].str.extract(r"^(.+_\d+ctg)")

# Step 4: Filter
geodia_annots = annotations[annotations["Plasmid"].isin(geodia_set)].copy()
agelas_annots = annotations[annotations["Plasmid"].isin(agelas_set)].copy()

# Step 5: Save results
geodia_annots.to_csv("C:/Users/hayat/Downloads/R_files/data/eggnog_annotations_Geodia.tsv", sep="\t", index=False)
agelas_annots.to_csv("C:/Users/hayat/Downloads/R_files/data/eggnog_annotations_Agelas.tsv", sep="\t", index=False)

print("Annotations saved for Geodia and Agelas.")
