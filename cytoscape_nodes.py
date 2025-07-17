# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 10:13:02 2025

@author: hayat
"""

import pandas as pd

# === Load main RPKM + metagenome data ===
main_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/top_widespread_plasmids_with_mean_total_rpkm.csv")  # Your file with Mean_RPKM, Total_RPKM, etc.

# === Load Mob-positive plasmid list ===
# Either CSV or TXT with one column of plasmid names
mob_file = "C:/Users/hayat/Downloads/R_files/data/mobc_plasmids.txt"
mob_df = pd.read_csv(mob_file, header=None, names=["Plasmid"])
mob_set = set(mob_df["Plasmid"])
main_df["Mob"] = main_df["Plasmid"].apply(lambda x: "yes" if x in mob_set else "no")

int_file = "C:/Users/hayat/Downloads/R_files/data/integrase_plasmids.txt"
int_df = pd.read_csv(mob_file, header=None, names=["Plasmid"])
int_set = set(mob_df["Plasmid"])
main_df["Integrase"] = main_df["Plasmid"].apply(lambda x: "yes" if x in mob_set else "no")

# === Load Host count file ===
host_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_host_counts.csv")  # Has columns: Plasmid, Host_Count
main_df = main_df.merge(host_df, on="Plasmid", how="left")

# Fill missing host counts with 0 (if any plasmids not in host list)
main_df["Host_Count"] = main_df["Host_Count"].fillna(0).astype(int)

# === Rename Plasmid column to ID for Cytoscape ===
main_df = main_df.rename(columns={"Plasmid": "ID"})

# === Save final node table ===
main_df.to_csv("C:/Users/hayat/Downloads/R_files/data/cytoscape_nodes.csv", index=False)

print("✅ Saved: cytoscape_nodes.csv")
