# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 17:52:37 2025

@author: hayat
"""

import pandas as pd

# Load integrase list
with open("C:/Users/hayat/Downloads/R_files/data/integrase_proteins_widepsread.txt") as f:
    integrase_ids = [line.strip().rsplit('_', 1)[0] for line in f]

# Make a set for quick lookup and remove duplicates
integrase_plasmids = set(integrase_ids)

# Load your main dataframe
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_sequence_lengths.csv")

# Assume your plasmid names are in a column named "plasmid"
df["Integrase"] = df["Plasmid"].apply(lambda x: "+" if x in integrase_plasmids else "-")

# Save to a new file
df.to_csv("C:/Users/hayat/Downloads/R_files/data/widespread_with_integrase_and_length.csv", index=False)
