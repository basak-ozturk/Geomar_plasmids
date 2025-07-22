# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 13:38:38 2025

@author: hayat
"""

import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

# Replace with your filename
filename = "C:/Users/hayat/Downloads/R_files/data/eggnog_annotations_agelas.tsv"

df = pd.read_csv(filename, sep="\t")

pfams_series = df['PFAMs'].dropna()

all_pfams = []
for pfam_list in pfams_series:
    pfams = pfam_list.split(",")
    # Filter out '-' entries
    pfams = [p for p in pfams if p != '-']
    all_pfams.extend(pfams)

pfam_counts = Counter(all_pfams)
top_10 = pfam_counts.most_common(10)

names, counts = zip(*top_10)

plt.figure(figsize=(10,6))
plt.barh(names[::-1], counts[::-1], color='skyblue')
plt.xlabel('Count')
plt.title('Top 10 PFAM domains (excluding "-")')
plt.tight_layout()
plt.show()