# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 14:41:42 2025

@author: hayat
"""

import pandas as pd

host_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/GEOLOCATION_HostGenus_Samples.txt", sep="\t")  # contains Run, biome_genus, Host
rpkm_matrix = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/sponge_plasmids_rpkm_submatrix_6945.tsv", sep="\t", index_col=0)  # plasmids x metagenomes

genus_to_type = {
    # HMA sponges
    'Agelas':           'HMA',
    'Aiolochroia':      'HMA',
    'Aplysina':         'HMA',
    'Chondrilla':       'HMA',
    'Coscinoderma':     'HMA',
    'Geodia':           'HMA',
    'Ircinia':          'HMA',
    'Petrosia':         'HMA',
    'Pseudoceratina':   'HMA',
    'Rhabdastrella':    'HMA',
    'Sarcotragus':      'HMA',
    'Smenospongia':     'HMA',
    'Theonella':        'HMA',
    'Thoosa':           'HMA',
    'Verongula':        'HMA',
    'Rhopaloeides':     'HMA',
    'Xestospongia':     'HMA',

    # LMA sponges
    'Amphimedon':       'LMA',
    'Axinella':         'LMA',
    'Baikalospongia':   'LMA',
    'Cinachyrella':     'LMA',
    'Clathria':         'LMA',
    'Cliona':           'LMA',
    'Crella':           'LMA',
    'Cymbastela':       'LMA',
    'Dysidea':          'LMA',
    'Ephydatia':        'LMA',
    'Eunapius':         'LMA',
    'Halichondria':     'LMA',
    'Haliclona':        'LMA',
    'Hymedesmia':       'LMA',
    'Ianthella':        'LMA',
    'Isodictya':        'LMA',
    'Lamellodysidea':   'LMA',
    'Leucetta':         'LMA',
    'Mycale':           'LMA',
    'Myxilla':          'LMA',
    'Niphates':         'LMA',
    'Phyllospongia':    'LMA',
    'Scopalina':        'LMA',
    'Spheciospongia':   'LMA',
    'Spongilla':        'LMA',
    'Stylissa':         'LMA',
    'Tedaniidae':       'LMA',
    'Pericharax':       'LMA',
    'Lophophysema':     'LMA',
    'Manihinea':        'LMA',
    'Haplosclerida':    'LMA',

    # Unknown status
     
    'Acarnus':          'N.D.',
    'Not_Defined':      'N.D.',
}

#Add HMA/LMA info to host data
host_df["Type"] = host_df["biome_genus"].map(genus_to_type)

# Filter out n.d. (unmapped)
host_df_filtered = host_df.dropna(subset=["Type"])

# Create lists of metagenomes by type
hma_mg = host_df_filtered.query("Type == 'HMA'")["Run"].tolist()
lma_mg = host_df_filtered.query("Type == 'LMA'")["Run"].tolist()

# Get only metagenomes present in the RPKM matrix
hma_mg = [m for m in hma_mg if m in rpkm_matrix.columns]
lma_mg = [m for m in lma_mg if m in rpkm_matrix.columns]

# Step 7: Define "presence"
rpkm_binary = rpkm_matrix >= 1  # True/False DataFrame

# Step 8: Find plasmids present in all HMA, absent in all LMA
hma_only = rpkm_binary[hma_mg].all(axis=1) & (~rpkm_binary[lma_mg].any(axis=1))

# Step 9: Find plasmids present in all LMA, absent in all HMA
lma_only = rpkm_binary[lma_mg].all(axis=1) & (~rpkm_binary[hma_mg].any(axis=1))

# Step 10: Extract plasmid names
plasmids_in_hma_only = rpkm_binary.index[hma_only].tolist()
plasmids_in_lma_only = rpkm_binary.index[lma_only].tolist()

# Optional: Save results
pd.Series(plasmids_in_hma_only).to_csv("C:/Users/hayat/Downloads/R_files/data/plasmids_in_hma_only.csv", index=False)
pd.Series(plasmids_in_lma_only).to_csv("C:/Users/hayat/Downloads/R_files/data/plasmids_in_lma_only.csv", index=False)
