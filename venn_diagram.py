# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 07:45:23 2025

@author: hayat
"""

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

# === 2-circle Venn (IMG/PR vs Sponge Plasmids) ===
only_imgpr = 699973
only_sponge = 6823
overlap = 122

fig, ax = plt.subplots()
v = venn2(subsets=(only_imgpr, only_sponge, overlap), set_labels=('IMG/PR', 'Sponge Plasmids'), set_colors=('skyblue', 'lightgreen'))
v.get_label_by_id('10').set_text(str(only_imgpr))
v.get_label_by_id('01').set_text(str(only_sponge))
v.get_label_by_id('11').set_text(str(overlap))
ax.set_title("IMG/PR vs Sponge Plasmids (Nucleotide)")
plt.show()
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/venn2_imgpr_sponge_not_scaled.svg")

# === 3-circle Venn (Sponge vs IMG/PR vs PLSDB) ===
subsets = (
    84516,      # only Sponge
    4835662,    # only IMG/PR
    3100,       # Sponge ∩ IMG/PR only
    701321,     # only PLSDB
    53,         # Sponge ∩ PLSDB only
    573374,     # IMG/PR ∩ PLSDB only
    742         # Sponge ∩ IMG/PR ∩ PLSDB
)

fig3, ax3 = plt.subplots()
v3 = venn3(subsets=subsets, set_labels=('Sponge Plasmids', 'IMG/PR', 'PLSDB'), set_colors=('lightcoral', 'skyblue', 'gold'))
for label, count in zip(['100', '010', '110', '001', '101', '011', '111'], subsets):
    if v3.get_label_by_id(label):
        v3.get_label_by_id(label).set_text(str(count))
ax3.set_title("Sponge Plasmids vs IMG/PR vs PLSDB (Protein)")
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/venn3_sponge_imgpr_plsdb_not_scaled.svg")
plt.show()
#print("Saved: venn2_imgpr_sponge_not_scaled.svg and venn3_sponge_imgpr_plsdb_not_scaled.svg")
