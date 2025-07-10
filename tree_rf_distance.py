# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 10:09:00 2025

@author: hayat
"""

from ete3 import Tree

# Load trees
tree_ani = Tree("C:/Users/hayat/Downloads/R_files/data/plasmid_ani_tree.nwk")
tree_mash = Tree("C:/Users/hayat/Downloads/R_files/data/plasmid_mash_tree.nwk")
tree_pscope = Tree("C:/Users/hayat/Downloads/R_files/data/plasmidscope_comparative_tree_cleaned.nwk")  

# Find common tips across all trees
common_tips = set(tree_ani.get_leaf_names()) & set(tree_mash.get_leaf_names()) & set(tree_pscope.get_leaf_names())
print(f"Common tips: {len(common_tips)}")

# Prune trees to common tips
for tree in [tree_ani, tree_mash, tree_pscope]:
    tree.prune(common_tips)

# Save pruned trees
tree_ani.write(outfile="C:/Users/hayat/Downloads/R_files/data/plasmid_ani_tree_pruned.nwk")
tree_mash.write(outfile="C:/Users/hayat/Downloads/R_files/data/plasmid_mash_tree_pruned.nwk")
tree_pscope.write(outfile="C:/Users/hayat/Downloads/R_files/data/plasmidscope_tree_pruned.nwk")

rf_ani_mash, max_rf, *_ = tree_ani.robinson_foulds(tree_mash, unrooted_trees=True)
rf_ani_pscope, _, *_ = tree_ani.robinson_foulds(tree_pscope, unrooted_trees=True)
rf_mash_pscope, _, *_ = tree_mash.robinson_foulds(tree_pscope, unrooted_trees=True)


print(f"RF distance ANI vs Mash: {rf_ani_mash} / {max_rf} (normalized {rf_ani_mash/max_rf:.3f})")
print(f"RF distance ANI vs PlasmidScope: {rf_ani_pscope} / {max_rf} (normalized {rf_ani_pscope/max_rf:.3f})")
print(f"RF distance Mash vs PlasmidScope: {rf_mash_pscope} / {max_rf} (normalized {rf_mash_pscope/max_rf:.3f})")
