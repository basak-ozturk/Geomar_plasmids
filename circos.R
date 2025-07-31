library(circlize)
library(dplyr)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)
# Load node metadata (includes Host_Group)
nodes <- read.csv("C:/Users/hayat/Downloads/R_files/data/nodes_with_host_groups.csv")

# Load host mapping file
host_map <- read.csv("C:/Users/hayat/Downloads/R_files/data/top_widespread_hosts.csv")

# Load MCL cluster info
mcl_df <- read.csv("C:/Users/hayat/Downloads/R_files/data/cytoscape_edges_with_mcl_clusters.csv")

# Extract unique node-cluster relationships
mcl_nodes <- mcl_df |>
  select(name, X__mclCluster) |>
  rename(ID = name) |>
  rename(Cluster = X__mclCluster) |>
  distinct()

# Merge nodes + host + cluster info
merged_df <- nodes |>
  inner_join(host_map, by = "ID") |>
  left_join(mcl_nodes, by = "ID")

# OPTIONAL: Filter to one MCL cluster 
 #merged_df <- merged_df |> filter(Cluster == 2)

# Count number of plasmids per Host ↔ Host_Group
chord_data <- merged_df |>
  group_by(Host, Host_Group) |>
  summarise(Count = n(), .groups = "drop")

# Convert to wide matrix format
chord_matrix <- pivot_wider(
  chord_data,
  names_from = Host_Group,
  values_from = Count,
  values_fill = 0
)
colnames(chord_matrix) <- gsub("[–—−]", "-", colnames(chord_matrix))  # all dash variants to hyphen-minus

desired_order <- c("1", "2-3", "4-7", "8-14", "15+")
all(desired_order %in% colnames(chord_matrix))  # should be TRUE

# Reorder columns
chord_matrix <- chord_matrix[, desired_order]
chord_matrix <- as.matrix(chord_matrix)  

all_hosts <- sort(unique(merged_df$Host)) 

# Set row names and convert to matrix
rownames(chord_matrix) <- all_hosts
# Define color palettes
# Define fixed host color map from ALL unique hosts in the nodes dataset
set.seed(42)
host_colors <- setNames(viridis(length(all_hosts), option = "C"), all_hosts)


# Likewise for Host_Group
all_host_groups <- sort(unique(nodes$Host_Group))
set.seed(24)
host_group_colors <- setNames(brewer.pal(length(all_host_groups), "Set3"), all_host_groups)

# Combine into one color mapping
grid.col <- c(host_colors, host_group_colors)

# Create Host legend
host_legend <- Legend(
  title = "Host",
  legend_gp = gpar(fill = host_colors),
  labels = names(host_colors),
  ncol = 1,
  grid_height = unit(4, "mm"),
  grid_width = unit(4, "mm"),
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 8)
)
# Plot
circos.clear()
#pdf("C:/Users/hayat/Downloads/R_files/graphs/circos_plot_host_count.pdf")
svg("C:/Users/hayat/Downloads/R_files/graphs/circos_plot_host_count.svg")

chordDiagram(chord_matrix, grid.col = grid.col, transparency = 0.2)
dev.off()

svg("C:/Users/hayat/Downloads/R_files/graphs/legend_only.svg", width = 5, height = 10)
draw(host_legend, x = unit(0.5, "npc"), just = "center")
dev.off()

svg("C:/Users/hayat/Downloads/R_files/graphs/circos_plot_host_count_slanted.svg")

# Prepare the plot with space for labels on the outer track
chordDiagram(
  chord_matrix,
  grid.col = grid.col,
  transparency = 0.2,
  annotationTrack = "grid",   # only draw grid tracks, no default labels
  preAllocateTracks = 1       # allocate one track for labels
)

# Draw labels manually with rotation
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector_name = get.cell.meta.data("sector.index")
  xcenter = get.cell.meta.data("xcenter")
  ycenter = get.cell.meta.data("ylim")[1] + mm_y(1)  # slightly outside sector
  
  # Draw vertical text (90 degrees)
  circos.text(
    xcenter, ycenter,
    sector_name,
    #facing = "clockwise",  # or "inside", "outside"
    niceFacing = TRUE,
    adj = c(0, 0.5),
    facing = "inside",
    cex = 0.7,
    srt = 45               # rotation angle in degrees, try 45, 90, etc.
  )
}, bg.border = NA)

dev.off()

