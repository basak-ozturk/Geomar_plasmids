# Load required libraries
library(tidyverse)
library(scales)  # For log scale and automatic breaks

# Read and process the file
data <- read.table("C:/Users/hayat/Downloads/R_files/data/KEGG_Pathway_per_plasmid_with_names.txt", header = FALSE, sep = "\t", col.names = c("Plasmid", "Pathways"))

data_expanded <- data |>
  separate_rows(Pathways, sep = ",") |>
  count(Pathways, name = "Count") |>
  filter(Count >= 100) |>  # Keep only pathways in ≥100 plasmids
  arrange(desc(Count))

ggplot(data_expanded, aes(x = log10(Count), y = reorder(Pathways, Count), size = Count)) +
  geom_point(alpha = 0.7, color = "lightblue") +
  scale_x_continuous(
    breaks = seq(1.5, 4, 0.5),
    labels = function(x) round(x, 1)
  ) +
  scale_size_continuous(trans = "log10", range = c(2, 10)) +  
  scale_y_discrete(expand = expansion(mult = 0.1)) +
  labs(
    x = "log10(Number of Plasmids)",
    y = "Pathway", 
    title = "Pathway Enrichment (≥100 Plasmids)", 
    size = "Plasmid Count"
  ) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),  
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )


# Save correctly formatted plot
ggsave("C:/Users/hayat/Downloads/R_files/graphs/KEGG_Pathway_BubblePlot_with_names.png", width = 10, height = 8, dpi = 300, bg = "white")
ggsave("C:/Users/hayat/Downloads/R_files/graphs/KEGG_Pathway_BubblePlot_with_names.svg", width = 10, height = 8, bg = "white")
