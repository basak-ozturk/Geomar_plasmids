library(VennDiagram)

# Define input files
file1 <- "C:/Users/hayat/Downloads/R_files/data/plasmids_with_plasmid_hallmark.txt"
file2 <- "C:/Users/hayat/Downloads/R_files/data/oriT_alloriT_blast_results_names.txt"
file3 <- "C:/Users/hayat/Downloads/R_files/data/unique_CONJscan_hit_ids.txt"

# Read content into sets
set1 <- readLines(file1)
set2 <- readLines(file2)
set3 <- readLines(file3)

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = list("Plasmid Hallmarks" = set1, "OriT" = set2, "CONJscan" = set3),
  filename = NULL,  # Store plot in variable
  fill = c("red", "green", "blue"),
  alpha = 0.5,
  cat.cex = 1.5,
  cex = 2,
  fontface = "bold",
  cat.fontface = "bold",
  main = "Venn Diagram of Sequences Identified",
  main.cex = 2
)

# Save the plot if it was created successfully
if (!is.null(venn.plot)) {
  png("C:/Users/hayat/Downloads/R_files/graphs/venn_diagram.png", width = 800, height = 800)
  grid.draw(venn.plot)
  dev.off()
} else {
  stop("Venn diagram creation failed!")
}
