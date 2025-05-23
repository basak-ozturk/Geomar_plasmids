library(ggplot2)
library(flextable)
# 
# 
# file.exists("../data/all_sponge_plasmids.tsv")
# # Read the TSV file
# df <- read.delim("C:/Users/hayat/Downloads/R_files/data/all_sponge_plasmids.tsv", header = TRUE, sep = "\t", na.strings = "NA")
# 
# 
# # Count total genes
# total_genes <- nrow(df)
# 
# # Count genes with NA in the last and penultimate columns
# na_last_col <- sum(is.na(df[[ncol(df)]]))
# na_penultimate_col <- sum(is.na(df[[ncol(df)-1]]))
# 
# # Create a dataframe for plotting
# plot_data <- data.frame(
#   Category = c("Total Genes", "NA in Annotation Description", "NA in Annotation Accession"),
#   Count = c(total_genes, na_last_col, na_penultimate_col)
# )
# 
# # Create the bar plot-------------
# ggplot(plot_data, aes(x = Category, y = Count, fill = Category)) +
#   geom_bar(stat = "identity") +
#   theme_minimal() +
#   labs(title = "Gene Counts and NA Annotations for All Sponge Plasmids",
#        x = "",
#        y = "Number of Genes")

# Count total genes
total_genes <- nrow(df)

# Count genes with NA in the last column (Annotation Description)
na_last_col <- sum(is.na(df[[ncol(df)]]))

# Compute genes with annotation
genes_with_annotation <- total_genes - na_last_col

# Calculate percentage
plot_data$Percentage <- round((plot_data$Count / sum(plot_data$Count)) * 100, 1)

# Create labels with percentages
plot_data$Label <- paste0(plot_data$Category, "\n", plot_data$Percentage, "%")

png("C:/Users/hayat/Downloads/R_files/graphs/pie_chart.png", width = 800, height = 800)
# Create the pie chart
ggplot(plot_data, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Proportion of Genes with and without Annotation") +
  scale_fill_manual(values = c("pink", "lightblue")) +  # Adjust colors if needed
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 5)

dev.off()

annotation_accessions <- df[[ncol(df) - 1]]
# Remove NA values
annotation_accessions <- annotation_accessions[!is.na(annotation_accessions)]

# Count frequency of each unique value
freq_table <- table(annotation_accessions)

# Convert to a sorted data frame
freq_df <- as.data.frame(freq_table)
freq_df <- freq_df[order(-freq_df$Freq), ]  # Sort by descending frequency

# Print the frequency table
print(freq_df)

library(ggplot2)

# Filter frequencies greater than or equal to 50
filtered_freq_df <- freq_df[freq_df$Freq >= 50, ]

# Rename columns for clarity
colnames(filtered_freq_df) <- c("Function", "Frequency")

# Create a histogram (bar plot)--------
ggplot(filtered_freq_df, aes(x = reorder(Function, -Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Functions with Frequency ≥ 50",
       x = "Function",
       y = "Frequency") +
  theme(
    axis.text.x = element_text(size = 3, angle = 45, hjust = 1),  # X-axis labels
    axis.text.y = element_text(size = 10),  # Y-axis labels
    axis.title.x = element_text(size = 10),  # X-axis title
    axis.title.y = element_text(size = 10),  # Y-axis title
    plot.title = element_text(size = 12, face = "bold")  # Title size & bold
  )

# Create an empty dataframe to store the results
annotation_table <- data.frame(Function = character(0), Frequency = integer(0), Annotation_Description = character(0))

# Loop through each function in filtered_freq_df
for(i in 1:nrow(filtered_freq_df)) {
  # Get the current function
  current_function <- filtered_freq_df$Function[i]

  # Subset the rows in df that match the current function in the penultimate column
  matching_rows <- df[df[[ncol(df)-1]] == current_function, ]

  # Get the first non-NA annotation description
  non_na_description <- matching_rows[[ncol(df)]][!is.na(matching_rows[[ncol(df)]])][1]

  # Add the result to the annotation_table if a non-NA description exists
  if(!is.na(non_na_description)) {
    annotation_table <- rbind(annotation_table, data.frame(
      Function = current_function,
      Frequency = filtered_freq_df$Frequency[i],
      Annotation_Description = non_na_description
    ))
  }
}

# Print the final table
print(annotation_table)


# Convert the data frame to a flextable
annotation_flextable <- flextable(annotation_table)

# Customize the flextable (optional)
annotation_flextable <- annotation_flextable |>
  compose(j = "Function", value = as_paragraph(as_chunk(Function))) |>
  compose(j = "Frequency", value = as_paragraph(as_chunk(Frequency))) |>
  compose(j = "Annotation_Description", value = as_paragraph(as_chunk(Annotation_Description))) |>
  set_table_properties(layout = "autofit") |>
  theme_vanilla()

# Print the flextable
annotation_flextable

# Filter frequencies smaller than 50
filtered_freq_df_small <- freq_df[freq_df$Freq < 50, ]

# Rename columns for clarity
colnames(filtered_freq_df_small) <- c("Function", "Frequency")

# Create an empty dataframe to store the results
annotation_table_small <- data.frame(Function = character(0), Frequency = integer(0), Annotation_Description = character(0))

# Loop through each function in filtered_freq_df_small
for(i in 1:nrow(filtered_freq_df_small)) {
  # Get the current function
  current_function <- filtered_freq_df_small$Function[i]

  # Subset the rows in df that match the current function in the penultimate column
  matching_rows <- df[df[[ncol(df)-1]] == current_function, ]

  # Get the first non-NA annotation description
  non_na_description <- matching_rows[[ncol(df)]][!is.na(matching_rows[[ncol(df)]])][1]

  # Add the result to the annotation_table if a non-NA description exists
  if(!is.na(non_na_description)) {
    annotation_table_small <- rbind(annotation_table_small, data.frame(
      Function = current_function,
      Frequency = filtered_freq_df_small$Frequency[i],  # Corrected here
      Annotation_Description = non_na_description
    ))
  }
}

print(annotation_table_small)


# Convert the data frame to a flextable
annotation_flextable_small <- flextable(annotation_table_small)

annotation_flextable_small <- annotation_flextable_small |>
  compose(j = "Function", value = as_paragraph(as_chunk(Function))) |>
  compose(j = "Frequency", value = as_paragraph(as_chunk(Frequency))) |>
  compose(j = "Annotation_Description", value = as_paragraph(as_chunk(Annotation_Description))) |>
  set_table_properties(layout = "autofit") |>
  theme_vanilla()

# Print the flextable
annotation_flextable_small

# Save the large annotation table as a TSV file
write.table(annotation_table, file = "C:/Users/hayat/Downloads/R_files/data/annotation_table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Save the small annotation table as a TSV file
write.table(annotation_table_small, file = "C:/Users/hayat/Downloads/R_files/data/annotation_table_small.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

