# Load necessary libraries
pacman::p_load(conflicted, tidyverse, wrappedtools, readxl, writexl)

# Set preference for dplyr functions
conflicts_prefer(dplyr::filter)

# Load data
file_path <- "C:/Users/hayat/Downloads/R_files/data/200228_NJ_Ingrid_Protein_Report_n6.xlsx"
proteome_data <- read_excel(file_path)

# Clean column names (replace spaces with underscores)
colnames(proteome_data) <- gsub(" ", "_", colnames(proteome_data)) 

# Define RF and RS columns
RF_columns <- c("RF1_Area", "RF2_Area", "RF3_Area")
RS_columns <- c("RS1_Area", "RS2_Area", "RS3_Area")

# Step 1: Log Transformation (log2)
proteome_data <- proteome_data |>
  mutate(across(c(all_of(RF_columns), all_of(RS_columns)), ~ log2(. + 10000), .names = "log2_{.col}"))

# Step 2: Normalize samples (by column median and multiplying by overall median)
normalize_samples <- function(data, columns) {
  # Calculate the median of each column (log2-transformed columns)
  column_medians <- data |>
    summarize(across(all_of(columns), ~ median(. , na.rm = TRUE), .names = "median_{.col}"))
  
  # Calculate the median of all column medians
  median_of_medians <- median(unlist(column_medians), na.rm = TRUE)
  
  # Normalize the columns by dividing by column median and multiplying by the overall median
  normalized_data <- data |>
    mutate(across(all_of(columns), ~ . / column_medians[[paste0("median_", cur_column())]] * median_of_medians, .names = "normalized_{.col}"))
  
  return(normalized_data)
}

# Normalize RF and RS columns
proteome_data <- normalize_samples(proteome_data, c(paste0("log2_", RF_columns)))
proteome_data <- normalize_samples(proteome_data, c(paste0("log2_", RS_columns)))

# Step 3: Remove proteins missing in more than 1 replicate in either RF or RS group
proteome_data_filtered <- proteome_data |>
  filter(
    rowSums(is.na(across(c("log2_RF1_Area", "log2_RF2_Area", "log2_RF3_Area")))) <= 1,  # No more than 1 missing in RF
    rowSums(is.na(across(c("log2_RS1_Area", "log2_RS2_Area", "log2_RS3_Area")))) <= 1   # No more than 1 missing in RS
  )

# Step 4: Calculate the average for each group (RF and RS)
proteome_data_filtered <- proteome_data_filtered |>
  mutate(
    RF_Avg = rowMeans(across(starts_with("normalized_log2_RF")), na.rm = TRUE),
    RS_Avg = rowMeans(across(starts_with("normalized_log2_RS")), na.rm = TRUE)
  )

# Step 5: Calculate the log2 fold change (RF vs RS)
proteome_data_filtered <- proteome_data_filtered |>
  mutate(
    log2_fold_change = RF_Avg - RS_Avg  # Subtract averages directly
  )

# Step 6: Create a separate table for rows that were excluded (rows that don't meet the criteria)
proteome_data_excluded <- proteome_data |>
  filter(
    rowSums(is.na(across(c("log2_RF1_Area", "log2_RF2_Area", "log2_RF3_Area")))) > 1 |  # More than 1 missing in RF
      rowSums(is.na(across(c("log2_RS1_Area", "log2_RS2_Area", "log2_RS3_Area")))) > 1   # More than 1 missing in RS
  )

# Step 7: Write both datasets to an Excel file with separate sheets
write_xlsx(
  list(
    "Filtered Data" = proteome_data_filtered,
    "Excluded Data" = proteome_data_excluded
  ),
  path = "proteome_data_filtered.xlsx"
)


