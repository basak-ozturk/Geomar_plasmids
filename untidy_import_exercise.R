pacman::p_load(conflicted,tidyverse, wrappedtools, 
               readxl, readODS, foreign, haven, here)

data1 <- read_excel("../R_files/data/UntidyImportChallenge.xlsx", 
                     sheet = "Sheet1", range = "B4:E11", col_names= T)

data2 <- read_excel("../R_files/data/UntidyImportChallenge.xlsx", 
                    sheet = "Sheet1", range = "H4:K12", col_names= T)

data3 <- read_excel("../R_files/data/UntidyImportChallenge.xlsx", 
                    sheet = "Sheet1", range = "N4:Q7", col_names= T)

data_combined <- bind_rows(data1, data2, data3, .id = "origin") |>distinct() 

data_arranged <- data_combined[order(data_combined$`Meas./Treatm.`), decreasing = TRUE]

print(data_arranged)

long_data <- pivot_longer(
  data=data_arranged,
  cols=contains("h"),
  names_to = "All times", 
  values_to = "Weight"
)
long_data

wide_data <- pivot_wider(
  data=long_data,
  names_from = "All times",
  values_from = "Weight")

wide_data
