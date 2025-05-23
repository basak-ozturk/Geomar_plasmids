library(tidyr)
library(dplyr)
library(readxl)
library(ggplot2)
library(paletteer)


# Specify the file path
file_path <- "C:/Users/hayat/OneDrive/Desktop/reddit_FFN_flair_counts.xlsx"

# Import the first sheet
flair_counts <- read_excel(file_path)

colnames(flair_counts)

library(tidyr)
library(dplyr)

# Rename the first column to 'flair'
colnames(flair_counts)[1] <- "flair"

# Reshape the data from wide to long format
long_data <- flair_counts %>%
  pivot_longer(
    cols = -flair,  # Select all columns except 'flair'
    names_to = "date",  # Name for the new date column
    values_to = "count"  # Name for the new count column
  )

# Convert the date column to Date type
long_data <- long_data %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d"))

# View the cleaned and reshaped data
head(long_data)

str(long_data)




ggplot(long_data, aes(x = date, y = count, color = flair, group = flair)) +
geom_line(size = 1) +
  scale_color_paletteer_d("ggthemes::Tableau_20") +  # Use the d3 palette
  scale_y_continuous(breaks = seq(0, max(long_data$count), by = 20)) +  # Adjust the y-axis ticks
    labs(
    title = "Trends in Post Counts by Flair Over Time",
    x = "Date",
    y = "Post Count",
    color = "Flair"
  ) +
  theme_minimal()


# Calculate the total posts for each date
total_data <- long_data |>
  group_by(date) |>
  summarise(count = sum(count)) |>
  mutate(flair = "Total")  # Add a new "flair" category for totals

# Combine the original data with the totals
long_data_with_total <- bind_rows(long_data, total_data)

# Check the updated data
head(long_data_with_total)


ggplot(long_data_with_total, aes(x = date, y = count, color = flair, group = flair)) +
  geom_line(size = 1) +
  scale_color_paletteer_d("ggthemes::Tableau_20") +  # Use the d3 palette
  #scale_y_continuous(breaks = seq(0, max(long_data$count), by = 20)) +  # Adjust the y-axis ticks
  labs(
    title = "Trends in Post Counts by Flair Over Time",
    x = "Date",
    y = "Post Count",
    color = "Flair"
  ) +
  theme_minimal()





