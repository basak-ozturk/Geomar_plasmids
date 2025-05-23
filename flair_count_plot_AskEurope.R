library(tidyr)
library(dplyr)
library(readxl)
library(ggplot2)
library(paletteer)


file_path <- "C:/Users/hayat/OneDrive/Desktop/reddit_AskEurope_flair_counts.xlsx"


flair_counts <- read_excel(file_path)

colnames(flair_counts)


# Rename the first column to 'flair'
colnames(flair_counts)[1] <- "flair"

# Reshape the data from wide to long format in order to feed it into ggplot. There should be columns flair, date and count. We don't touch the "flair" column
long_data <- flair_counts |>
  pivot_longer(
    cols = -flair,  # Select all columns except 'flair'
    names_to = "date", 
    values_to = "count"  
  )

# Convert the date column to date type so that it can be understood as such
long_data <- long_data |>
  mutate(date = as.Date(date, format = "%Y-%m-%d"))

# Calculate the total posts for each date. First we group by date and then count all the posts across all flairs for that given date.
total_data <- long_data |>
  group_by(date) |>
  summarise(count = sum(count)) |>
  mutate(flair = "Total")  # Add a new "flair" category for totals


long_data_with_total <- bind_rows(long_data, total_data) # Combine the new category with the previous table

#head(long_data_with_total)

ggplot(long_data, aes(x = date, y = count, color = flair, group = flair)) +
  geom_line(linewidth = 1) +
  scale_color_paletteer_d("ggthemes::Tableau_20") + 
  scale_y_continuous(breaks = seq(0, max(long_data$count), by = 20)) +  
  labs(
    title = "Trends in Post Counts by Flair Over Time",
    x = "Date",
    y = "Post Count",
    color = "Flair"
  ) +
  theme_minimal()

ggplot(long_data_with_total, aes(x = date, y = count, color = flair, group = flair)) +
  geom_line(linewidth = 1) +
  scale_color_paletteer_d("ggthemes::Tableau_20") +  
  scale_y_continuous(breaks = seq(0, max(long_data$count), by = 20)) + 
  labs(
    title = "Trends in Post Counts by Flair Over Time",
    x = "Date",
    y = "Post Count",
    color = "Flair"
  ) +
  theme_minimal()

ggplot(long_data_with_total, aes(x = date, y = count, color = flair, group = flair)) +
  geom_line(linewidth = 1) +
  scale_color_paletteer_d("ggthemes::Tableau_20") +  
  scale_y_continuous() + 
  labs(
    title = "Trends in Post Counts by Flair Over Time",
    x = "Date",
    y = "Post Count",
    color = "Flair"
  ) +
  facet_wrap(~flair, scales = "free_y") +  # scale freely per facet
  theme_minimal()
