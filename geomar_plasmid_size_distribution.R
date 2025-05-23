pacman::p_load(tidyverse, wrappedtools
               )
data <- read.csv("data/sequence_sizes.csv", sep="\t")  # Use `\t` for tab-separated files


data$Sequence_Length <- as.numeric(gsub(",", "", data$Sequence_Length))

# Convert Sequence_Length to numeric
data$Sequence_Length <- as.numeric(gsub(",", "", data$Sequence_Length))

# Debugging: Check for NA values and print summary statistics
if (any(is.na(data$Sequence_Length))) {
  print("There are NA values in Sequence_Length. Check your data!")
}
print(summary(data$Sequence_Length))

# Add a small constant to avoid issues with log(0)
data <- data %>% mutate(Log_Length = log10(Sequence_Length + 1))

# Plot the histogram with a log-transformed x-axis
ggplot(data, aes(x = Log_Length)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Log-Transformed Sequence Lengths",
       x = "Log10(Sequence Length + 1)",
       y = "Frequency") +
  theme_minimal() +
  scale_x_continuous(labels = scales::math_format(10^.x)) +  # Logarithmic scale labels
  scale_y_continuous(labels = scales::comma)                # Format y-axis labels

# Filter for sequence lengths up to 10,000 bp
filtered_data <- data |> filter(Sequence_Length <= 10000)

# Plot histogram for filtered data
ggplot(filtered_data, aes(x = Sequence_Length)) +
  geom_histogram(binwidth = 500, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Sequence Lengths (≤ 10 kbp)",
       x = "Sequence Length (bp)",
       y = "Frequency") +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(labels = scales::comma)
