pacman::p_load(conflicted,
               tidyverse,
               wrappedtools,
               readxl,
               here,
               flextable)
conflicts_prefer(dplyr::filter)
datasampleweight <- read_excel(
  "../R_files/data/DOC-20230130-WA0000_.xlsx",
  sheet = "Sheet1",
  range = "A1:E33",
  col_names = T
)

# Convert column names to title case and correct typos
datasampleweight <- datasampleweight |>
  rename_with(.fn = str_to_title) |>
  rename_all(~ str_replace_all(.x, c("Weigth" = "Weight", "Alminium" = "Aluminum")))

colnames(datasampleweight)

# Calculate new column for weight of sample before drying
datasampleweight <- datasampleweight |>
  mutate(
    `Weight Of Sample Before Drying` = `Weight Of Aluminum Cup + Sample (Wt + S)` - `Weight Of Empty Aluminum(Wt)`
  )

# Calculate new column for weight of sample after drying
datasampleweight <- datasampleweight |>
  mutate(
    `Weight Of Sample After Drying` = `Weight Of Aluminium Cap + Sample  After Drying (Wt-Al +S+D)` - `Weight Of Empty Aluminum(Wt)`
  )

# Calculate new column for percent moisture
datasampleweight <- datasampleweight |>
  mutate(
    Percentage_Moisture = ((`Weight Of Sample Before Drying` - `Weight Of Sample After Drying`) / `Weight Of Sample Before Drying`) * 100
  )

# Replace problem names in the 'Code Of Cup' column
datasampleweight <- datasampleweight |>
  mutate(`Code Of Cup` = str_replace_all(`Code Of Cup`, "/", "_"))|>
  mutate(`Code Of Cup` = str_replace_all(`Code Of Cup`, "\\*", "Asterix")) |>
  mutate(`Code Of Cup` = str_replace_all(`Code Of Cup`, "\\#", "Hash")) |>
  mutate(`Code Of Cup` = replace_na(`Code Of Cup`, "none"))

measurements <- c("Weight Of Sample Before Drying", "Weight Of Sample After Drying", "Percentage_Moisture")


shapiro_p_values <- numeric(length(measurements))


# loop through each measurement
for (m in measurements) {
  
  datatotest <- datasampleweight[[m]] # Extract data to be tested 
  
  counter=1 #index for inserting p-values into a vector
  

  # Shapiro-Wilk test
  shapiro_p_value <- shapiro.test(datatotest)$p.value
  print(cat("Shapiro-Wilk p-value for", m, ":", shapiro_p_value, "\n\n"))
  shapiro_p_values[counter]<-shapiro_p_value
  
  counter<-counter+1
}

print(shapiro_p_values)

ggplot(datasampleweight,aes(x = `Weight Of Sample Before Drying`))+
  geom_density(fill="pink")

ggplot(datasampleweight,aes(x = `Weight Of Sample After Drying`))+
  geom_density(fill="pink")
ggplot(datasampleweight,aes(x = Percentage_Moisture))+
  geom_density(fill="pink")

data_summary <- datasampleweight |>
  group_by(`Code Of Cup`) |>
  summarize(
    count = n(),
    median_weight_before = median(`Weight Of Sample Before Drying`),
    median_weight_after = median(`Weight Of Sample After Drying`),
    IQR_before = IQR(`Weight Of Sample Before Drying`),
    IQR_after = IQR(`Weight Of Sample After Drying`),
    min_before = min(`Weight Of Sample Before Drying`),
    max_before = max(`Weight Of Sample Before Drying`),
    min_after = min(`Weight Of Sample After Drying`),
    max_after = max(`Weight Of Sample After Drying`)
  )
table_1 <- flextable(data_summary)

(table_1)

data_summary2 <- datasampleweight |>
  group_by(`Code Of Sample`) |>
  summarize(
    count = n(),
    median_weight_before = median(`Weight Of Sample Before Drying`),
    median_weight_after = median(`Weight Of Sample After Drying`),
    median_moisture = median(Percentage_Moisture),
    median_weight_after = median(`Weight Of Sample After Drying`),
    IQR_before = IQR(`Weight Of Sample Before Drying`),
    IQR_after = IQR(`Weight Of Sample After Drying`),
    IQR_moisture = IQR(Percentage_Moisture),
    min_before = min(`Weight Of Sample Before Drying`),
    max_before = max(`Weight Of Sample Before Drying`),
    min_after = min(`Weight Of Sample After Drying`),
    max_after = max(`Weight Of Sample After Drying`),
    min_moisture = min(Percentage_Moisture),
    max_moisture = max(Percentage_Moisture)
  )



ggplot(datasampleweight, aes(x = `Weight Of Sample Before Drying`, y = `Weight Of Sample After Drying`)) +
  geom_point() +
  labs(title = "Scatter plot of weight before vs after drying",
       x = "Weight Of Sample Before Drying",
       y = "Weight Of Sample After Drying") +
  theme_minimal()

# Extract rows with Code Of Sample = "A"
df_A <- datasampleweight |>
  filter(`Code Of Sample` == "A")

# Extract rows with Code Of Sample = "D"
df_D <- datasampleweight |>
  filter(`Code Of Sample` == "D")

ggplot(df_A,aes(x = `Weight Of Sample Before Drying`))+
  geom_density(fill="pink")

ggplot(df_D,aes(x = `Weight Of Sample Before Drying`))+
  geom_density(fill="pink")

ggplot(df_A,aes(x = `Weight Of Sample After Drying`))+
  geom_density(fill="pink")

ggplot(df_D,aes(x = `Weight Of Sample After Drying`))+
  geom_density(fill="pink")

ggplot(df_A,aes(x = Percentage_Moisture))+
  geom_density(fill="pink")

ggplot(df_D,aes(x = Percentage_Moisture))+
  geom_density(fill="pink")

shapiro_p_values_A <- numeric(length(measurements))


# loop through each measurement
for (m in measurements) {
  
  datatotest <- df_A[[m]] # Extract data to be tested 
  
  counter=1 #index for inserting p-values into a vector
  
  
  # Shapiro-Wilk test
  shapiro_p_value <- shapiro.test(datatotest)$p.value
  print(cat("Shapiro-Wilk p-value for", m, ":", shapiro_p_value, "\n\n"))
  shapiro_p_values[counter]<-shapiro_p_value
  
  counter<-counter+1
}



shapiro_p_values_D <- numeric(length(measurements))


# loop through each measurement
for (m in measurements) {
  
  datatotest <- df_D[[m]] # Extract data to be tested 
  
  counter=1 #index for inserting p-values into a vector
  
  
  # Shapiro-Wilk test
  shapiro_p_value <- shapiro.test(datatotest)$p.value
  print(cat("Shapiro-Wilk p-value for", m, ":", shapiro_p_value, "\n\n"))
  shapiro_p_values[counter]<-shapiro_p_value
  
  counter<-counter+1
}



