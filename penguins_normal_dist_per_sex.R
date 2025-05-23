pacman::p_load(conflicted,
               tidyverse,
               wrappedtools,
               palmerpenguins,
               conflicted, 
               here,
               ggbeeswarm,
               hrbrthemes)


conflicts_prefer(dplyr::filter)
data(penguins)
head(penguins)

penguins <- na.omit(penguins) # Remove NA values


measurements <- c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g")

# initialize lists to store p-values for each sex
ks_p_values <- list(male = numeric(length(measurements)), female = numeric(length(measurements)))
shapiro_p_values <- list(male = numeric(length(measurements)), female = numeric(length(measurements)))

# Loop through each measurement
for (m in measurements) {
#  print(which(measurements == m))
  
  # split the dataset by sex
  data_by_sex <- split(penguins[[m]], penguins$sex)  
  #print(data_by_sex)
  for (sex in names(data_by_sex)) {
    datatotest <- data_by_sex[[sex]]  # Extract data for the current sex
    
    # Kolmogorov-Smirnov test
    ks_p_value <- ks.test(datatotest, "pnorm", mean(datatotest), sd(datatotest))$p.value
    ks_p_values[[sex]][which(measurements == m)] <- ks_p_value   
    
    # Shapiro-Wilk test
    shapiro_p_value <- shapiro.test(datatotest)$p.value
    shapiro_p_values[[sex]][which(measurements == m)] <- shapiro_p_value  
  }
}

print(ks_p_values)

print(shapiro_p_values)