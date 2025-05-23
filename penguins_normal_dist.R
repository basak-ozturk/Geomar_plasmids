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

penguins <- na.omit(penguins) # remove NA values

measurements <- c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g")

ks_p_values <- numeric(length(measurements))
shapiro_p_values <- numeric(length(measurements))


# loop through each measurement
for (m in measurements) {
  
  datatotest <- penguins[[m]] # Extract data to be tested
  
    counter=1 #index for inserting p-values into a vector
    
  # Kolmogorov-Smirnov test
  #kstest <- ks.test(datatotest, "pnorm", mean(datatotest), sd(datatotest))
  ks_p_value <- ks.test(datatotest, "pnorm", mean(datatotest), sd(datatotest))$p.value
  ks_p_values[counter]<-ks_p_value
  #print(cat("Kolmogorov-Smirnov p-value for", m, ":", ks_p_value, "\n"))
  
  # Shapiro-Wilk test
  shapiro_p_value <- shapiro.test(datatotest)$p.value
  #print(cat("Shapiro-Wilk p-value for", m, ":", shapiro_p_value, "\n\n"))
  shapiro_p_values[counter]<-shapiro_p_value
  
  counter<-counter+1
}

print(ks_p_values)
print(shapiro_p_values)

ggplot(penguins,aes(x = `bill_length_mm`))+
  geom_density(fill="pink")
ggplot(penguins,aes(x = `bill_length_mm`,fill=sex))+
  geom_density(alpha=.4)
ggplot(penguins,aes(x = `bill_length_mm`,fill=species))+
  geom_density(alpha=.4)

ggplot(penguins,aes(x = `bill_depth_mm`))+
  geom_density(fill="pink")
ggplot(penguins,aes(x = `bill_depth_mm`,fill=sex))+
  geom_density(alpha=.4)
ggplot(penguins,aes(x = `bill_depth_mm`,fill=species))+
  geom_density(alpha=.4)

ggplot(penguins,aes(x = `flipper_length_mm`))+
  geom_density(fill="pink")
ggplot(penguins,aes(x = `flipper_length_mm`,fill=sex))+
  geom_density(alpha=.4)
ggplot(penguins,aes(x = `flipper_length_mm`,fill=species))+
  geom_density(alpha=.4)

ggplot(penguins,aes(x = `body_mass_g`))+
  geom_density(fill="pink")
ggplot(penguins,aes(x = `body_mass_g`,fill=sex))+
  geom_density(alpha=.4)
ggplot(penguins,aes(x = `body_mass_g`,fill=species))+
  geom_density(alpha=.4)

# density plot with facet_grid
ggplot(penguins, aes(x = bill_length_mm, fill = species)) +
  geom_density(alpha = 0.4) +
  labs(title = "Density Plot of Bill Length by Sex and Species") +
  theme_minimal() +
  facet_grid(sex ~ species)

ggplot(penguins, aes(x = bill_depth_mm, fill = species)) +
  geom_density(alpha = 0.4) +
  labs(title = "Density Plot of Bill Depth by Sex and Species") +
  theme_minimal() +
  facet_grid(sex ~ species)

ggplot(penguins, aes(x = flipper_length_mm, fill = species)) +
  geom_density(alpha = 0.4) +
  labs(title = "Density Plot of Flipper Length by Sex and Species") +
  theme_minimal() +
  facet_grid(sex ~ species)

ggplot(penguins, aes(x = body_mass_g, fill = species)) +
  geom_density(alpha = 0.4) +
  labs(title = "Density Plot of Body Mass by Sex and Species") +
  theme_minimal() +
  facet_grid(sex ~ species)
