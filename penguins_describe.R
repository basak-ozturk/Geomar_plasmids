pacman::p_load(conflicted,tidyverse,wrappedtools,
               flextable, here, palmerpenguins)
set_flextable_defaults(font.size = 9, 
                       padding.bottom = 1, 
                       padding.top = 3,
                       padding.left = 3,
                       padding.right = 4
)
data(penguins)
penguins <- na.omit(penguins) # Remove NA values

ggplot(penguins,aes(bill_length_mm, sex))+
  geom_point()+
  geom_smooth(se=F)+
  geom_smooth(method="lm",color="red", 
              fill="gold", alpha=.15)

ggplot(penguins,aes(flipper_length_mm, sex))+
  geom_point()+
  geom_smooth(se=F)+
  geom_smooth(method="lm",color="red", 
              fill="gold", alpha=.15)

(mean_size_flipper <- mean(penguins$flipper_length_mm))
(sd_size_flipper <- sd(penguins$flipper_length_mm))
min(penguins$flipper_length_mm)
SEM(penguins$flipper_length_mm)

(mean_size_body <- mean(penguins$body_mass_g))
(sd_size_body <- sd(penguins$body_mass_g))
min(penguins$body_mass_g)
SEM(penguins$body_mass_g)

# meansd(penguins$body_mass_g, roundDig = 4,
#        range = TRUE,
#        add_n = TRUE)
# meansd(penguins$flipper_length_mm, roundDig = 4,
#        range = TRUE,
#        add_n = TRUE)
# meanse(penguins$body_mass_g, roundDig = 4)
# meanse(penguins$flipper_length_mm, roundDig = 4)

cat(
  "Summary of the normally distributed measurements in the penguins dataset:\n",
  
  "Body Mass (g):\n",
  "Mean (SD, range, n): ", meansd(penguins$body_mass_g, roundDig = 4, range = TRUE, add_n = TRUE), "\n",
  "Mean (SE): ", meanse(penguins$body_mass_g, roundDig = 4), "\n\n",
  
  "Flipper Length (mm):\n",
  "Mean (SD, range, n): ", meansd(penguins$flipper_length_mm, roundDig = 4, range = TRUE, add_n = TRUE), "\n",
  "Mean (SE): ", meanse(penguins$flipper_length_mm, roundDig = 4), "\n"
)


# median(penguins$bill_length_mm)
# quantile(penguins$bill_length_mm,probs = c(.25,.75))
# median_quart(penguins$bill_length_mm)
# 
# median(penguins$bill_depth_mm)
# quantile(penguins$bill_depth_mm,probs = c(.25,.75))
# median_quart(penguins$bill_depth_mm,range = T)

cat(
  "Summary of the ordinal measurements in the penguins dataset:\n\n",
  
  "Bill Length (mm):\n",
  "Median: ", median(penguins$bill_length_mm), "\n",
  "IQR (25th and 75th percentiles): ", quantile(penguins$bill_length_mm, probs = c(0.25, 0.75)), "\n",
  "Median and Quartiles: ", median_quart(penguins$bill_length_mm), "\n\n",
  
  "Bill Depth (mm):\n",
  "Median: ", median(penguins$bill_depth_mm), "\n",
  "IQR (25th and 75th percentiles): ", quantile(penguins$bill_depth_mm, probs = c(0.25, 0.75)), "\n",
  "Median, Quartiles, and Range: ", median_quart(penguins$bill_depth_mm, range = TRUE), "\n"
)

penguins |> group_by(species, sex) |>
  
  summarise(
    Count = n(),
    WeightSummary = paste(
      round(mean(body_mass_g), 2), " ± ", 
      round(sd(body_mass_g), 2)
    ),
    .groups = "drop"
  )|>
  flextable()




penguins_summary <- penguins |>
  group_by(species, sex) |>
  summarise(
    Count = n(),
    WeightSummary = paste(round(mean(body_mass_g), 2), " ± ", round(sd(body_mass_g), 2)),
    .groups = "drop"
  )


ft <- flextable(penguins_summary) |>
  set_header_labels(
    species = "Species",
    sex = "Sex",
    Count = "Number",
    WeightSummary = "Body Mass (g) Mean ± SD"
  ) |>
  theme_vanilla() |>
  align(align = "center", part = "all") |>
  autofit()


(ft)

summary_table_2 <- penguins |>
  
  group_by(species, sex) |>
  summarise(Count = n()) |>
  ungroup() |>
  group_by(species) |>
  mutate(Percentage = round(Count / sum(Count) * 100, 1)) |>
  ungroup()


ft2 <- flextable(summary_table_2) |>
  set_header_labels(
    species = "Species",
    sex = "Sex",
    Count = "Number",
    Percentage = "Percentage (%)"
  )|>
  theme_vanilla() |>
  align(align = "center", part = "all") |>
  autofit()


(ft2)