pacman::p_load(conflicted,
               tidyverse,
               wrappedtools,
               palmerpenguins,
               ggbeeswarm,
               hrbrthemes)


conflicts_prefer(dplyr::filter)
data(penguins)
head(penguins)

penguins_filtered <- penguins %>% 
  filter(!is.na(sex))

# count sex within species

ggplot(penguins, aes(x = species, fill = sex)) +
  geom_bar() +
  labs(title = "Number of Penguins in Each Species per Sex")


# boxplot+beeswarm weight vs. species

ggplot(penguins, aes(x = species, y = body_mass_g)) +
  geom_boxplot(coef = 3) +
  ggbeeswarm::geom_beeswarm(cex = 2, size = 1, alpha = .25)

# boxplot+beeswarm weight vs. species AND sex

ggplot(penguins, aes(x = species, y = body_mass_g, fill = sex)) +
  geom_boxplot(coef = 3, position = position_dodge(0.8)) +
  ggbeeswarm::geom_beeswarm(cex = 2, size = 1, alpha = 0.25, dodge.width = 0.8)

# scatterplot flipper length vs. body mass

ggplot(data = penguins, aes(x = flipper_length_mm, y = body_mass_g, color = species)) +
  geom_point()

# group by species (and sex?)

ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g)) +
  geom_point(alpha = 0.5) +
  facet_wrap( ~ species + sex, scales = "free") 



ggplot(data=penguins_filtered,aes(flipper_length_mm, body_mass_g, color=sex, shape=species))+

  geom_smooth(linewidth = 1, color="red") +
  geom_smooth(method= "lm", linewidth=1, color="black") +
  geom_point(size=1 )+
  scale_color_manual(
    values = c("male" = "aquamarine2", "female" = "chocolate")
  ) +
  labs(
    title = "Scatterplot of Flipper Length vs Body Mass by Sex",
    x = "Flipper Length (mm)",
    y = "Body Mass (g)",
    color = "Sex"
  )
# now use bill depth instead of flipper length,
# facet by sex and species and compare linear/non-linear regression
ggplot(data = penguins_filtered, aes(bill_depth_mm, body_mass_g, color = sex)) +
  geom_point(size = 1) +
  geom_smooth(linewidth = 1, color = "red", aes(group = 1)) + 
  geom_smooth(method = "lm", linewidth = 1, aes(linetype = species), color = "black") + 
  facet_wrap(~ species + sex, scales = "free") +
  scale_color_manual(
    values = c("male" = "aquamarine2", "female" = "chocolate")
  ) +
  labs(
    title = "Scatterplot of Bill Depth vs Body Mass by Sex",
    x = "Bill Depth (mm)",
    y = "Body Mass (g)",
    color = "Sex"
  )

ggplot(data = penguins_filtered, aes(bill_depth_mm, body_mass_g, color = sex)) +
  geom_point(size = 1) +
  geom_smooth(linewidth = 1, color="red") +                                 
  geom_smooth(method = "lm", linewidth = 1, aes(linetype = species), color = "black") + 
  facet_grid(species ~ sex, scales = "free", margins = TRUE) +
  scale_color_manual(
    values = c("male" = "aquamarine2", "female" = "chocolate")
  ) +
  labs(
    title = "Scatterplot of Bill Depth vs Body Mass by Species and Sex (with Margins)",
    x = "Bill Depth (mm)",
    y = "Body Mass (g)",
    color = "Sex"
  )