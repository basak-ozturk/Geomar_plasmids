pacman::p_load(conflicted,tidyverse,multcomp,wrappedtools, palmerpenguins, broom, car)
conflicts_prefer(dplyr::select)

data(penguins)
head(penguins)

penguins <- na.omit(penguins) # remove NA values


ggplot(penguins, aes(x = body_mass_g, y = flipper_length_mm, fill = sex)) +
  ggbeeswarm::geom_beeswarm(aes(shape = sex), alpha = .5, dodge.width = .15) +
  geom_smooth(method = "lm") +
  scale_x_continuous("Body Mass (g)")

ggplot(penguins, aes(x = body_mass_g, y = flipper_length_mm)) +
  ggbeeswarm::geom_beeswarm(alpha = .5, dodge.width = .15) +
  geom_smooth(method = "lm") +
  scale_x_continuous("Body Mass (g)")

ggplot(penguins, aes(x = body_mass_g, y = bill_depth_mm, fill = sex)) +
  ggbeeswarm::geom_beeswarm(aes(shape = sex), alpha = .5, dodge.width = .15) +
  geom_smooth(method = "lm") +
  scale_x_continuous("Body Mass (g)")

ggplot(penguins, aes(x = body_mass_g, y = bill_depth_mm, fill=species)) +
  ggbeeswarm::geom_beeswarm(alpha = .5, dodge.width = .15) +
  geom_smooth(method = "lm") +
  scale_x_continuous("Body mass(g)")

ggplot(penguins, aes(x = body_mass_g, y = bill_depth_mm, fill = species)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_beeswarm(alpha = .5) +
  scale_y_continuous()c


# Regression like lm
# flipper length by body mass
# bill depth by body mass

#regression flipper_length vs body mass

(regression_flipper_length<-lm(flipper_length_mm~body_mass_g,data=penguins))
summary(regression_flipper_length)
regression_flipper_length$coefficients

# computation of SSQs and p-values, flipper length

(anova_out_flipper<-anova(regression_flipper_length))
anova_out_flipper$`Pr(>F)`
tidy(anova_out_flipper)

#regression bill_depth vs body mass

(regression_bill_depth<-lm(bill_depth_mm~body_mass_g,data=penguins))
summary(regression_bill_depth)
regression_bill_depth$coefficients


# computation of SSQs and p-values, bill_depth
(anova_out_bill<-anova(regression_bill_depth))
anova_out_bill$`Pr(>F)`
tidy(anova_out_bill)

# ANOVA like lm
# body mass by species
# body mass by species and sex (interaction?)

ggplot(penguins, aes(x = species, y = body_mass_g, fill = sex, color = sex)) +
  ggbeeswarm::geom_beeswarm(aes(shape = sex), alpha = .25, dodge.width = .30) +
  geom_smooth(method = "lm") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "red",
               position = position_dodge(width = 0.30)) +  # Adding mean points
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", 
               position = position_dodge(width = 0.30))

# body mass by species  
(anova_mass_by_species<-lm(body_mass_g ~ species,data=penguins))
tidy(anova_mass_by_species)

(t <- anova(anova_mass_by_species))
t$`Pr(>F)`
tidy(t)

#post-hoc test

pairwise.t.test(penguins$body_mass_g, penguins$species, p.adjust.method = "bonferroni")

# body mass by species and sex (interaction?)

#with interaction
(lm_out_species_sex <- lm(body_mass_g~species*sex,
              data=penguins))

Anova(lm_out_species_sex,type = 3) |> 
  tidy() |> slice(1:4) |> 
  mutate(p.value=formatP(p.value,ndigits=3, mark=TRUE))

#without interaction

(lm_out_species_sex <- lm(body_mass_g~species+sex,
                          data=penguins))

Anova(lm_out_species_sex,type = 2) |> 
  tidy() |> slice(1:2) |> 
  mutate(p.value=formatP(p.value,ndigits=3, mark=TRUE))

# General  lm
# flipper length by body mass and species and sex (interaction?)

ggplot(penguins, aes(x = body_mass_g, y = flipper_length_mm, color = species, shape = sex)) +
  geom_point(alpha = 0.6) +  
  geom_smooth(method = "lm", aes(color = species), se = FALSE) +  #regression line by species
  labs(x = "Body Mass (g)", y = "Flipper Length (mm)", title = "Flipper Length by Body Mass, Species, and Sex") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("Adelie" = "#F8766D", "Chinstrap" = "#00BFC4", "Gentoo" = "#7CAE00"))  


#no interaction

(regression_flipper_length_2<-lm(flipper_length_mm~body_mass_g+species+sex,data=penguins))

Anova(regression_flipper_length_2,type = 2) |> 
  tidy() |> slice(1:2) |> 
  mutate(p.value=formatP(p.value,ndigits=3, mark=TRUE))


#interaction

(regression_flipper_length_3<-lm(flipper_length_mm~body_mass_g*species*sex,data=penguins))

Anova(regression_flipper_length_3,type = 3) |> 
  tidy() |> slice(1:4) |> 
  mutate(p.value=formatP(p.value,ndigits=3, mark=TRUE))

#post-hoc pairwise comparisons for species and sex
posthoc_species <- glht(regression_flipper_length_2, linfct = mcp(species = "Tukey"))
posthoc_sex <- glht(regression_flipper_length_2, linfct = mcp(sex = "Tukey"))


summary(posthoc_species)
summary(posthoc_sex)


# General  lm
# bill depth by body mass and species and sex (interaction?)

#no interaction

(regression_bill_depth_2<-lm(bill_depth_mm~body_mass_g+species+sex,data=penguins))
regression_bill_depth_2$coefficients

Anova(regression_bill_depth_2,type = 2) |> 
  tidy() |> slice(1:2) |> 
  mutate(p.value=formatP(p.value,ndigits=3, mark=TRUE))


#interaction

(regression_bill_depth_3<-lm(bill_depth_mm~body_mass_g*species*sex,data=penguins))
regression_bill_depth_3$coefficients
Anova(regression_bill_depth_3,type = 3) |> 
  tidy() |> slice(1:4) |> 
  mutate(p.value=formatP(p.value,ndigits=3, mark=TRUE))

# Post-hoc pairwise comparisons for species and sex
posthoc_species_2 <- glht(regression_bill_depth_2, linfct = mcp(species = "Tukey"))
posthoc_sex_2 <- glht(regression_bill_depth_2, linfct = mcp(sex = "Tukey"))


# Summary of results
summary(posthoc_species_2)
summary(posthoc_sex_2)
