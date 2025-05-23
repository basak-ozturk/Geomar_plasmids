pacman::p_load(conflicted,tidyverse,wrappedtools,
               flextable, here, palmerpenguins, ggbeeswarm, patchwork)
set_flextable_defaults(font.size = 9, 
                       padding.bottom = 1, 
                       padding.top = 3,
                       padding.left = 3,
                       padding.right = 4
)
data(penguins)
penguins <- na.omit(penguins)

#variable prep-------------

#count the number of species in the dataset

species_counts <- table(penguins$species)

species_counts

#make subtable with only Adelie penguins

adelie_penguins <- penguins[penguins$species == "Adelie", ]

adelie_penguins


numeric_cols <- sapply(penguins, is.numeric)  # vector of numeric columns
numeric_col_names <- names(penguins)[numeric_cols]  # names of numeric columns

# empty list to store plots
plot_list <- list()

# loop over each numeric column and create density plots -------------
for (col in numeric_col_names) {
  p <- ggplot(adelie_penguins, aes_string(x = col, fill = "sex")) +  # Use aes_string for variable column names
    geom_density(alpha = 0.4) +
    labs(title = paste("Density Plot of", col), x = col, y = "Density") +
    theme_minimal()
  
  plot_list[[col]] <- p  # Store the plot in the list
}

# combine all plots using patchwork
combined_plot <- plot_list[[1]]
for (i in 2:length(plot_list)) {
  combined_plot <- combined_plot + plot_list[[i]] + plot_layout(ncol = 2)
}
combined_plot


# stat summary for Gaussian body_mass_g and flipper_length_m -------------

ggplot(adelie_penguins,aes(x=sex,y=body_mass_g))+
  geom_beeswarm(size=3)+
  stat_summary(color="red",size=1.2,alpha=.7, 
               fun.data="mean_se",fun.args=list(mult=2))+
  ylab("size (mean \u00B1 2*SEM)")

ggplot(penguins,aes(x = flipper_length_mm,fill=sex))+
  geom_density(alpha=.4)

ggplot(adelie_penguins,aes(x=sex,y=flipper_length_mm))+
  geom_beeswarm(size=3)+
  stat_summary(color="blue",size=1.2,alpha=.7, 
               fun.data="mean_se",fun.args=list(mult=2))+
  ylab("size (mean \u00B1 2*SEM)")

adelie_penguins |> 
  group_by(sex) |> 
  summarize(MeanSE=meanse(body_mass_g),
            SD=sd(body_mass_g))

var_test_body_mass <- var.test(body_mass_g ~ sex, data = adelie_penguins)

var_test_body_mass

var_test_flipper_length <- var.test(flipper_length_mm ~ sex, data = adelie_penguins)

var_test_flipper_length

tOut_mass<-t.test(adelie_penguins$body_mass_g~adelie_penguins$sex)
tOut_mass$p.value

adelie_penguins |> 
  group_by(sex) |> 
  summarize(MeanSE=meanse(flipper_length_mm),
            SD=sd(flipper_length_mm))


tOut_flipper<-t.test(adelie_penguins$flipper_length_mm~adelie_penguins$sex)
tOut_flipper$p.value

#Wilcoxon test for ordinal measurements bill_length_mm and bill_depth_mm -------------

wilcox.test(bill_length_mm~sex,exact=F,correct=F,
            data=adelie_penguins)
wilcox.test(bill_depth_mm~sex,exact=F,correct=F,
            data=adelie_penguins)

ggplot(adelie_penguins,aes(sex,bill_length_mm))+
  geom_boxplot()+
  geom_beeswarm(alpha=.5)

ggplot(adelie_penguins,aes(sex,bill_depth_mm))+
  geom_boxplot()+
  geom_beeswarm(alpha=.5)

