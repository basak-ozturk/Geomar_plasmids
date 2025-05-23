pacman::p_load(car, wrappedtools, tidyverse, ggbeeswarm, 
                pROC, palmerpenguins, broom, rpart, rpart.plot)

data(penguins)
penguins <- na.omit(penguins) # Remove NA values

str(penguins)

# Fit the model-----------
logreg_out <- glm(sex ~ bill_depth_mm + bill_length_mm + flipper_length_mm + body_mass_g + species, 
                  data = penguins, 
                  family = binomial())

# extract and transform model parameters to get odds ratios (ORs)
ORs <- exp(logreg_out$coefficients)
CIs <- exp(confint(logreg_out)) #transform body mass to kg for better results

CIs
# test model with ANOVA
Anova_out <- Anova(logreg_out, type = 2) |> 
  broom::tidy() |> 
  mutate(p.value = formatP(p.value, ndigits = 5))
Anova_out
#test each OR
sum_out <- summary(logreg_out)
sum_out
# data for forest plot---------------

# tidy tibble of ORs and CIs with labels
forest_data <- tibble(
  term = names(ORs)[-1],  # exclude the intercept term
  OR = ORs[-1],
  CI_lower = CIs[-1,1],
  CI_upper = CIs[-1,2],
  p = summary(logreg_out)$coefficients[-1, 4]  # extract p-values excluding the intercept
)

forest_data <- forest_data |>
  mutate(

    # format species labels 
    term = str_replace_all(term, "speciesGentoo", "Species (Gentoo)"),
    term = str_replace_all(term, "speciesChinstrap", "Species (Chinstrap)"),
    term = str_replace_all(term, "speciesAdelie", "Species (Adelie)"),
    term = str_replace_all(term, "flipper_length_mm", "Flipper Length (mm)"),
    term = str_replace_all(term, "bill_depth_mm", "Bill Depth (mm)"),
    term = str_replace_all(term, "bill_length_mm", "Bill Length (mm)"),
    term = str_replace_all(term, "body_mass_g", "Body Mass (g)"),
  )

# add significance labels based on p-value threshold ( p < 0.05)
forest_data <- forest_data |>
  mutate(
    Significance = if_else(p < 0.05, "*", "ns"),  # "*" for significant, "ns" for non-significant
    Label = paste(term, "(", Significance, ")")  
  )

forest_data

# plot version 1 with log scale-----------
baseplot <- ggplot(forest_data, aes(x = Label, y = OR)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper)) + 
  geom_hline(yintercept = 1, linewidth = 0.2, linetype = 2) + 
  coord_flip()

# log scale
baseplot + 
  scale_y_log10(breaks=logrange_1, 
                minor_breaks=logrange_123456789 )+ 
  geom_text(aes(label=Significance), vjust=1.5,color='red')+ 
  ggtitle('OddsRatios shown on log-scale')+ 
  xlab(NULL)

#plot version 2 with linear scale

baseplot2 <- ggplot(forest_data, aes(x = Label, y = OR)) +
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper)) + 
  geom_hline(yintercept = 1, linewidth = 0.2, linetype = 2) + 
  coord_flip()

# Remove scale_y_log10 and just use the default linear scale
baseplot2 + 
  geom_text(aes(label = Significance), vjust = 1.5, color = 'red') + 
  ggtitle('Odds Ratios with 95% Confidence Intervals') + 
  xlab(NULL)

penguins$pGLM <- predict(logreg_out, type = 'response')  # Predict probability 0-1

#ROC-----------
roc_out <- roc(response = penguins$sex, predictor = penguins$pGLM)

# determine the optimal cutoff using Youden's index
youden <- coords(roc_out, x = 'best', best.method = 'youden')

print(youden)

# plot predictions 
penguins |> 
  mutate(`prediction quality`=
           sex_when(sex=="female" &
                       pGLM<youden$threshold ~ 
                       "false negative",
                     sex=="male" & 
                       pGLM>=youden$threshold 
                     ~ "false positive", 
                     .default = 'correct' )) |>
  ggplot(aes(sex,pGLM))+ 
  geom_boxplot(outlier.alpha = 0)+ 
  scale_y_continuous(breaks=seq(0,1,.1))+ 
  geom_beeswarm(alpha=.75, 
                aes(color=`prediction quality`))+ 
  scale_color_manual(values=c("seagreen","firebrick","magenta"))+ 
  geom_hline(yintercept = c(.5, youden$threshold,.8),
             color='red',
             linetype=2:4)+
  annotate(geom = "label",
           x = 1,y=youden$threshold, 
           label=paste("Youden-cutoff:",
                       roundR(youden$threshold)),
           hjust=1.2,vjust=0.25)+
  theme(legend.position="bottom")

#machine learning-------------


#check columns for predictors
predvars <- ColSeeker(namepattern = c("mm", "body", "species"), data = penguins)

#create formula for rpart
rtformula <- paste("sex ~", paste(predvars$bticked, collapse = "+"))

#run the regression tree model
regtree_out <- rpart(rtformula, minsplit = 10, cp = 0.001, data = penguins)

#plot the tree
rpart.plot(regtree_out, type = 2, tweak = 2.0, varlen = 4, faclen = 5, leaf.round = 0)

#calculate importance and plot
importance <- as_tibble(regtree_out$variable.importance, rownames = 'Predictor') |>
  dplyr::rename('Importance' = 2) |>
  mutate(Predictor = fct_reorder(.f = Predictor, .x = Importance, .fun = min)) |>
  arrange(desc(Importance))

#plot importance
importance |>
  ggplot(aes(Predictor, Importance)) + 
  geom_col() + 
  coord_flip()

penguins$pRT <- predict(regtree_out)[,2]

#pROC 
roc_out_rt <- roc(response=penguins$sex,
                  predictor=penguins$pRT ) 
youden <- pROC::coords(roc_out_rt,x='best',
                       best.method='youden') 
youden 
ggroc(roc_out_rt,legacy.axes = T)+ 
  geom_abline(slope = 1,intercept = 0)+ 
  geom_point(x=1-youden$specificity,
             y=youden$sensitivity, color='red', size=2 )

ggroc(list(RTbased=roc_out_rt,GLM_based=roc_out),legacy.axes = T)+ 
  geom_abline(slope = 1,intercept = 0)

ggplot(penguins,aes(x=sex,y=pRT))+ 
  geom_boxplot(coef=3)+ 
  scale_y_continuous(breaks = seq(from = 0,to = 1,by = .1))+ 
  geom_hline(yintercept = c(.5,youden$threshold), 
             color=c('red',"blue"), linetype=2)+ 
  ggbeeswarm::geom_beeswarm(alpha=0.25) 
ggplot(penguins,aes(pGLM,pRT, color=sex,shape=sex))+ 
  geom_point(size=2)+ 
  scale_color_manual(values = c('darkgreen','red'))+ 
  scale_shape_manual(values = c(0,6))+ 
  stat_summary(fun.data=mean_cl_boot) 
ggplot(penguins,aes(x=sex,y=pRT))+ 
  geom_violin()+ 
  scale_y_continuous(breaks = seq(from = 0,to = 1,by = .1))+ 
  geom_hline(yintercept = .5,color='red')

#jackknife---------------

# pGLM: The original predicted probability of sex using the logistic regression model (glm).
# pGLM_JK: The jackknife prediction probability of sex from the logistic regression model, excluding each observation one at a time.
# pRT_JK: The jackknife prediction probability of sex from the regression tree model (rpart), also excluding each observation one at a time.
# pRT: Predictions from the regression tree 

penguins$pRT_JK <- NA_real_
penguins$pGLM_JK <- NA_real_

for(pat_i in 1:nrow(penguins)) { 
  tempdata <- penguins[-pat_i,]
  regtree_out_tmp <- rpart(rtformula, minsplit = 5, cp = 0.001, data = tempdata)
  penguins$pRT_JK[pat_i] <- predict(regtree_out_tmp, newdata = penguins[pat_i,])[, 2]
  
  glm_out_tmp <- glm(rtformula, family = binomial(), data = tempdata)
  penguins$pGLM_JK[pat_i] <- predict(glm_out_tmp, newdata = penguins[pat_i,], type = "response")
}

#visualize regression tree predictions
ggplot(penguins, aes(sex, pRT_JK)) + 
  geom_boxplot()

#pivot and plot to compare all predictions
penguins |> 
  dplyr::select(sex, pGLM, pGLM_JK, pRT_JK, pRT) |>
  pivot_longer(cols = c(pGLM, pGLM_JK, pRT_JK, pRT), 
               names_to = 'Analysis', 
               values_to = 'pAffected') |>
  ggplot(aes(sex, pAffected, fill = Analysis)) + 
  geom_boxplot()
