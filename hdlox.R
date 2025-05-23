# Load necessary libraries
pacman::p_load(conflicted, 
               tidyverse, 
               car,
               readxl,
               ggplot2,
               flextable,
               corrplot,
               broom,
               rpart,
               rpart.plot,
               gt,
               beeswarm
               )

file_path <- "data/DatensatzHDLox.xlsx"
hd_data <- read_excel(file_path)

print(colnames(hd_data))

colnames(hd_data) <- colnames(hd_data) |>
  str_replace_all("\\(.*?\\)", "") |> # Remove anything in parentheses
  str_replace_all("[^[:alnum:]\\s]", "") |> # Remove non-alphanumeric characters
  str_trim() |> # Trim whitespace
  str_replace_all("\\s+", "_") |> # Replace spaces with underscores
  tolower() # Convert to lowercase

print(colnames(hd_data))

colnames(hd_data)[duplicated(colnames(hd_data))]

colnames(hd_data)[duplicated(colnames(hd_data))] <- "hbalc_duplicate"


# Remove empty columns


hd_filtered <- hd_data |>
  select(where(~ !all(is.na(.))))

hd_filtered <- hd_filtered[!is.na(hd_filtered$hdlox), ]

hd_filtered <- hd_filtered |>
  mutate(glukose = ifelse(glukose == 500, 50, glukose))

hd_filtered <- hd_filtered |>
  mutate(größe = ifelse(größe == 176.00, 1.76, größe))

# # Density plot for HDLox by Health Status
# ggplot(hd_filtered, aes(x = hdlox_no_unit, fill = factor(healthy_subjects))) +
#   geom_density(alpha = 0.5) +
#   labs(
#     title = "Density Plot of HDLox by Health Status",
#     x = "HDLox (no unit)",
#     y = "Density",
#     fill = "Health Status"
#   ) +
#   scale_fill_manual(
#     values = c("0" = "lightblue", "1" = "lightcoral"),
#     labels = c("0" = "Not Healthy", "1" = "Healthy")  
#   ) +
#   theme_minimal()
# 
# # Density plot for ldl by health status
# ggplot(hd_filtered, aes(x = ldl_in_mgdl_7, fill = factor(healthy_subjects))) +
#   geom_density(alpha = 0.5) +
#   labs(
#     title = "Density Plot of LDL by Health Status",
#     x = "LDL in mg/dl (Umrechnungsfaktor 38,66)",
#     y = "Density",
#     fill = "Health Status"
#   ) +
#   scale_fill_manual(
#     values = c("0" = "lightblue", "1" = "lightcoral"),
#     labels = c("0" = "Not Healthy", "1" = "Healthy")  
#   ) +
#   theme_minimal()
# 
# # Density plot for BMI by Health Status
# 
# ggplot(hd_filtered, aes(x = bmi, fill = factor(healthy_subjects))) +
#   geom_density(alpha = 0.5) +
#   labs(
#     title = "Density Plot of BMI by Health Status",
#     x = "BMI",
#     y = "Density",
#     fill = "Health Status"
#   ) +
#   scale_fill_manual(
#     values = c("0" = "lightblue", "1" = "lightcoral"),
#     labels = c("0" = "Not Healthy", "1" = "Healthy")  
#   ) +
#   theme_minimal()
# 
# # Density plot for Triglycerides by Health Status
# ggplot(hd_filtered, aes(x = triglyceride_in_mmoll8, fill = factor(healthy_subjects))) +
#   geom_density(alpha = 0.5) +
#   labs(
#     title = "Density Plot of Triglycerides by Health Status",
#     x = "Triglyceride in mmol",
#     y = "Density",
#     fill = "Health Status"
#   ) +
#   scale_fill_manual(
#     values = c("0" = "lightblue", "1" = "lightcoral"),
#     labels = c("0" = "Not Healthy", "1" = "Healthy")  
#   ) +
#   theme_minimal()
# 
# # Density plot for Albumin by Health Status
# ggplot(hd_filtered, aes(x = albumin_im_spoturin_mggkrea, fill = factor(healthy_subjects))) +
#   geom_density(alpha = 0.5) +
#   labs(
#     title = "Density Plot of Albumin by Health Status",
#     x = "Albumin im Spoturin. mg/gKrea",
#     y = "Density",
#     fill = "Health Status"
#   ) +
#   scale_fill_manual(
#     values = c("0" = "lightblue", "1" = "lightcoral"),
#     labels = c("0" = "Not Healthy", "1" = "Healthy")  
#   ) +
#   theme_minimal()
# 
# # Density plot for egfr by Health Status
# ggplot(hd_filtered, aes(x = egfr_mlmin173m2, fill = factor(healthy_subjects))) +
#   geom_density(alpha = 0.5) +
#   labs(
#     title = "Density Plot of EGFR by Health Status",
#     x = "eGFR, ml/min/1,73m2",
#     y = "Density",
#     fill = "Health Status"
#   ) +
#   scale_fill_manual(
#     values = c("0" = "lightblue", "1" = "lightcoral"),
#     labels = c("0" = "Not Healthy", "1" = "Healthy")  
#   ) +
#   theme_minimal()

# # Density plot for hdl by Health Status
# ggplot(hd_filtered, aes(x = hdl_in_mmoll, fill = factor(healthy_subjects))) +
#   geom_density(alpha = 0.5) +
#   labs(
#     title = "Density Plot of HDL by Health Status",
#     x = "HDL",
#     y = "Density",
#     fill = "Health Status"
#   ) +
#   scale_fill_manual(
#     values = c("0" = "lightblue", "1" = "lightcoral"),
#     labels = c("0" = "Not Healthy", "1" = "Healthy")
#   ) +
#   theme_minimal()

# Define the columns to plot and test for automatization. This will be used for density plots and normality tests. Reshape the data to long format. We will need this for plotting.

# Define the columns to be tested
cols <- c("hdlox_no_unit", "ldl_in_mgdl_7", "bmi", 
          "triglyceride_in_mmoll8", "albumin_im_spoturin_mggkrea", 
          "egfr_mlmin173m2", "hscrp", "pla2_47", "glukose", "größe", "age", "cholesterin_in_mmoll9", "hdl_in_mmoll")

hd_filtered[cols] <- lapply(hd_filtered[cols], function(x) as.numeric(as.character(x)))


# Loop through the columns to perform the Shapiro-Wilk test
shapiro_results <- data.frame(Column = character(), P_Value = numeric())

for (col in cols) { # col is a function name, better to change. Don't put spaces between function names and round bracket
  
  
  test_result <- shapiro.test(hd_filtered[[col]])
  

  shapiro_results <- rbind(shapiro_results, data.frame(Column = col, P_Value = test_result$p.value))
}

print(shapiro_results)

#Facet plot (add box plots)

hd_long <- hd_filtered |>
  select(healthy_subjects, all_of(cols)) |>
  pivot_longer(cols = all_of(cols), names_to = "Variable", values_to = "Value")


head(hd_long)

# Density plot for variables by Health Status
ggplot(hd_long, aes(x = Value, fill = factor(healthy_subjects))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Density Plots of  Variables by Health Status",
    x = "Value",
    y = "Density",
    fill = "Health Status"
  ) +
  scale_fill_manual(
    values = c("0" = "lightblue", "1" = "lightcoral"),
    labels = c("0" = "Not Healthy", "1" = "Healthy")
  ) +
  facet_wrap(~ Variable, scales = "free", ncol = 3) +
  theme_minimal()



ggplot(hd_long, aes(x = factor(healthy_subjects), y = Value, fill = factor(healthy_subjects))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Suppress outliers to avoid overlap with jitter
  geom_jitter(aes(color = factor(healthy_subjects)), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              size = 1, alpha = 0.3) +
  labs(
    title = "Box Plots of Variables with Beeswarm by Health Status",
    x = "Health Status",
    y = "Value",
    fill = "Health Status",
    color = "Health Status"
  ) +
  scale_fill_manual(
    values = c("0" = "lightblue", "1" = "lightcoral"),
    labels = c("0" = "Not Healthy", "1" = "Healthy")
  ) +
  scale_color_manual(
    values = c("0" = "blue", "1" = "red"),
    labels = c("0" = "Not Healthy", "1" = "Healthy")
  ) +
  facet_wrap(~ Variable, scales = "free", ncol = 3) +
  theme_minimal()


# None of the columns are normally distributed. Log transform the relevant columns, adding a small constant to avoid issues with zeros. Write the log-transformed values out into new columns

# hd_filtered$log_hdlox_no_unit <- log(hd_filtered$hdlox_no_unit + 0.0001)
# hd_filtered$log_ldl_in_mgdl_7 <- log(hd_filtered$ldl_in_mgdl_7 + 0.0001)
# hd_filtered$log_bmi <- log(hd_filtered$bmi + 1)
# hd_filtered$log_triglyceride_in_mmoll8 <- log(hd_filtered$triglyceride_in_mmoll8 + 0.0001)
# hd_filtered$log_albumin_im_spoturin_mggkrea <- log(hd_filtered$albumin_im_spoturin_mggkrea + 0.0001)
# hd_filtered$log_egfr_mlmin173m2 <- log(hd_filtered$egfr_mlmin173m2 + 0.0001)
# hd_filtered$log_hscrp <- log(hd_filtered$hscrp + 0.0001)
# hd_filtered$log_pla2_47 <- log(hd_filtered$pla2_47 + 0.0001)
# hd_filtered$log_glukose <- log(hd_filtered$glukose + 0.0001)
# 
# List of log-transformed columns

log_cols <- c("hdlox_no_unit", "ldl_in_mgdl_7", "bmi", 
              "triglyceride_in_mmoll8", "albumin_im_spoturin_mggkrea", 
              "egfr_mlmin173m2", "hscrp", "pla2_47")

# Add log-transformed columns
hd_filtered <- hd_filtered |>
  mutate(across(all_of(log_cols), ~ log(. + 0.0001), .names = "log_{.col}"))
log_shapiro_results <- data.frame(Column = character(), P_Value = numeric())

# Loop through the columns and run the Shapiro-Wilk test
for (col in log_columns) {
  
  log_test_result <- shapiro.test(hd_filtered[[col]])
  

  log_shapiro_results <- rbind(log_shapiro_results, data.frame(Column = col, P_Value = log_test_result$p.value))
}

print(log_shapiro_results)

#Facet plot of log-transformed data

hd_long_log <- hd_filtered |>
  select(healthy_subjects, all_of(log_columns)) |>
  pivot_longer(cols = all_of(log_columns), names_to = "Variable", values_to = "Value")


head(hd_long_log)

# Density plot for log-transformed variables by Health Status
ggplot(hd_long_log, aes(x = Value, fill = factor(healthy_subjects))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Density Plots of Log-Transformed Variables by Health Status",
    x = "Value",
    y = "Density",
    fill = "Health Status"
  ) +
  scale_fill_manual(
    values = c("0" = "lightblue", "1" = "lightcoral"),
    labels = c("0" = "Not Healthy", "1" = "Healthy")
  ) +
  facet_wrap(~ Variable, scales = "free", ncol = 3) +
  theme_minimal()


# HDlox outcome group measurements for different diseases-----------

# Density plot for HDLox by sex
ggplot(hd_filtered, aes(x = hdlox_no_unit, fill = factor(gender))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Density Plot of HDLox by Sex",
    x = "HDLox (no unit)",
    y = "Density",
    fill = "Sex"
  ) +
  scale_fill_manual(
    values = c("1" = "lightblue", "2" = "lightcoral"),
    labels = c("1" = "Male", "2" = "Female")
  ) +
  theme_minimal()


# Density plot for HDLox by pavkcavk_aktuell
ggplot(hd_filtered, aes(x = hdlox_no_unit, fill = factor(pavkcavk_aktuell))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Density Plot of HDLox by pAVK/cAVK",
    x = "HDLox (no unit)",
    y = "Density",
    fill = "pAVK/cAVK"
  ) +
  scale_fill_manual(
    values = c("0" = "lightblue", "1" = "lightcoral"),
    labels = c("0" = "No", "1" = "Yes")
  ) +
  theme_minimal()

# Density plot for HDLox by koro_im_letzten_jahr
ggplot(hd_filtered, aes(x = hdlox_no_unit, fill = factor(koro_im_letzten_jahr))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Density Plot of HDLox by koro_im_letzten_jahr",
    x = "HDLox (no unit)",
    y = "Density",
    fill = "koro_im_letzten_jahr"
  ) +
  scale_fill_manual(
    values = c("0" = "lightblue", "1" = "lightcoral"),
    labels = c("0" = "No", "1" = "Yes")
  ) +
  theme_minimal()

# Density plot for HDLox by herzinfarkt
ggplot(hd_filtered, aes(x = hdlox_no_unit, fill = factor(herzinfarktsteminstemi))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Density Plot of HDLox by Heart Attack",
    x = "HDLox (no unit)",
    y = "Density",
    fill = "Heart attack"
  ) +
  scale_fill_manual(
    values = c("0" = "lightblue", "1" = "lightcoral"),
    labels = c("0" = "No", "1" = "Yes")
  ) +
  theme_minimal()

# Density plot for HDLox by koronarsklerose_aktuell
ggplot(hd_filtered, aes(x = hdlox_no_unit, fill = factor(koronarsklerose_aktuell))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Density Plot of HDLox by Koronarsklerose",
    x = "HDLox (no unit)",
    y = "Density",
    fill = "Koronarsklerose"
  ) +
  scale_fill_manual(
    values = c("0" = "lightblue", "1" = "lightcoral"),
    labels = c("0" = "No", "1" = "Yes")
  ) +
  theme_minimal()

# Density plot for HDLox by herzinsuffizienz
ggplot(hd_filtered, aes(x = hdlox_no_unit, fill = factor(herzinsuffizienz))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Density Plot of HDLox by Herzinsuffizienz",
    x = "HDLox (no unit)",
    y = "Density",
    fill = "herzinsuffizienz"
  ) +
  scale_fill_manual(
    values = c("0" = "lightblue", "1" = "lightcoral"),
    labels = c("0" = "No", "1" = "Yes")
  ) +
  theme_minimal()

# Plot of HDLox by age and health status
ggplot(hd_filtered, aes(x = age, y = hdlox, color = factor(healthy_subjects))) +
  geom_point(alpha = 0.7) +
  labs(
    title = "Scatter Plot of Age vs. HDLox by Health Status",
    x = "Age",
    y = "HDLox (no unit)",
    color = "Health Status"
  ) +
  scale_color_manual(
    values = c("0" = "lightblue", "1" = "lightcoral"),
    labels = c("0" = "Not Healthy", "1" = "Healthy")
  ) +
  theme_minimal()

# Density plot for HDLox by peripheralcarotid_artery_disease
ggplot(hd_filtered, aes(x = hdlox_no_unit, fill = factor(peripheralcarotid_artery_disease))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Density Plot of HDLox by peripheral carotid artery disease",
    x = "HDLox (no unit)",
    y = "Density",
    fill = "peripheral carotid artery disease"
  ) +
  scale_fill_manual(
    values = c("0" = "lightblue", "1" = "lightcoral"),
    labels = c("0" = "No", "1" = "Yes")
  ) +
  theme_minimal()



# Statistical tests on peripheral carotid artery disease as outcome -------------

hd_filtered$peripheralcarotid_artery_disease <- as.factor(hd_filtered$peripheralcarotid_artery_disease) #we need factors for Levene's test

# Check variance equality. Take the log transformed HDLox values as they're more normal
leveneTest(log_hdlox_no_unit ~ peripheralcarotid_artery_disease, data = hd_filtered)

# Perform t-test with equal variences set to true as per the outcome of Levene test
t_test_result<- t.test(log_hdlox_no_unit ~ peripheralcarotid_artery_disease, data = hd_filtered, var.equal = TRUE)


print(t_test_result)

# Extract results from t-test in order to put into a table
t_stat <- t_test_result$statistic
df <- t_test_result$parameter
p_value <- t_test_result$p.value
conf_interval_log <- t_test_result$conf.int
means_log <- t_test_result$estimate

# Extrapolated results to the original scale (since I used log-transformed results for the t-test)
means_original <- exp(means_log)
conf_interval_original <- exp(conf_interval_log)

# Create a data frame to store the results
results_df <- data.frame(
  Metric = c("t-statistic", "Degrees of Freedom", "p-value", 
             "Log Mean (Group 0)", "Log Mean (Group 1)",
             "Original Mean (Group 0)", "Original Mean (Group 1)",
             "Log CI Lower", "Log CI Upper",
             "Original CI Lower", "Original CI Upper"),
  Value = c(
    t_stat, df, p_value,
    means_log[1], means_log[2],
    means_original[1], means_original[2],
    conf_interval_log[1], conf_interval_log[2],
    conf_interval_original[1], conf_interval_original[2]
  )
)

#results_table <- flextable(results_df)

results_table <- flextable(results_df) |>
  set_header_labels(Metric = "Metric", Value = "Value") |>
  autofit()


results_table

# Relation between HDLox and other measurements--------------

#I added some more measurements to those I explored before

gaussiancols <- c("log_hdlox_no_unit", "größe", "glukose", "bmi", "pla2_47", "egfr_mlmin173m2" )

ordinalcols <- c("hdlox_no_unit", "ldl_in_mgdl_7",  
                 "triglyceride_in_mmoll8", "albumin_im_spoturin_mggkrea", 
                  "hscrp", "cholesterin_in_mmoll9", "hdl_in_mmoll")

# Compute correlations for Gaussian variables. use = "complete.obs" tells R to only use rows that do not contain any missing values for the specified columns. Pairwise complete would be another one (maybe better one, use="pair"). Check cor.test to get P values for the pairwise comparisons of HDLox and other measurements. With complete.obs, all measurements have exactly the same sample size.

gaussian_correlations <- cor(hd_filtered[, gaussiancols], use = "complete.obs")
gaussian_correlations


corrplot(gaussian_correlations, method = "circle", type = "upper", 
         tl.col = "black", tl.srt = 45, # Adjust the text angle for easy reading
         title = "Correlation Matrix of Gaussian Measurements")


# Compute Spearman correlations for ordinal variables
ordinal_correlations <- cor(hd_filtered[, ordinalcols], method = "spearman", use = "complete.obs")


png("data/spearmann_corrplot.png", width = 800, height = 600)

corrplot(ordinal_correlations, 
         method = "number", 
         main = "Spearman Correlation Plot of Ordinal Measurements",
         tl.col = "black",
         tl.cex = 0.8, 
         tl.srt = 45,
         cex.main = 1,
         mar=c(0, 0, 1, 1),
         )

dev.off()

# Logistic regression model with log_HDlox--------------
glmmodel1 <- glm(peripheralcarotid_artery_disease ~ log_hdlox_no_unit, 
              data = hd_filtered, family = binomial)
summary(glmmodel1)

glmmodel_with_age <- glm(peripheralcarotid_artery_disease ~ log_hdlox_no_unit + age, 
                 data = hd_filtered, family = binomial)
summary(glmmodel_with_age)


# Summarize models
model1_summary <- tidy(glmmodel1)
model2_summary <- tidy(glmmodel_with_age)

results_table2 <- rbind(
  cbind(model = "Model 1 (HDlox only)", model1_summary),
  cbind(model = "Model 2 (HDlox + Age)", model2_summary)
)

#Trying gt for making tables instead of flextable

results_table2 |>
  gt() |>
  tab_header(
    title = "Logistic Regression Results",
    subtitle = "Effect of HDlox on Peripheral Carotid Artery Disease"
  )

results_table2

#Anova Chisq test for goodness of fit
anova_hdlox <- anova(glmmodel1, glmmodel_with_age, test = "Chisq")

anova_df <- as.data.frame(anova_hdlox)

anova_df <- cbind(
  Model = c("peripheral carotid_artery disease ~ log HDLox", "peripheral carotid_artery disease ~ log HDLox+age"),
  anova_df
)

anova_table <- anova_df |>
  gt() |>
  tab_header(
    title = "Chi-Squared Test for Goodness of Fit",
    
  ) |>
  fmt_number(
    columns = c("Resid. Df", "Resid. Dev", "Df", "Deviance", "Pr(>Chi)"),
    decimals = 4
  ) |>
  cols_label(
    Model = "Model",
    `Resid. Df` = "Residual DF",
    `Resid. Dev` = "Residual Deviance",
    Df = "Degrees of Freedom",
    Deviance = "Deviance",
    `Pr(>Chi)` = "p-value"
  ) |>
  cols_align(align = "center", columns = everything())


anova_table

#regression tree for machine learning-----------

# Merge columns to be tested
selected_columns <- c(gaussiancols, ordinalcols) #need to check again for double HDLox

# Fit the model using only the selected columns
formula <- as.formula(paste("peripheralcarotid_artery_disease ~", paste(selected_columns, collapse = " + ")))

tree_model_pcad <- rpart(formula, 
                    data = hd_filtered, 
                    method = "class", 
                    control = rpart.control(minsplit = 20, cp = 0.01)) #(minsplit = 20, cp = 0.01 seem like default values


summary(tree_model_pcad)

png("data/decision_tree_pcad.png", width = 800, height = 600)

rpart.plot(tree_model_pcad, extra = 101, cex = 0.8, mar = c(1, 1, 1, 1), main="Regression Tree")

dev.off()

