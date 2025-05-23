pacman::p_load(conflicted,
               tidyverse,
               wrappedtools,  # just tools 
               #palmerpenguins, # data
               #ggforce, # for cluster plots, hulls, zoom etc
               ggbeeswarm,
               flextable,
               caret, # Classification and Regression Training
               preprocessCore,  # pre-processing functions
               gmodels, # tools for model fitting
               easystats, 
               yardstick,
               pheatmap
               )  

# conflict_scout()
conflict_prefer('slice','dplyr')
conflict_prefer("filter", "dplyr")
conflict_prefer("anti_join", "dplyr")
conflict_prefer("group_by", "dplyr")

rawdata <- readRDS("../R_files/data/cervical_unstared.RDS")

# Remove columns with only zeros--------------
rawdata <- rawdata[, colSums(rawdata != 0) > 0]

# Keep SampleID and Tissuetype, remove PatID
numeric_data <- rawdata[, c("SampleID", "Tissuetype", names(rawdata)[4:ncol(rawdata)])]

# Separate SampleID and Tissuetype
sampleID <- numeric_data$SampleID
tissuetype <- numeric_data$Tissuetype

# Create a dataset with only numeric columns for preprocessing
numeric_data <- numeric_data[, 4:ncol(numeric_data)]


# contol prints
print(dim(numeric_data))

print(head(numeric_data))

numeric_data <- na.omit(numeric_data)

# Shapiro-Wilk test to all columns and extract only p-values--------------
normality_results <- numeric_data |>
  summarize(across(everything(), 
                   ~shapiro.test(.)$p.value,
                   .names = "{.col}_shapiro_p_value")) |>
  pivot_longer(everything(), 
               names_to = "variable", 
               values_to = "shapiro_p_value") |>
  mutate(variable = sub("_shapiro_p_value$", "", variable))

ggplot(normality_results, aes(x = shapiro_p_value)) +
  geom_histogram(binwidth = 0.001, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Shapiro-Wilk p-values",
       x = "Shapiro-Wilk p-value",
       y = "Count") +
  scale_x_log10()  # Use log scale for x-axis

summary(normality_results$shapiro_p_value)

# preProcess object for scaling
preprocessed <- caret::preProcess(numeric_data, method = c("nzv",
                       "YeoJohnson",
                      "corr",
                       "range"))

#apply the scaling to the data
scaled <- predict(preprocessed, numeric_data)

target <- as.factor(rawdata$Tissuetype)

#heatmap to see range 

subset_data <- scaled[, 1:50]  # first 50 columns

pheatmap(subset_data,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Heatmap of Cervical Cancer Data")

# combine preprocessed features with target
preprocessed_data <- cbind(SampleID = sampleID, Tissuetype = tissuetype, scaled)

# training and testing sets--------------
set.seed(444)  # for reproducibility


# create training set
traindata <- preprocessed_data|> 
  select(Tissuetype, SampleID)|>
  dplyr::group_by(Tissuetype)|>
  slice_sample(prop = 2/3) |>
  ungroup()#|>
  #select(-SampleID)

# #method1 to create test data
# testdata <- rawdata|>
#   filter(!SampleID %in% traindata$SampleID) |>
#   select(SampleID, Tissuetype, everything()) |>
#   select(-SampleID)

# Create test set. The anti-join operation in R is used to identify observations that exist in the first dataset but not in the second dataset.
testdata <- preprocessed_data |>
  dplyr::anti_join(traindata, by = "SampleID")

# verify the distribution
cat("Training set distribution:\n", table(traindata$Tissuetype), "\n")
cat("Test set distribution:\n", table(testdata$Tissuetype), "\n")

# drop unnecessary stuff
train_features <- traindata |> select(-SampleID, -Tissuetype)
train_target <- traindata$Tissuetype
train_sampleID <- traindata$SampleID  #keep this separate just in case

test_features <- testdata |> select(-SampleID, -Tissuetype)
test_target <- testdata$Tissuetype
test_sampleID <- testdata$SampleID 

# Perform KNN-------------

train_out <- knn3Train(
  train = train_features,   #matrix (training)
  test = test_features,     #matrix (testing)
  cl = train_target,        #target labels for training set
  k = 3                     #number of neighbors
)
str(train_out) 
head(train_out) 
# train_res <- 
#   attr(x = train_out,which = 'prob') |> 
#   as_tibble() |> 
#   mutate(predicted=factor(train_out)) |> 
#   bind_cols(test_features)

train_res <- 
  attr(x = train_out, which = 'prob') |> 
  as_tibble() %>% 
  mutate(predicted = factor(train_out)) |>  
  bind_cols(test_features) |>
  mutate(actual = test_target)  #  mutate actual after binding test_features 

# pivot the probabilities into long format for ggplot2
train_res_long <- train_res |>
  pivot_longer(
    cols = c("Control", "Tumor"),  # columns with probabilities (control or tumor)
    names_to = "Class",           # class labels
    values_to = "Probability"     # probability values
  )

sensitivity <- cm$byClass['Sensitivity']
print(sensitivity)
# create a violin plot with beeswarm
ggplot(train_res_long, aes(x = Class, y = Probability)) +
  geom_violin(fill = "coral", alpha = 0.7) +  
  geom_beeswarm(size = 1.5, alpha = 0.5) +     
  theme_minimal() +
  labs(
    title = "Predicted Probabilities for Each Class",
    x = "Predicted Class",
    y = "Probability"
  ) +
  facet_wrap(~ predicted, ncol = 1, labeller = "label_both") 

CrossTable(train_res$predicted,
           train_res$actual, prop.chisq = FALSE, prop.t = FALSE, 
           format = 'SPSS')

