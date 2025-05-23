pacman::p_load(conflicted,
               tidyverse,
               wrappedtools, 
               #palmerpenguins,
               ggfortify, GGally, nnet,
               caret,randomForest,kernlab,naivebayes,
               mlbench)
conflict_scout()
#conflict_prefer(name = "filter", winner = "dplyr")
conflicts_prefer(
  dplyr::filter(),
  ggplot2::alpha,
  dplyr::combine,
  dplyr::slice)

rawdata <- readRDS("../R_files/data/cervical_unstared.RDS")

# Preprocess the data
rawdata <- rawdata[, c("Tissuetype", names(rawdata)[4:ncol(rawdata)])]

print(head(rawdata))

#apply preprocessing (output will be a preProcess object with transformed data)
preProc <- caret::preProcess(rawdata, method = c("nzv", "YeoJohnson", "range"))

# Extract the preprocessed data from the preProcess object (otherwise slice doesn't work)
rawdata <- predict(preProc, rawdata)

# Set seed for reproducibility
set.seed(2001)

# Split the data into training and test sets
trainindex <- createDataPartition(
  y = rawdata$Tissuetype,  # Use the processed data for partitioning
  times = 1,
  p = 0.75,
  list = FALSE
)

# Train and test data
traindata <- rawdata[trainindex, ]
testdata <- rawdata[-trainindex, ]

# Train control parameters
ctrl <- trainControl(method = "repeatedcv", 
                     number = 5,  # 5-fold cross-validation
                     repeats = 25)  # 25 repeats



# define method-specific options ####
tune <- expand.grid(k=seq(1,9,2))
# tune a model ######
knnfit <- train(form = Tissuetype~.,
                data = traindata,
                method='knn',
                metric='Accuracy',
                trControl = ctrl,
                tuneGrid=tune)

tune <- expand.grid(
  nrounds = c(5, 10, 50),   # Number of boosting rounds (trees)
  max_depth = seq(5, 15, 5),    # Maximum depth of trees
  eta = 1,                      # Learning rate
  gamma = c(0.01, 0.001),       # Regularization parameter
  colsample_bytree = 1,         # Fraction of columns to sample for each tree
  min_child_weight = 1,         # Minimum sum of instance weight for child
  subsample = 1) # Fraction of data to sample for each tree

xgbfit <- train(form = Tissuetype~., # Formula: predicting Tissuetype
                data = rawdata,
                 preProcess = c(#"nzv",
                #                "YeoJohnson",
                                "corr"),
                #                "range"),
                method='xgbTree',
                objective = "binary:logistic", # Binary logistic regression (2-class classification)
                metric='Accuracy',
                trControl = ctrl,
                tuneGrid=tune,
                verbosity = 0)

ldafit <- train(form = Tissuetype~.,
                data = rawdata,
                 preProcess = c(
                #                "YeoJohnson",
                                "corr"),
                #                "range"),
                method='lda',
                metric='Accuracy',
                trControl = ctrl)

# define method-specific options ####
tune <- expand.grid(mtry=seq(2,3,1))
rffit <- train(form = Tissuetype~.,
               data = rawdata,
               #preProcess = c('center','scale'),
               method='rf',
               metric='Accuracy',
               trControl = ctrl,
               tuneGrid=tune)

svmfit <- train(form = Tissuetype~.,
                data = rawdata,
                #preProcess = c('center','scale'),
                method='svmLinear',
                metric='Accuracy',
                trControl = ctrl)

bayesfit <- train(form = Tissuetype~.,
                  data = rawdata,
                  #preProcess = c('center','scale'),
                  method='naive_bayes',
                  metric='Accuracy',
                  trControl = ctrl)

nnfit <- train(form=Tissuetype~.,
               data=rawdata, 
               method="nnet",
               #preProcess=c("center","scale"),
               metric="Accuracy",
               trControl=ctrl)

#combine models for comparison-----


resamps <- resamples(list(knn = knnfit, 
                          lda = ldafit,
                          rf=rffit,
                          xgb=xgbfit,
                          svm=svmfit,
                          bayes=bayesfit,
                          nn=nnfit))
summary(resamps)
resamps$values |> head()
resamps$values |> 
  pivot_longer(contains('~'),
               names_to = c('Model','Measure'),
               names_sep = '~') |> 
  ggplot(aes(Model,value))+
  geom_boxplot(outlier.alpha = .6)+
  facet_wrap(facets = vars(Measure),scales = 'free')
diffs <- diff(resamps)
summary(diffs)


#access individual models------
xgbfit
xgbfit[["bestTune"]]
knnfit[["bestTune"]]

# use xgbfit to predict test data ####
xgbpred <- predict(xgbfit, newdata = testdata)
confusionMatrix(xgbpred, testdata$Tissuetype)
# Importance
importance_xgb <- 
  xgboost::xgb.importance(model=xgbfit[["finalModel"]]) |> 
  arrange(Gain) |> 
  mutate(Feature=fct_inorder(Feature))
importance_xgb
importance_xgb |> 
  ggplot(aes(Feature,Gain))+
  geom_col(aes(fill=Gain))+
  coord_flip()+
  guides(fill="none")

# define method-specific options ####
tune <- expand.grid(mtry=seq(2,3,1))
rffit <- train(form = Tissuetype~.,
               data = traindata,
               preProcess = c("nzv",
                              "YeoJohnson",
                              "corr",
                              "range"),
               method='rf',
               metric='Accuracy',
               trControl = ctrl,
               tuneGrid=tune)

svmfit <- train(form = Tissuetype~.,
                data = traindata,
                preProcess = c("nzv",
                               "YeoJohnson",
                               "corr",
                               "range"),
                method='svmLinear',
                metric='Accuracy',
                trControl = ctrl)

bayesfit <- train(form = Tissuetype~.,
                  data = traindata,
                  preProcess = c("nzv",
                                 "YeoJohnson",
                                 "corr",
                                 "range"),
                  method='naive_bayes',
                  metric='Accuracy',
                  trControl = ctrl)

nnfit <- train(form=Tissuetype~.,
               data=traindata, 
               method="nnet",
               preProcess=c("nzv",
                            "YeoJohnson",
                            "corr",
                            "range"),
               metric="Accuracy",
               trControl=ctrl)

