#Load required libraries
library(data.table)
library(xgboost)
library(caret)
library(dplyr)
library(ggplot2)

dir_out <- "Centralized_train"

# Combine training data from 11 dataset
train_data_dir <- "Training_set/"
data_files <- list.files(path = train_data_dir, pattern = "*.csv", full.names = TRUE)
full_train_data <- bind_rows(lapply(data_files, read_csv))
print(dim(full_train_data))  # Rows: 530 Columns: 54 


#########################################################
## Train model with centralized data to estimate hyper parameters
#########################################################
set.seed(135)

trainVar <- full_train_data[, -ncol(full_train_data)]
trainResponse <- full_train_data$response

# Tuning hyper-parameters
hyperparam_grid <- expand.grid(
  nrounds = seq(from = 10, to = 300, by = 50),  # Number of boosting iterations
  eta = seq(0.01, 0.3, length.out = 5),       # Learning rate (10 equally spaced values)
  max_depth = c(3, 6, 9),                      # Maximum tree depth
  gamma = c(0, 0.1, 0.5, 1),                # Minimum loss reduction for splits
  colsample_bytree = c(0.6, 0.8, 1),           # Fraction of columns sampled
  min_child_weight = 1,               # Minimum sum of instance weights
  subsample = c(0.6, 0.8, 1)                   # Subsampling ratio
)

# Define control function for train
control <- trainControl(method = "cv", 
                        number = 5, 
                        summaryFunction = twoClassSummary, 
                        classProbs = TRUE,
                        verboseIter = FALSE,
                        allowParallel = FALSE) 

# Train the model using caret's train function with xgboost method
bst_model <- train(x = trainVar, 
                   y = trainResponse, 
                   method = "xgbTree", 
                   trControl = control,
                   metric = "ROC",
                   tuneGrid = hyperparam_grid)

####################################################################
# Using Caret for the best hyper-parameters to get final model
####################################################################
group <-  ifelse(trainResponse == "R", 0, 1) 
xgb_model <- xgboost(data = as.matrix(trainVar),
                     label = group,
                     booster = "gbtree",
                     objective = "binary:logistic",
                     nrounds = bst_model$bestTune$nrounds,
                     max_depth= bst_model$bestTune$max_depth,
                     colsample_bytree = bst_model$bestTune$colsample_bytree,
                     min_child_weight = bst_model$bestTune$min_child_weight,
                     subsample = bst_model$bestTune$subsample,
                     eta = bst_model$bestTune$eta,
                     gamma = bst_model$bestTune$gamma,
                     verbose = 0)

final_model <- xgb_model
importance_matrix <- xgb.importance(colnames(trainVar), model = final_model)

xgb.plot.importance(importance_matrix)

mv.res <- list(final_model = final_model,
               importance_matrix = importance_matrix,
               bestTune = bst_model$bestTune)

saveRDS(final_model, file.path(dir_out, "Centralized_trainModel.rds"))


