# Script to train XGBoost using gene signature scores and clinical data for each study

library(reticulate)
library(xgboost)
library(dplyr)
library(sparklyr)
library(pROC)
library(caret)
library(pROC)
library(PRROC)

# Load the .Rda file
#data.path <- "/Users/nicole/PHD-Univerity of Toronto/PhD-BHK Work Folder/Insight/data" 

#data.name <- 'ICB_Gide__Melanoma__PD-(L)1.rda'

#study.name <- "ICB_Gide__Melanoma__PD-(L)1"


args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
study.name <- tools::file_path_sans_ext(basename(data_file))

# Load the data
load("Output/train_set_ICB_Gide__Melanoma__PD-(L)1.rda")
study.name <- "train_set_ICB_Gide__Melanoma__PD-(L)1"
if (!exists("train_dat")) {
  stop("Error: 'train_dat' object not found in the loaded file.")
}

model_out <- Sys.getenv("SPARK_LOCAL_DIRS", unset = "Model_2")


# conert 
#########################################################
## Train model to estimate hyper parameters
#########################################################
set.seed(135)

train_var <- train_dat[, -ncol(train_dat)]
train_response <- train_dat[, ncol(train_dat)]

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
bst_model <- train(x = train_var, 
                   y = train_response, 
                   method = "xgbTree", 
                   trControl = control,
                   metric = "ROC",
                   tuneGrid = hyperparam_grid)


# Find the optimal threshold for the best model
train_probs <- predict(bst_model, train_var, type = "prob")[, "R"] # Probability of response "R"
roc_curve <- roc(train_response, train_probs, levels = c("NR", "R"))

# threshold using Youden's J statistic
best_threshold <- mean(coords(roc_curve, "best", best.method = "youden")[, "threshold"])

# Evaluate the model performance at the optimal threshold
train_pred_labels <- ifelse(train_probs >= best_threshold, "R", "NR")
train_pred_labels <- factor(train_pred_labels, levels = c("NR", "R"))
conf_matrix <- confusionMatrix(train_pred_labels, train_response)

# AUROC
roc_obj <- roc(train_response, train_probs, levels = c("NR", "R"), direction = "<")
auroc <- as.numeric(auc(roc_obj))

# Compute AUPRC
pr_obj <- pr.curve(scores.class0 = train_probs[train_response == "R"], 
                   scores.class1 = train_probs[train_response == "NR"], 
                   curve = FALSE)
auprc <- pr_obj$auc.integral

accuracy <- conf_matrix$overall["Accuracy"] 
precision <- conf_matrix$byClass["Precision"] # specificity 
recall <- conf_matrix$byClass["Recall"] # sensitivity 
f1_score <- 2 * (precision * recall) / (precision + recall)

evaluation_metrics <- list(
  Accuracy = round(accuracy, 4),
  Precision = round(precision, 4),
  Recall = round(recall, 4),
  F1_score = round(f1_score, 4),
  AUROC = round(auroc, 4),
  AUPRC = round(auprc, 4)
)

# Using Caret for the best hyper-parameters to get final model
group <-  ifelse(train_response == "R", 1, 0) 
xgb_model <- xgboost(data = as.matrix(train_var),
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

local_model <- xgb_model
importance_matrix <- xgb.importance(colnames(train_var), model = local_model)

local_model <- list(local_model = local_model,
               importance_matrix = importance_matrix,
               bestTune = bst_model$bestTune,
               evaluation_metrics = evaluation_metrics)

model.filename <- paste0("trainModel_", study.name, ".rda")

save(local_model, file = file.path(model_out, model.filename))
