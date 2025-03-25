# Method Description
Multivariable model was developed with a distributed pipeline. Eleven datasets, each containing more than 25 patient samples and representing various cancer types including melanoma, lung, gastric, kidney, bladder were selected for analysis.  Distributed model training was conducted using Apache Spark, a cluster computing framework that distributes data across the cluster and enables parallelized operations, and SparkR, an R API for Spark. Independent XGBoost models were trained in parallel for each dataset without data pooling. Hyperparameter tuning for each model was performed using a grid search approach.Hyperparameter tuning for each local model was performed through grid search. *All data from each dataset was used for training due to small sample size*. 

A tree-based aggregation approach was used to collect and integrate independently trained XGBoost models into a single global model. The aggregation processes followed a tree-bagging strategy, where individual decision trees from each local model were sequentially merged into the global model without modifying the learned structures, ensuing efficient model integration while maintaining local dataset-specific patterns. 

The global model was then validated on four private dataset. The validation datasets were selected based on their diverse composition including pan-cancer and single cancer samples, to ensure a robust assessment of the global model's predictive capacity across multiple cancer types. Receiver Operating Characteristic (ROC) curve analysis with 95% confidence intervals (CIs) was performed using the ‘pROC’ R package (v1.18.5).

## Directories and Scripts

### data/
contains 11 original dataset for each study in .rda files 

### Training_set/
training set created for each study, including all selected signature scores and response columns, stored in .csv and .rda files

### Scripts/
Code scripts for

1 - creating training set 

2 - Train_Distributed_XGBoost.r : traing model with distributed ML

3 - Aggregate models into global model (.py)

4 - validation of global model. 

### Local_model/
stores trained xgboost model for each study

### Global_model/
stores aggregated model in different formats

### Validation/
Saves the raw test dataset, signature scores and result figures 

### common_feature.txt
53 signatures that are shared among 11 dataset

### selected_signatures.txt
62 signatures that are most associated with IO responses

### training_set_summary.csv
a summary including sample size (patients with expression and response data), number of selected signatures and number of missing signatures in the training set

## To Run Distributed Training
First, download spark-3.2.1-bin-hadoop3.2.tgz from the Apache Spark Archive: https://archive.apache.org/dist/spark/spark-3.2.1/
and install it on your local path. 

In scripts/2-Train_Distributed_XGBoost.r, modify the path of Spark to the file path on your local computer: 

`Sys.setenv(SPARK_HOME = "/Users/nicole/spark/spark-3.2.1-bin-hadoop3.2/")`

`.libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))`

Then in your R script: 

`install.packages("SparkR")`

This pipeline uses:

R version 4.4.1 
Python version 3.8.19



