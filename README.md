data/: all original dataset for each study in .rda files 

Training_set/: training set created for each study, including all selected signature scores and response columns, stored in .csv and .rda files

Scripts/: R scripts for creating training set, traing model and distributed ML

Model/: stores trained xgboost model for each study

training_set_summary.csv: a summary including sample size (patients with expression and response data), number of selected signatures and number of missing signatures in the training set
