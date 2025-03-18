
data/: contains 11 original dataset for each study in .rda files 

Training_set/: training set created for each study, including all selected signature scores and response columns, stored in .csv and .rda files

Scripts/: R scripts for 1-creating training set, 2-traing model with distributed ML, and 3-Aggregate models into global model (.py) and 4 - validation of global model. 

Local_model/: stores trained xgboost model for each study

Global_model/: stores aggregated model in different formats

common_feature.txt: 53 signatures that are shared among 11 dataset

selected_signatures.txt:  62 signatures that are most associated with IO responses

training_set_summary.csv: a summary including sample size (patients with expression and response data), number of selected signatures and number of missing signatures in the training set

