# Script to calculate gene signature scores, find common features across study dataset, create training set for each study dataset 

library(PredictioR)
library(MultiAssayExperiment)
library(reticulate)
library(dplyr)
library(GSVA)
library(pROC)
library(stringr)
library(readr)

source("Scripts/Compute_GeneSigScore.r")

data.path <- "data/"
file_list <- list.files(path = data.path, pattern = "\\.rda$", full.names = TRUE)
dir_out <- "Training_set/"


Prefiltered_TrainSet <- list()
study.icb <- list()
# Loop through each file to Calculate gene signature score 
for (file in file_list) {
    load(file)
    patient_response <- data.frame(
        patientid = colData(dat$ICB)$patientid,
        response = colData(dat$ICB)$response)
    patient_response <- na.omit(patient_response)
    valid_patients <- patient_response$patientid
    # Skip if valid patient count is less than 25
    if (length(valid_patients) < 25) {
        next
    }

    # Keep only the samples (columns) that have valid response data
    dat$ICB <- dat$ICB[, colnames(assay(dat$ICB)) %in% valid_patients]
    # Update colData to remove NA response patients
    colData(dat$ICB) <- colData(dat$ICB)[colData(dat$ICB)$patientid %in% valid_patients, ]

    # Compute gene signature score
    expr <- dat$ICB
    signature <- dat$signature
    signature_info <- dat$sig.info
    geneSig.score <- compute_gene_signature_scores(expr, signature, signature_info, study.icb)

    ## COMMENT: selected_signature hasn't loaded or defined. 
    trainVar <- geneSig.score[rownames(geneSig.score) %in% selected_signature, , drop = FALSE]

    # Add patient response data 
    trainClin <- patient_response
    trainClin$response <- factor(trainClin$response, levels = c("NR", "R"))

    # transpose and add response column 
    trainVar <- t(trainVar)
    trainVar_df <- as.data.frame(trainVar)
    trainVar_df$patientid <- rownames(trainVar_df)
    TrainSet <- merge(trainVar_df, trainClin[, c("patientid", "response")], by = "patientid")
    TrainSet$patientid <- NULL

    #add each gene signature score and response data to the list 
    Prefiltered_TrainSet[[file]] <- TrainSet
    study.icb[[file]] <- colnames(TrainSet)
  }



# Find common gene signature across selected datasets
common_var <- Reduce(intersect, study.icb)
#writeLines(common_var, "common_features.txt")

# Create and save training set
for (file in names(Prefiltered_TrainSet)) {
  filtered_scores <- Prefiltered_TrainSet[[file]] %>% select(all_of(common_var))
  output_file <- file.path(dir_out, paste0("train_set_", tools::file_path_sans_ext(basename(file)), "_filtered.csv"))
  write_csv(filtered_scores, output_file)
  save(train_dat, file=file.path(dir_out, paste0("train_set_", study.name, ".rda")))
} 


########## CHECK ##########
# Check the column names for each training set file and see they have the same gene signature#
#####

train_set <- list.files(path = dir_out, pattern = "\\.csv$", full.names = TRUE)

all_colnames <- lapply(train_set, function(file) {
  data <- read.csv(file)
  num_col = ncol(data)
  column_names = colnames(data)
})
