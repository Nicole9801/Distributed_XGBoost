
# Script to clean up data, calculate gene signature scores and create training set for each study dataset 

library(PredictioR)
library(MultiAssayExperiment)
library(reticulate)
library(dplyr)
library(GSVA)
library(pROC)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) ==0) {
  stop("Error: Missing arguments. Please provide the path to the data file and the study name.")
}

data.path <- "../data/"
data.name <- args[1]  # Get the file name from command line arguments
study.name <- tools::file_path_sans_ext(basename(data.name))

# Load the .rda file
load(file.path(data.path, data.name))

if (!exists("dat")) {
  stop("Error: 'dat' object not found in the loaded file.")
}

dir_out <- "../Output/"

# extract the clinical response data in the SummarizedExperiment object
selected_signature <- c(
  "TMEscoreA_plus", "APM_Thompson", "CD39-CD8Tcell_Chow", "CD8_Rooney",
  "CYT_Rooney", "Cytotoxic-cells_Danaher", "PassON_Du", "Teff-IFNG_Fehrenbacher",
  "Tcells_Danaher", "CCL5-CXCL9_Dangaj", "Co-Inhibition-Tcell_Rooney", "CYT_Davoli",
  "Tcells_Bindea", "VIGex_Hernando-Calvo", "Inflammatory_Thompson", "CD8_Tcells_Danaher",
  "Cytotoxic-cells_Bindea", "Tcell-inflamed_Ayers", "ImmunoScore_Hao", "Teff_McDermott",
  "MHC-I_Rooney", "Bcells_Rooney", "MHC-I_Liu", "peri-Tcell_Hwang",
  "STAT1_Care", "M1_Hwang", "PDL1_Nishino", "APM_Wang",
  "ImmuneScore_Roh", "TIS_Damotte", "Response_Chen", "NK-cells_Danaher",
  "MHC-II_Liu", "IFN_Ayers", "Co-Inhibition-APC_Rooney", "Co-Stimulation-Tcell_Rooney",
  "NK-CD56-cells_Danaher", "IRG_Ayers", "Bcell_Budczies", "Chemokine_Messina",
  "Macrophage_Danaher", "TypeI-IFNReponse_Rooney", "Bcells_Danaher", "T-helper-cells_Bindea",
  "TLS_Cabrita", "TMEscoreB_plus", "ICB-resistance_Peng", "NK-cells_Rooney",
  "DNA_replication", "Homologous_recombination", "CIN25_Carter", "Bcell_Helmink",
  "CellCycle", "EMT2", "EMT3", "Base_excision_repair",
  "CIN70_Carter", "IMPRES_Auslander", "SW480-cells_Bindea", "IPSOV_Shen",
  "FGFR3_related", "proliferation_Pabla"
)

patient_response <- data.frame(
  patientid = colData(dat$ICB)$patientid,
  response = colData(dat$ICB)$response
)

patient_response <- na.omit(patient_response)

valid_patients <- patient_response$patientid

# Keep only the samples (columns) that have valid response data
dat$ICB <- dat$ICB[, colnames(assay(dat$ICB)) %in% valid_patients]
# Update colData to remove NA response patients
colData(dat$ICB) <- colData(dat$ICB)[colData(dat$ICB)$patientid %in% valid_patients, ]


expr <- dat$ICB
signature <- dat$signature
signature_info <- dat$sig.info
trainClin <- patient_response

# Compute gene signature scores with existing pipeline
geneSig.score <- lapply(1:length(signature), function(i){
  
  print(paste(i , names(signature)[i], sep="/"))
  sig_name <- names(signature)[i]
  
  # check which method is associated with the signature in signature_info and applies the corresponding function.
  if(signature_info[signature_info$signature == sig_name, "method"] == "GSVA"){
    
    geneSig <- geneSigGSVA(dat.icb = expr,
                           sig = signature[[i]],
                           sig.name = sig_name,
                           missing.perc = 0.5,
                           const.int = 0.001,
                           n.cutoff = 15,
                           sig.perc = 0.8,
                           study = study.icb)
    
    
    if(sum(!is.na(geneSig)) > 0){
      geneSig <- geneSig[1,]
    }     
    
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Weighted Mean"){
    
    geneSig <- geneSigMean(dat.icb = expr,
                           sig = signature[[i]],
                           sig.name = sig_name,
                           missing.perc = 0.5,
                           const.int = 0.001,
                           n.cutoff = 15,
                           sig.perc = 0.8,
                           study = study.icb)
    
  }
  
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "ssGSEA"){
    
    geneSig <- geneSigssGSEA(dat.icb = expr,
                             sig = signature[[i]],
                             sig.name = sig_name,
                             missing.perc = 0.5,
                             const.int = 0.001,
                             n.cutoff = 15,
                             sig.perc = 0.8,
                             study = study.icb)
    
    if(sum(!is.na(geneSig)) > 0){
      geneSig <- geneSig[1,]
    }     
    
    
  }
  
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "COX-IS_Bonavita"){
    
    geneSig <- geneSigCOX_IS(dat.icb = expr,
                             sig = signature[[i]],
                             sig.name = sig_name,
                             missing.perc = 0.5,
                             const.int = 0.001,
                             n.cutoff = 15,
                             sig.perc = 0.8,
                             study = study.icb)
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPS_Charoentong"){
    
    geneSig <- geneSigIPS(dat.icb = expr,
                          sig = signature[[i]],
                          sig.name = sig_name,
                          missing.perc = 0.5,
                          const.int = 0.001,
                          n.cutoff = 15,
                          study = study.icb)
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PredictIO_Bareche"){
    
    geneSig <- geneSigPredictIO(dat.icb = expr,
                                sig = signature[[i]],
                                sig.name = sig_name,
                                missing.perc = 0.5,
                                const.int = 0.001,
                                n.cutoff = 15,
                                sig.perc = 0.8,
                                study = study.icb)
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPRES_Hugo"){
    
    geneSig <- geneSigIPRES(dat.icb = expr,
                            sig = signature[[i]],
                            sig.name = sig_name,
                            missing.perc = 0.5,
                            const.int = 0.001,
                            n.cutoff = 15,
                            sig.perc = 0.8,
                            study = study.icb)
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PassON_Du"){
    
    geneSig <- geneSigPassON(dat.icb = expr,
                             sig = signature[[i]],
                             sig.name = sig_name,
                             missing.perc = 0.5,
                             const.int = 0.001,
                             n.cutoff = 15,
                             sig.perc = 0.8,
                             study = study.icb)
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPSOV_Shen"){
    
    geneSig <- geneSigIPSOV(dat.icb = expr,
                            sig = signature[[i]],
                            sig.name = sig_name,
                            missing.perc = 0.5,
                            const.int = 0.001,
                            n.cutoff = 15,
                            sig.perc = 0.8,
                            study = study.icb)
    
  }
  
  
  if(sum(!is.na(geneSig)) > 0){
    
    geneSig <- geneSig
    
  }     
  
  if(sum(!is.na(geneSig)) == 0){
    
    geneSig <- rep(NA, ncol(expr))
    
  }
  
  geneSig
  
})

# bind gene signature scores and remove rows with NA values
geneSig.score <- do.call(rbind, geneSig.score)
rownames(geneSig.score) <- names(signature)
remove <- which(is.na(rowSums(geneSig.score)))
if(length(remove) > 0){
  
  geneSig.score <- geneSig.score[-remove, ]
  
}

###############################
# Create training set
###############################

train_data <- geneSig.score[rownames(geneSig.score) %in% selected_signature, ]
train_data <- t(train_data)
trainClin$response <- factor(trainClin$response, levels = c("NR", "R"))

# merge clinical response with gene signature scores to create complete training set
train_data_df <- as.data.frame(train_data)
train_data_df$patientid <- rownames(train_data_df)
full_train_set <- merge(train_data_df, trainClin[, c("patientid", "response")], by = "patientid")
full_train_set$patientid <- NULL
train_dat <- full_train_set

save(train_dat, file=file.path(dir_out, paste0("train_set_", study.name, ".rda")))
write.csv(train_dat, file.path(dir_out, paste0("train_set_", study.name, ".csv")), row.names = FALSE)


data.path <- "Output/"
file_list <- list.files(path = data.path, pattern = "\\.csv$", full.names = TRUE)

results <- list()
for (file in file_list) {
  df <- read.csv(file)
  missing_signatures <- setdiff(sort(selected_signature), sort(colnames(df)))
  results[[basename(file)]] <- list(
    num_signature = ncol(df)-1,
    missing_signatures = paste(missing_signatures, collapse = ",")
  )
  results_df <- do.call(rbind, lapply(names(results), function(x) {
    c(File = x, results[[x]])
  }))
  results_df <- as.data.frame(results_df)
}

results_df <- data.frame(
  num_signature = as.numeric(unlist(lapply(results, `[[`, "num_signature"))),
  missing_signatures = sapply(results, function(x) paste(x$missing_signatures, collapse = ", ")),
  row.names = names(results),
  stringsAsFactors = FALSE
)

write.csv(results_df, file.path("Output/", "train_set_summary.csv"), row.names = TRUE)

print(results_df)


df <- read.csv("Output/train_set_ICB_Van_Allen__Melanoma__CTLA4.csv")
missing_signatures <- setdiff(selected_signature, colnames(df))
