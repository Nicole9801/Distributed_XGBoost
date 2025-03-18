
# Script to test the global model on four additional dataset, compare with PredictIO and model trained with centralized data

library(data.table)
library(xgboost)
library(dplyr)
library(ggplot2)
library(readr)

source("Scripts/Compute_GeneSigScore.r")


val_dir <- "Validation_data"

# Load the global model 
global_model <- xgb.load("Global_model/global_xgb_model.model")

#################################################
## load validation data and signature data
#################################################
load(file.path("Validation_data/testDataIMmotion150.rda"))
selected_signature <- scan("common_features.txt", what = character(), sep = "\n")

testClin <- data.frame(colData(dat$testData))
testClin <- testClin[!is.na(testClin$response), ]
group <- ifelse(testClin$response == "R", 0, 1) 
expr <- assay(dat$testData)
expr <- expr[, testClin$patientid]

signature <- dat$signatureData$signature
signature_info <- dat$signatureData$sig.info

sig <- dat$PredictIOSig

##############################################
## compute signature scores
##############################################
study.icb <- "IMmotion150"

geneSig.score <- compute_gene_signature_scores(expr, signature, signature_info, study.icb)

save(geneSig.score, file=file.path(dir, paste(study.icb, "sig_score.rda", sep="_")))


load("Validation_data/IMmotion_sig_score.rda")


############################################
## Processing signature score data
############################################

testExpr <- t(geneSig.score)
#testExpr <- testExpr[testClin$patientid, ]
testExpr <- testExpr[, colnames(testExpr) %in% selected_signature]

missing_features <- setdiff(selected_signature, colnames(testExpr))

if(length(missing_features) != 0){
  
  # Add missing features to testExpr and set values to 0
  testExpr_missing <- lapply(1:length(missing_features), function(k){
    
    rep(0, nrow(testExpr))
    
  })
  
  testExpr_missing <- do.call(cbind, testExpr_missing)
  colnames(testExpr_missing) <- missing_features
  rownames(testExpr_missing) <- rownames(testExpr)
  testExpr <- cbind(testExpr, testExpr_missing)
  
}

# Ensure the order is correct
testExpr <- testExpr[, selected_signature]


###########################################################
## Model1: Validate MV model prediction
###########################################################
predMV <- predict(global_model, newdata = testExpr)
rocMV <- roc(testClin$response, predMV)
ciMV <- ci.auc(rocMV, method = "bootstrap", boot.n = 1000)

###########################################################
## Model2: Validate PredictIO signature prediction
###########################################################
# compute PredictIO signature score
geneSig <- geneSigPredictIO(dat.icb = expr,
                            sig = sig,
                            sig.name = "PredictIO",
                            missing.perc = 0.5,
                            const.int = 0.001,
                            n.cutoff = 15,
                            sig.perc = 0.8,
                            study = study.icb)

group <- ifelse(testClin$response %in% "R", 0, ifelse(testClin$response %in% "NR", 1, NA))
fit <- glm(group ~ geneSig, family = binomial(link = "logit"))
predPredictIO <- predict(fit, type = "response")
rocPredictIO <- roc(group, predPredictIO)
ciPredictIO <- ci.auc(rocPredictIO, method = "bootstrap", boot.n = 1000)


###########################################################
## Model3: Centralized model prediction
###########################################################


###########################################################
## Model4: 11 local trained models: run each model on validation dataset
###########################################################
model_files <- list.files(path = "Local_model/", pattern = "^XGB.*\\.rds$", full.names = TRUE)

roc_results <- list()

for (i in seq_along(model_files)) {
  # Load the trained XGBoost model
  local_model <- readRDS(model_files[i])
  # Predict with the model
  predMV <- predict(local_model, newdata = testExpr)  # Ensure testExpr is a matrix
  rocMV <- roc(testClin$response, predMV)
  roc_results[[i]] <- rocMV
  print(paste("Local Model:", as.character(i), "AUC:", auc(rocMV)))
}


###########################################################
## visualization
###########################################################

roc_obj <- list("XGBoost" = rocMV, "PredictIO" = rocPredictIO)
roc_obj_withLocal <- c(roc_obj, roc_results) # a list of ROC containing 11 local model predictions on validation data
pred <- list('predPredictIO' = predPredictIO, 'predMV' = predMV)
res <- list('roc_obj' = roc_obj, 'prediction' = pred)

#save(res, file = file.path(dir, paste(study.icb, "validation.RData", sep="_")))

# Create the plot with ggplot2
jpeg(file = file.path(val_dir, paste(study.icb, "roc_curve.jpeg", sep="_")), width = 650, height = 650, res = 150)

# Figure of global xgboost model vs. predictIO
ggroc(roc_obj, linewidth = 1.3, legacy.axes = TRUE) +
  scale_colour_manual(values = c("#733a4d", "#204035FF", "black")) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color="#8E949FFF", linetype="dashed", linewidth = 1.3) +
  xlab("Specificity") + 
  ylab("Sensitivity") +
  ggtitle("IMmotion150: Kidney") +
  theme(axis.text.x = element_text(size = 10, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), plot.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.position = "none", 
        #legend.text = element_text(size = 7, face = "bold"), 
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  annotate("text", x = 0.7, y = 0.15, size = 3.5, label = paste(paste("PredictIO AUC:", round(auc(roc_obj$PredictIO), 2), sep=""), 
                                                              paste("[", round(ciPredictIO[1], 2), "-", round(ciPredictIO[2], 2), "]", sep=""), sep=""), color = "#204035FF") +
  annotate("text", x = 0.7, y = 0.1, size = 3.5, label = paste(paste("XGBoost AUC:", round(auc(roc_obj$XGBoost), 2), sep=""), 
                                                             paste("[", round(ciMV[1], 2), "-", round(ciMV[2], 2), "]", sep=""), sep=""), color = "#733a4d")


# Figure of decentralized vs. centralized model 

dev.off()




