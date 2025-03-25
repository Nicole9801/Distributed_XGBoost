
# Script to test the global model on four additional dataset, compare with PredictIO and model trained with centralized data

library(data.table)
library(xgboost)
library(dplyr)
library(ggplot2)
library(readr)
library(pROC)

source("Scripts/Compute_GeneSigScore.r")

# Load the global model 
global_model <- xgb.load("Global_model/global_xgb_model.model")

#################################################
## load validation data and signature data
#################################################
val_dir <- "Validation"
val_files <- list.files(val_dir, pattern = "\\.rda$", full.names = TRUE)

selected_signature <- scan("common_features.txt", what = character(), sep = "\n")

## compute signature scores for each validation dataset and save the scores
for (file in val_files){
  load(file)
  study.icb <- tools::file_path_sans_ext(basename(file))
  testClin <- data.frame(colData(dat$testData))
  testClin <- testClin[!is.na(testClin$response), ]
  group <- ifelse(testClin$response == "R", 0, 1)
  expr <- assay(dat$testData)
  expr <- expr[, testClin$patientid]

  signature <- dat$signatureData$signature
  signature_info <- dat$signatureData$sig.info

  sig <- dat$PredictIOSig
  geneSig.score <- compute_gene_signature_scores(expr, signature, signature_info, study.icb)

  # Model validation 
  testExpr <- t(geneSig.score)
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

  # Save computed signature scores
  output_dir <- "Validation/Scores/"
  if (!dir.exists(output_dir)) dir.create(output_dir)

  save(testExpr, file = file.path(output_dir, paste0(study.icb, "_sig_score.rda")))

}



############################################
## Load test data - signature score for model validation 
############################################
data_dir <- "Validation/"
score_dir <- "Validation/Scores"

results_dir <- "Validation/Results"
data_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE) # dat


if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}


# Run the global model, predictIO and local models on four validation data
for (data_file in data_files) {
  study.icb <- tools::file_path_sans_ext(basename(data_file))
  print(study.icb)
  score_file <- file.path(score_dir, paste0(study.icb, "_sig_score.rda"))
  load(data_file) 
  load(score_file)

  testClin <- data.frame(colData(dat$testData))
  testClin <- testClin[!is.na(testClin$response), ]
  expr <- assay(dat$testData)
  expr <- expr[, testClin$patientid]
  sig <- dat$PredictIOSig

  # Model 1: Global model
  predMV <- predict(global_model, newdata = testExpr)
  rocMV <- roc(testClin$response, predMV)
  ciMV <- ci.auc(rocMV, method = "bootstrap", boot.n = 1000)

  # Model2: Validate PredictIO signature prediction
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


  # Figure 1 - Generate ROC curve for PredictIO and global model 
  roc_obj <- list("XGBoost" = rocMV, "PredictIO" = rocPredictIO)
  pred <- list('predPredictIO' = predPredictIO, 'predMV' = predMV)
  res <- list('roc_obj' = roc_obj, 'prediction' = pred)

  output_file <- file.path(results_dir, paste0(study.icb, "_roc_curve.jpeg"))

  roc_plot <- ggroc(roc_obj, linewidth = 1.3, legacy.axes = TRUE) +
    scale_colour_manual(values = c("#733a4d", "#204035FF", "black")) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                color="#8E949FFF", linetype="dashed", linewidth = 1.3) +
    xlab("Specificity") + 
    ylab("Sensitivity") +
    ggtitle(paste(study.icb)) +
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
                                                              
    ggsave(filename = output_file, plot = roc_plot, width = 6.5, height = 6.5, dpi = 150)


  ## Model 3: 11 local trained models: run each model on validation dataset
  model_files <- list.files(path = "Local_model/", pattern = "^XGB.*\\.rds$", full.names = TRUE)

  roc_results <- list()

  for (i in seq_along(model_files)) {
    # Load the trained XGBoost model
    local_model <- readRDS(model_files[i])
    # Predict with the model
    predMV <- predict(local_model, newdata = testExpr)  # Ensure testExpr is a matrix
    rocMV_local <- roc(testClin$response, predMV)
    roc_results[[i]] <- rocMV_local
    print(paste("Local Model:", as.character(i), "AUC:", auc(rocMV_local)))
  }

  spec_seq <- seq(0, 1, length.out = 150)
  sens_matrix <- matrix(NA, nrow = length(spec_seq), ncol = length(roc_results))

  for (i in seq_along(roc_results)) {
    roc_obj_i <- roc_results[[i]]
    interp_sens <- coords(roc_obj_i, x = spec_seq, input = "specificity", ret = "sensitivity", transpose = FALSE)[, 1]
    sens_matrix[, i] <- interp_sens
  }

  mean_sens <- rowMeans(sens_matrix, na.rm = TRUE)
  sd_sens <- apply(sens_matrix, 1, sd, na.rm = TRUE)

  auc_local_mean <- mean(sapply(roc_results, auc), na.rm = TRUE)
  ci_local <- quantile(sapply(roc_results, auc), probs = c(0.025, 0.975), na.rm = TRUE)

  auc_xgb <- auc(rocMV)
  auc_predictio <- auc(rocPredictIO)

  label_local <- paste0("Local XGBoost Mean AUC: ", round(auc_local_mean, 2))
  label_xgb <- paste0("Distributed XGBoost AUC: ", round(auc_xgb, 2))
  label_predictio <- paste0("PredictIO AUC: ", round(auc_predictio, 2))

  df_mean_roc <- data.frame(
    specificity = spec_seq,
    sensitivity = mean_sens,
    sens_upper = pmin(mean_sens + sd_sens, 1),
    sens_lower = pmax(mean_sens - sd_sens, 0),
    model = label_local
  )

  df_xgb <- data.frame(
    specificity = rev(rocMV$specificities),
    sensitivity = rev(rocMV$sensitivities),
    model = label_xgb
  )

  df_predictio <- data.frame(
    specificity = rev(rocPredictIO$specificities),
    sensitivity = rev(rocPredictIO$sensitivities),
    model = label_predictio
  )

  df_all_roc <- bind_rows(df_mean_roc, df_xgb, df_predictio)

  # Generate figure 2 - ROC curve with local model, distributed model and PredictIO
  color_map <- setNames(c("#295F85", "#B03C2B", "black"),
                        c(label_local, label_xgb, label_predictio))

  output_file <- file.path(results_dir, paste0(study.icb, "_roc_curve_withLocal.jpeg"))

  roc_plot <- ggplot() +
    geom_ribbon(data = df_mean_roc,
                aes(x = 1 - specificity, ymin = sens_lower, ymax = sens_upper),
                fill = "grey80", alpha = 0.5) +
    geom_line(data = df_all_roc,
              aes(x = 1 - specificity, y = sensitivity, color = model),
              linewidth = 1.3) +
    scale_color_manual(values = color_map) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                color = "#8E949FFF", linetype = "dashed", linewidth = 1.3) +
    xlab("Specificity") +
    ylab("Sensitivity") +
    ggtitle(study.icb) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      panel.background = element_blank(), plot.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) +
    annotate("text", x = 0.7, y = 0.20, size = 3.5,
            label = paste0("Local XGBoost Mean AUC: ", round(auc_local_mean, 2), " [", round(ci_local[1], 2), "–", round(ci_local[2], 2), "]"),
            color = "#295F85") +
    annotate("text", x = 0.7, y = 0.15, size = 3.5,
            label = paste0(label_xgb, " [", round(ciMV[1], 2), "-", round(ciMV[2], 2), "]"),
            color = "#B03C2B") +
    annotate("text", x = 0.7, y = 0.10, size = 3.5,
            label = paste0(label_predictio, " [", round(ciPredictIO[1], 2), "-", round(ciPredictIO[2], 2), "]"),
            color = "black")

  ggsave(filename = output_file, plot = roc_plot, width = 6.5, height = 6.5, dpi = 150)
  }






