
library(sparklyr)

library(dplyr)
library(xgboost)
library(parallel)

Sys.setenv(PATH = paste("/usr/local/bin", Sys.getenv("PATH"), sep = ":"))

# Connect to spark 
sc <- spark_connect(master = "local", version = "3.2.1") # successful connection. Modify version if necessary


# set data path and model path
data.path <- "Output/"
data.files <- list.files(path = data.path, pattern = "*.rda", full.names = TRUE)   
model_path <- normalizePath("Scripts/train_local_model.r")

rscript_path <- Sys.which("Rscript")

#rscript_path <- "/usr/local/bin/Rscript"

#convert data into spark dataframe
data_tbl <- data.frame(dataset_path = data.files) %>%
    sdf_copy_to(sc, ., overwrite = TRUE)


# function to distibute training for each dataset on Spark
distribute_train <- function(data_tbl, model_path, rscript_path){
    # apply training function with Spark
    data_tbl %>%
        spark_apply(
            function(df, context){
            data_path <- df$dataset_path
            output_path <- Sys.getenv("SPARK_LOCAL_DIRS", unset = "Model_2")

            temp_data_path <- file.path(tempdir(), basename(data_path))
            file.copy(data_path, temp_data_path, overwrite = TRUE)

            # Copy the training script to a temporary location on the worker node
            temp_model_path <- file.path(tempdir(), "train_local_model.r")
            file.copy(context$script_path, temp_model_path, overwrite = TRUE)

            #record start time for each training set
            start_time <- Sys.time()
            message(paste("Starting training for:", data_path, "at", start_time))
            if (!file.exists(temp_model_path)) stop("train_local_model.r not found")

            # Run the training script
            command <- paste(shQuote(rscript_path),
                             shQuote(temp_model_path), 
                             shQuote(temp_data_path),
                             shQuote(context$output_path))
    
            output <- system(command)
            message("Command Output for ", data_path, ": ", paste(output, collapse = "\n"))

            # record end time for this training set
            end_time <- Sys.time()
            message(paste("Completed training for:", data_path, "at", end_time))
            duration <- difftime(end_time, start_time, units = "secs")

            return(data.frame(status = paste("Processed: ", data_path),time_taken = as.numeric(duration)))
        },
        context = list(script_path = model_path,
                       rscript_path = rscript_path,
                       output_path = output_path),

        columns = list(status = "character", time_taken = "numeric") # define output column
    ) %>% collect()
}

global_start_time <- Sys.time()
message("Starting all training jobs at:", global_start_time)


distribute_train(data_tbl, model_path, rscript_path)

model_output <- tryCatch(
  distribute_train(data_tbl, model_path, rscript_path),
  error = function(e) {
    message("Training failed: ", e$message)
    return(NULL)
  }
)


global_end_time <- Sys.time()
message("All training jobs completed at:", global_end_time)
total_duration <- difftime(global_end_time, global_start_time, units = "secs")
message("Total execution time:", total_duration, "seconds")


spark_disconnect(sc)


# function to aggregate local models using federated learning
aggregate_models <- function(local_models) { 
    num_models <- length(local_models)
    avg_model <- local_models[[1]]
    for (i in 2:num_models) {
        model <- readRDS(local_models[[i]])
        avg_model <- avg_model + model
    }
}


if (!dir.exists(model_out)) {
  message("⚠️ Model output directory does not exist. Creating: ", model_out)
  dir.create(model_out, recursive = TRUE)
}
