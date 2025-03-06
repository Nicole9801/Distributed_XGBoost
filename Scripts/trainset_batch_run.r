data.path <- "../data/"
output.path <- "../Output/"
script <- "Create_train_set.r"

# list all .rda files in data directory
file_list <- list.files(path = data.path, pattern = "\\.rda", full.names = TRUE)
cat("Found files:\n", paste(file_list, collapse = "\n"), "\n")

#loop through each file and call Create_train_set.r script
for (file in file_list) {
  cat("Processing file: ", file, "\n")
  system(paste("Rscript", script, shQuote(file), shQuote(output.path)))
}

cat("Training set creation complete\n")

