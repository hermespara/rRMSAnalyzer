# Manually load data for local testing (since package not installed)
# Try to find the file relative to where tests run (usually tests/testthat)

data_path <- "../../data/ribo_toy.rda"
if (file.exists(data_path)) {
  load(data_path, envir = .GlobalEnv)
  message("Loaded ribo_toy from ", data_path)
} else {
  # Try from project root
  data_path <- "data/ribo_toy.rda"
  if (file.exists(data_path)) {
     load(data_path, envir = .GlobalEnv)
     message("Loaded ribo_toy from ", data_path)
  } else {
     warning("Could not find ribo_toy.rda")
  }
}
