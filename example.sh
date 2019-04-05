#!/usr/local/bin/Rscript

# generate dataset with certain seed
set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/raceid_stemid",
  num_cells = 99,
  num_features = 101,
  model = "tree",
  normalise = FALSE
)

# add method specific args (if needed)
data$parameters <- list(pdishuf = 50)

data$seed <- 1L

# write example dataset to file
file <- commandArgs(trailingOnly = TRUE)[[1]]
dynutils::write_h5(data, file)
