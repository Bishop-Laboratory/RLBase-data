#' Cleans up old datasets which need to be re-run

suppressPackageStartupMessages(library(tidyverse))

# Get args
args <- commandArgs(trailingOnly = TRUE)

# Location of the catalog excel file
CATALOG=args[1]
message(CATALOG)
CONFIG=args[2]
message(CONFIG)

# get the config
config <- read_tsv(CONFIG)

# Get the catalog, condition map, and modes
redo <- readxl::read_excel(CATALOG, sheet = "redo") %>% pull(toRedo)

# Change config genome
srxToRedo <- config %>%
  filter(experiment_original %in% redo | experiment %in% redo) %>%
  pull(experiment)

if (length(srxToRedo) > 0) {
  # Delete all entries pertaining to these samples
  fls <- list.files("rlbase-data/rlpipes-out/", full.names = TRUE, recursive = TRUE, pattern = srxToRedo)
  sapply(fls, file.remove)
}


# Echo out
message("\nDone\n")
