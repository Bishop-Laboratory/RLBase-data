#' This script uses RLSeq to classify R-loop-mapping samples

library(RLSeq)
library(tidyverse)

# Load the rlfsRes
load("rlbase-data/misc/rlfsRes.rda")

# Load the difficult-to-classify ones
load("misc-data/todiscard.rda")

# Get the predictions
pred <- lapply(rlfsRes, predictCondition)

# Extract label predicted
verd <- sapply(pred, purrr::pluck, "Verdict")

# Combine results with sample data
manifest <- readr::read_tsv("rlbase-data/rlbase_manifest_final.tsv")
manifest <- left_join(manifest, 
                      tibble(
                        experiment = names(verd),
                        verdict = verd
                      ))

manifest %>%
  filter(! is.na(verdict)) %>%
  select(experiment, mode, condition, lab, tissue, verdict) %>%
  mutate(outlier = experiment %in% discard) %>%
  View()

