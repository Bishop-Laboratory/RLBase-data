#' This script uses RLSeq to classify R-loop-mapping samples

library(RLSeq)
library(tidyverse)

# Magic
MANIFEST <- "rlbase-data/rlbase_manifest_final.tsv"
RLFSRDA <- "rlbase-data/misc/rlfsRes.rda"
TODISCARD <- "misc-data/todiscard.rda"

# Load the rlfsRes
load(RLFSRDA)

# Load the difficult-to-classify ones
load(TODISCARD)

# Get the predictions
pred <- pbapply::pblapply(rlfsRes, predictCondition)

# Extract label predicted
verd <- sapply(pred, purrr::pluck, "Verdict")

# Combine results with sample data
manifest <- readr::read_tsv(MANIFEST)
if ("verdict" %in% colnames(manifest)) {
  manifest <- select(manifest, -verdict)
}
manifest <- left_join(manifest, 
                      tibble(
                        experiment = names(verd),
                        verdict = verd
                      )) %>%
  mutate(discarded = experiment %in% discard)

readr::write_tsv(manifest, MANIFEST)
