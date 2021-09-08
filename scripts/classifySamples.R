#' This script uses RLSeq to classify R-loop-mapping samples
library(RLSeq)
library(tidyverse)
library(pbapply)
pbo <- pboptions(type="txt") 

# Magic
MANIFEST <- "rlbase-data/rlpipes-out/config.tsv"
MANIFEST_FINAL <- "rlbase-data/rlbase_samples.tsv"
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
  mutate(discarded = experiment %in% discard) %>%
  filter((mode == "RNA-Seq" & is.na(verdict)) | mode != "RNA-Seq")

readr::write_tsv(manifest, MANIFEST_FINAL)
