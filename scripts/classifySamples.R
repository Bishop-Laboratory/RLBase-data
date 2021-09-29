#' This script uses RLSeq to classify R-loop-mapping samples
library(RLSeq)
library(tidyverse)
library(pbapply)
pbo <- pboptions(type="txt") 

# Magic
MANIFEST <- "rlbase-data/rlpipes-out/config.tsv"
MANIFEST_FINAL <- "rlbase-data/rlbase_samples.tsv"
RLFSRDA <- "rlbase-data/misc/rlfsRes.rda"
TODISCARD <- "misc-data/model/todiscard.rda"
PREPMODEL <- "misc-data/model/prepFeatures.rda"
FFTMODEL <- "misc-data/model/fftModel.rda"

# Load the rlfsRes
load(RLFSRDA)

# Load the difficult-to-classify ones
load(TODISCARD)

# Load all the models
load(PREPMODEL)
load(FFTMODEL)

# Get the predictions
pred <- pbapply::pblapply(rlfsRes, function(x) {
  predictCondition(
    rlfsRes = x,
    prepFeatures=prepFeatures,
    fftModel=fftModel
  )
})

# Extract label predicted
verd <- sapply(pred, purrr::pluck, "prediction")

# Combine results with sample data
manifest <- readr::read_tsv(MANIFEST, show_col_types = FALSE)
if ("prediction" %in% colnames(manifest)) {
  manifest <- select(manifest, -prediction)
}
manifest <- left_join(manifest, 
                      tibble(
                        experiment = names(verd),
                        prediction = verd
                      )) %>%
  mutate(discarded = experiment %in% {{ discard }}) %>%
  filter((mode == "RNA-Seq" & is.na(prediction)) | mode != "RNA-Seq")

readr::write_tsv(manifest, MANIFEST_FINAL)

message("Done")
