#' Script which builds the model

library(tidyverse)

# Magic
RLFSRDAFILE <- "rlbase-data/misc/rlfsRes.rda"
MANIFEST_FINAL <- "rlbase-data/rlpipes-out/config.tsv"
TODISCARD <- 'misc-data/todiscard.rda'
RLFSDIR <- "rlbase-data/rlpipes-out/rlfs_rda/"

# Get samples to discard
load(TODISCARD)
load(RLFSRDAFILE)

# Get the manifest
manifest <- read_tsv(MANIFEST_FINAL)
manifestModel <- manifest %>% 
  filter(experiment %in% names(rlfsRes),
         ! duplicated(experiment),
         ! experiment %in% discard) %>%
  arrange(condType, experiment) %>%
  mutate(
    group = case_when(
      condType == "POS" ~ "CASE",
      TRUE ~ "CTRL"
    ),
    filename = paste0(experiment, "_", genome, ".rlfs.rda")
  ) %>%
  select(id=experiment, group, filename) %>% unique()

write_csv(manifestModel, file = file.path(RLFSDIR, "manifest.csv"))

# Make the model!
rmarkdown::render(input = "scripts/FFT-clasiffier.Rmd", output_file = "../misc-data/FFT-classifier.html")

