#' Script which prepares peaks for R-loop consensus building workflow
#' It randomly subsets to a NUM_SELECT to ensure same # peaks 
#' are provided by each samples. It also used the ZData prediction
#' to ensure samples are of good quality. 

library(tidyverse)
library(magrittr)
library(pbapply)
pbo <- pboptions(type="txt") 

# Magic
MANIFEST <- "rlbase-data/rlbase_samples.tsv"
PEAKSORIG <- "rlbase-data/rlpipes-out/peaks/"
NUM_SELECT <- 5000
PCUTOFF <- 1.3

# New dirs for this analysis
RLCONDIR <- "rlbase-data/rlregions"
RL38DIR <- file.path(RLCONDIR, "peaksHG38/")
RLFINDIR <- file.path(RLCONDIR, "peaksFinal/")

# Make/clean dir tree
dir.create(RLCONDIR, showWarnings = FALSE)
system(paste0("rm -rf ", RLFINDIR)); dir.create(RLFINDIR, showWarnings = TRUE)
system(paste0("rm -rf ", RL38DIR)); dir.create(RL38DIR, showWarnings = TRUE)

# Get the manifest
manifest <- read_tsv(MANIFEST, show_col_types = FALSE) %>%
  unique()

# Get the number of peaks meeting a cutoff of p<1E5
if ("numPeaks" %in% colnames(manifest)) {
  manifest <- select(manifest, -numPeaks)
}
peaks <- list.files(PEAKSORIG, pattern = "\\.broadPeak", full.names = TRUE)
peaks <- tibble(
  experiment = gsub(peaks, pattern = ".+/([ES]{1}RX[0-9]+)_([a-zA-Z0-9]+)\\.broadPeak", replacement = "\\1"),
  genome = gsub(peaks, pattern = ".+/([ES]{1}RX[0-9]+)_([a-zA-Z0-9]+)\\.broadPeak", replacement = "\\2"),
  peakfile = peaks
)
# Get the number of peaks below the pval cutoff
lens <- pbapply::pbsapply(peaks$peakfile, function(x) {
  system(paste0("awk '$9>", PCUTOFF, "' ", x, " | wc -l "), intern = TRUE)
}) %>% as.numeric()
nPeaks <- peaks %>%
  select(experiment) %>%
  mutate(numPeaks = {{ lens }})
manifest <- left_join(manifest, nPeaks, by = "experiment") 
if ("run" %in% colnames(manifest)) {
  manifest <- unique(select(manifest, -run))
}
# Force unique experiment
manifest <- manifest %>%
  distinct(experiment, .keep_all = TRUE) 
write_tsv(manifest, file = MANIFEST)


# Get the final peakset for R-loop consensus
samplesForConsensus <- manifest %>%
  filter(
    genome == "hg38",
    group == "rl",
    condType == "POS",
    verdict == "Case",
    numPeaks > NUM_SELECT
  )
write_tsv(samplesForConsensus, file = "misc-data/samplesForConsensus.tsv")

# Move the peaks 
consSamp <- samplesForConsensus %>%
  dplyr::mutate(peak = paste0(PEAKSORIG, experiment, "_hg38.broadPeak")) %>%
  dplyr::filter(file.exists(peak)) %T>% {
    pull(., experiment) %>%
      pbsapply(function(x) {
        system(paste0("awk '$9>", PCUTOFF, "' ",paste0(PEAKSORIG, x, "_hg38.broadPeak"), 
                      ' | shuf -n ', NUM_SELECT, 
                      " > ", RLFINDIR, x, "_hg38.broadPeak"))
      })
  }

# Get the HG38 peaks (will be used for final RLRegions Table metadata)
hg38Mani <- manifest %>%
  dplyr::filter(
    genome == "hg38",
    condType == "POS",
    numPeaks > 0
  ) %>%
  dplyr::mutate(peak = paste0(PEAKSORIG, experiment, "_hg38.broadPeak")) %>%
  dplyr::filter(file.exists(peak)) %T>% {
    pull(., experiment) %>%
      pbsapply(function(x) {
        system(paste0('cp ', paste0(PEAKSORIG, x, "_hg38.broadPeak "), 
                      RL38DIR, x, "_hg38.broadPeak"))
      })
  }

# Write manifests for snakemake
consSamp %>% 
  select(peak) %>%
  mutate(peak = normalizePath(peak)) %>%
  unique() %>%
  write_csv("rlbase-data/rlregions/rlregion_manifest_All.csv", col_names = FALSE)
hg38Mani %>% 
  select(peak) %>%
  mutate(peak = normalizePath(peak)) %>%
  unique() %>%
  write_csv("rlbase-data/rlregions/hg38_manifest_All.csv", col_names = FALSE)

# Write manifests for S96 and RNH
consSamp %>%
  filter(ip_type == "S9.6") %>%
  select(peak) %>%
  mutate(peak = normalizePath(peak)) %>%
  unique() %>%
  write_csv("rlbase-data/rlregions/rlregion_manifest_S96.csv", col_names = FALSE)
consSamp %>%
  filter(ip_type == "dRNH") %>%
  select(peak) %>%
  mutate(peak = normalizePath(peak)) %>%
  unique() %>%
  write_csv("rlbase-data/rlregions/rlregion_manifest_dRNH.csv", col_names = FALSE)
hg38Mani %>% 
  filter(ip_type == "S9.6") %>%
  select(peak) %>%
  mutate(peak = normalizePath(peak)) %>%
  unique() %>%
  write_csv("rlbase-data/rlregions/hg38_manifest_S96.csv", col_names = FALSE)
hg38Mani %>% 
  filter(ip_type == "dRNH") %>%
  select(peak) %>%
  mutate(peak = normalizePath(peak)) %>%
  unique() %>%
  write_csv("rlbase-data/rlregions/hg38_manifest_dRNH.csv", col_names = FALSE)

