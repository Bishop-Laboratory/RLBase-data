library(tidyverse)
library(magrittr)

RMAPFFTSMALL <- "analyses/fft_analysis/rmapfftsmall.rda"
load(RMAPFFTSMALL)

manifest <- read_csv("misc/datasets_for_fft_testing/manifest_full_data.csv")
rmapfftsmall <- left_join(rmapfftsmall, manifest)


# This code find the peaks to include in the RLoop Zone analysis
# It randomly subsets to a NUM_SELECT to ensure same # peaks 
# are provided by each samples. It also used the ZData prediction
# to ensure samples are of good quality. 
RLDIR <- "analyses/Prepare-RMapDB-Tables/RLoops/data/peaksFinal/"
PEAKSDIR <- "analyses/Prepare-RMapDB-Tables/RLoops/data/peaks/"
RL38DIR <- "analyses/Prepare-RMapDB-Tables/RLoops/data/peaksHG38/"
# High value = fewer samples but greater # peaks to analyze and higher precision in results
NUM_SELECT <- 15000
system(paste0("rm -rf ", RLDIR)); dir.create(RLDIR, showWarnings = TRUE)
system(paste0("rm -rf ", RL38DIR)); dir.create(RL38DIR, showWarnings = TRUE)
a_ <- rmapfftsmall %>%
  dplyr::filter(
         genome == "hg38",
         Condition %in% c("S9.6", "D210N", "FLAG", "delta-HC"),
         pval <  .01,
         MACS2__total_peaks > NUM_SELECT,
         `5UTR__LogP enrichment (+values depleted)` < -1.3,
         prediction == "Case"
         ) %>%
  dplyr::mutate(peak = paste0(PEAKSDIR, sample_name, "_hg38.unstranded.broadPeak")) %>%
  dplyr::filter(file.exists(peak)) %T>% {
    pull(., sample_name) %>%
    sapply(function(x) {
      system(paste0('shuf -n ', NUM_SELECT, " ", paste0(PEAKSDIR, x, "_hg38.unstranded.broadPeak"), 
                    " > ", RLDIR, x, "_hg38.unstranded.broadPeak"))
    })
  }

# Get the HG38 peaks
b_ <- rmap %>%
  dplyr::filter(
    genome == "hg38",
    Condition %in% c("S9.6", "D210N", "FLAG", "delta-HC")
  ) %>%
  dplyr::mutate(peak = paste0(PEAKSDIR, sample_name, "_hg38.unstranded.broadPeak")) %>%
  dplyr::filter(file.exists(peak)) %T>% {
    pull(., sample_name) %>%
      sapply(function(x) {
        system(paste0('cp ', paste0(PEAKSDIR, x, "_hg38.unstranded.broadPeak "), 
                      RL38DIR, x, "_hg38.unstranded.broadPeak"))
      })
  }


## Get the RLFS download and fix the width/score issue
RLFS_PEAKS <- "hg38.rlfs.bed"
RLFS_PEAKS_URI <- "s3://rmapdb-data/misc/rlfs.hg38.fixed.bed"
download.file("https://rmapdb-data.s3.us-east-2.amazonaws.com/rlfs-beds/hg38.rlfs.bed", 
              destfile = "misc/hg38.rlfs.bed")
rlfs <- rtracklayer::import("misc/hg38.rlfs.bed")
names(rlfs) <- paste0("hg38_", seq(width(rlfs)))
score(rlfs) <- 0
rtracklayer::export(rlfs, con = RLFS_PEAKS)
system(paste0("aws s3 --region us-west-2 cp ", RLFS_PEAKS, " ", RLFS_PEAKS_URI, " --acl public-read"))

