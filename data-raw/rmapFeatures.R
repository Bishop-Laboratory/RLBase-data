# Script for making the rmapFeatures object
# This object contains the quality characteristics for each entry in rmapSamps
# Currently, only human is support for predicting case/control
# Mouse is also supported for typical annotations
library(tidyverse)
library(RSeqR)
library(parallel)
library(regioneR)

# Output locations
RLFS_RDAS <- "data-raw/rlfs_rda/"

# For parallel operations
CORES_USE <- 44

# Make output dir
dir.create(RLFS_RDAS, showWarnings = FALSE)

# S3 Peaks Bucket
peaks_url <- RSeqR:::RMAP_PEAKS_URL
peaks_s3 <- RSeqR:::RMAP_PEAKS_S3

# Get peak locations for every sample
rmapSampsPeaks <- rmapSamps %>%
  filter(genome %in% c("hg38", "mm10")) %>%
  mutate(peakLoc = paste0(peaks_url, id, "_", genome, ".broadPeak"))
  
# Get the available peaks
peaksAvail <- system(paste0("aws s3 ls ", peaks_s3), intern = TRUE)
peaksAvail <- paste0(peaks_url, gsub(peaksAvail[-1], pattern = ".+ ([ES]{1}RX[0-9}+_[a-zA-Z0-9]+\\.broadPeak)", replacement = "\\1"))

# Cross-reference to the rmap samples
rmapSampsPeaks <- filter(rmapSampsPeaks, peakLoc %in% peaksAvail)

# For each sample, get the peaks and then do the RLFS test, prediction, and feature test
featuresPrelim <- mclapply(seq(rmapSampsPeaks$id), function(i) {
  
  # Get current values
  peaks <- rmapSampsPeaks$peakLoc[i]
  id <- rmapSampsPeaks$id[i]
  genome <- rmapSampsPeaks$genome[i]
  output <- paste0(RLFS_RDAS, id, "_", genome, ".rda")
  
  # Print
  message(id, " -- ", i, " / ", length(rmapSampsPeaks$id))
  
  # Read in the data
  data <- readr::read_tsv(peaks, col_names = c("seqnames", "start", "end",
                                               "name", "score", "strand", 
                                               "signalVal", "pVal", "qVal"), 
                          show_col_types = FALSE, progress = FALSE)
  
  # Get GRanges
  gr <- GenomicRanges::makeGRangesFromDataFrame(data)
  
  # Number of peaks
  num_peak <- length(width(gr))
  
  # Run if 
  if (num_peak > 0) {
    if (! file.exists(output)) {
      # Get the RLFS Res
      rlfsRes <- RSeqR::analyzeRLFS(gr, genome = genome, force.parallel=FALSE)
      # Reduce object size
      rlfsRes$perTestResults$`regioneR::numOverlaps`$randomize.function <- NULL
      rlfsRes$perTestResults$`regioneR::numOverlaps`$evaluate.function <- NULL
      # Save
      save(rlfsRes, file = output, compress = "xz")
      
    } else {
      load(output)
    }
    
    # Compile decion-making results
    resdf <- tibble(
      value = c(rlfsRes$perTestResults$`regioneR::numOverlaps`$pval,
                rlfsRes$perTestResults$`regioneR::numOverlaps`$zscore,
                num_peak),
      annotation = c("RLFS__pval", "RLFS__Zmax", "num_peaks"),
      id = id
    )
    
    # Now test for feature enrichment
    featRes <- RSeqR::featureEnrich(gr, genome = genome)
    
    # Pivot and join with resdf -- then return
    featRes %>%
      pivot_longer(cols = -annotation) %>%
      mutate(id = id, annotation = paste0(annotation, "__", name)) %>%
      select(-name) %>%
      bind_rows(resdf) %>%
      pivot_wider(id_cols = id, names_from = annotation, values_from = value)
    
  }
  
}, mc.cores = CORES_USE) 

# Compile list into Tbl and convert to -log10 pval and cleanup NAs
problems <- c("Mt_rRNA", "Mt_tRNA")
featuresPrelimTbl <- bind_rows(featuresPrelim) %>%
  select(! contains(!! problems)) %>%
  mutate(RLFS__pval = -log10(RLFS__pval)) %>%
  # Convert the small number of NaNs to NAs 
  mutate(across(contains("__"), function(x) {ifelse(is.nan(x), NA, x)}))

# Join with rmapSamps 
rmapJoin <- inner_join(rmapSamps, featuresPrelimTbl, by = "id")

# TODO: Include fastp and bam stats reports

## Create manifest for model-building notebook ##
# First, we define eligible cases and controls based on:
# 1. Whether they are labeled as Normal-like (case) or RNH-like (control)
# 2. Cases are refined as those which have significant RLFS pval, > 500 peaks,
#    ZMax > 0, fiveUTR enrichment significant obs/exp > 0.
# The rationale is:
#   a. Cases should be significantly enriched with ZMax > 0 in RLFS as these are known to associate with R-loops found through DRIP-qPCR
#   b. Cases should have > 500 peaks. While this is an arbitrary number, it represents the balance between 
#      s/n ratio (which num peaks measures) and the reality that some modalities simply find fewer peaks (e.g., R-ChIP)
#   c. Cases should have significant over-representation of 5'UTR regions. This is because R-loops form preferenetially
#      during both nascent and elongating transcription. While DRIP can find elongating transcription, R-ChIP is much better
#      at finding nascent transcription. In either case, the 5'UTR will have an R-loop. Therefore, this is a good broad-base
#      sanity check of whether R-loops are really being found. 
samps <- rmapJoin %>%
  mutate(group = case_when(pval < .01 & Condition %in% CASE_CONDS & MACS2__total_peaks > 500 &
                             zscore > 0 &
                             pct_aligned > 60 & `5UTR__LogP enrichment (+values depleted)` < 0 ~ "Case",
                           Condition %in% CTRL_CONDS ~ "Control",
                           TRUE ~ "Unknown")) %T>%
  {
    group_by(., group, mode) %>%
      tally() %>%
      pivot_wider(id_cols = mode,
                  names_from = group,
                  values_from = n,
                  values_fill = 0) %>%
      print()
  }

# Pull out the samples that do not resemble their labels
TO_CENSOR <- c(
  "SRX1025907", "SRX1025909", "SRX1025911", "SRX1025913",
  "SRX113812", "SRX113813", "SRX113814", "SRX2455189",
  "SRX2455191", "SRX2455192", "SRX2455195", "SRX2481503",
  "SRX2642935", "SRX2642936", "SRX2642951", "SRX2642955",
  "SRX2642956", "SRX2642959", "SRX2642960", "SRX2642961",
  "SRX2642963", "SRX2642962", "SRX2642964", "SRX2642970",
  "SRX2733011", "SRX2733012", "SRX2733013", "SRX2733014",
  "SRX2918366", "SRX4457337", "SRX4671914", "SRX4671915",
  "SRX4732954", "SRX4732955", "SRX4732956", "SRX4732957",
  "SRX4776628", "SRX4776630", "SRX534096",
  "SRX550078", "SRX5696400", "SRX5696404","SRX6427715",
  "SRX6427722", "SRX6686232", "SRX6686233", "SRX6686234",
  "SRX6686235", "SRX6686236", "SRX7299431", "SRX7299435",
  "SRX7583977", "SRX7583979", "SRX7583980", paste0("SRX7805", 791:802),
  "SRX8110928", "SRX8110929", "SRX8110930", "SRX8110931"
)


write_csv(samps, "../RSeq-supplemental-analysis/misc/datasets_for_fft_testing/manifest_full_data.csv")

samps %>%
  filter(! id %in% TO_CENSOR) %>%
  filter(group != "Unknown") %T>%
  select(id, group, filename) %>%
  write_csv("../RSeq-supplemental-analysis/misc/datasets_for_fft_testing/manifest.csv")

