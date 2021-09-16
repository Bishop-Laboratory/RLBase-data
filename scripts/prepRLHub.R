#' Script for compiling latest RLBase data for RLHub
dir.create("../RLBase-data/misc-data/rlhub/", showWarnings = FALSE)
library(tidyverse)
library(RLSeq)
library(pbapply)
pbo <- pboptions(type="txt") 

AWS_HTTPS_URL <- "https://rlbase-data.s3.amazonaws.com/"

## 1. Curated Genomic Annotations
message("- 1. Annotations")
dir.create("../RLBase-data/misc-data/rlhub/annotations", showWarnings = FALSE)
file.copy("misc-data/annotations/annotations_hg38.rda",
          "misc-data/rlhub/annotations/annotations_hg38.rda", overwrite = TRUE)
file.copy("misc-data/annotations/annotations_mm10.rda",
          "misc-data/rlhub/annotations/annotations_mm10.rda", overwrite = TRUE)

## 2. RL Regions
message("- 2. RL-Regions")
dir.create("misc-data/rlhub/rlregions/", showWarnings = FALSE)
opts <- c("All", "dRNH", "S96")
pblapply(opts, function(opt) {
  # Get the rlregions table
  fl <- paste0("rlbase-data/misc/rlregions_", opt, "_table.tsv")
  outtbl <- paste0("misc-data/rlhub/rlregions/", opt, "_table.rda")
  rlregions_table <- read_tsv(fl, show_col_types=FALSE, progress=FALSE)
  
  expcor <- paste0("rlbase-data/misc/", opt, "_rlExpCorr.tsv")
  if (opt != "dRNH") {
    exp <- read_tsv(expcor, show_col_types = FALSE, progress = FALSE)
    rlregions_table <- left_join(rlregions_table, exp) %>%
      dplyr::select(-chrom) 
    save(rlregions_table, compress = "xz", file = outtbl)
  } else {
    rlregions_table <- rlregions_table %>%
      mutate(
        corrR = NA,
        corrPVal = NA,
        corrPAdj = NA
      ) 
    save(rlregions_table, compress = "xz", file = outtbl)
      
  }
  
  # Get the annotated rlregions table
  outtbl <- paste0("misc-data/rlhub/rlregions/", opt, "_annotations.rda")
  fl <- paste0("rlbase-data/misc/rlregions_", opt, "_annotations.csv")
  rlregions_annotated <- read_csv(fl, show_col_types = FALSE, progress = FALSE)
  save(rlregions_annotated, compress = "xz", file = outtbl)
})

## 3. RL Samples
message("- 3. RL-Samples")
dir.create("misc-data/rlhub/rlsamples", showWarnings = FALSE)
# Get the condition names for matching exp to rl
expToCond <- read_csv("misc-data/rl_to_exp.csv", show_col_types = FALSE) %>%
  right_join(read_tsv("rlbase-data/rlbase_samples.tsv", show_col_types = FALSE), 
             by = c("rl" = "experiment")) %>%
  mutate(toChoose = gsub(toChoose, pattern = "_", replacement =  "|"))
matchConds <- pbsapply(seq(nrow(expToCond)), function(i) {
  if (! is.na(expToCond$toChoose[i])) {
    cols <- paste0(colnames(expToCond)[grep(colnames(expToCond[i,]),
                                            pattern = expToCond$toChoose[i])])
    dplyr::select(expToCond[i,], all_of(cols[order(cols)])) %>% 
      pivot_longer(everything()) %>% 
      summarise(paste0(value, collapse = "_")) %>% 
      pull(1)
  } else {
    NA
  }
})
expToCond$matchCond <- matchConds

# Summarise to RL-sample level and build links
rlsamples <- expToCond %>%
  dplyr::select(-toChoose, rlsample = rl, expsamples=exp) %>%
  group_by(rlsample) %>%
  mutate(expsamples = paste0(expsamples, collapse=",")) %>%
  ungroup() %>%
  dplyr::filter(group == "rl") %>%
  distinct() %>%
  dplyr::rename(exp_matchCond = matchCond) %>%
  relocate(expsamples, .after = numPeaks) %>%
  mutate(
    coverage_s3 = paste0("coverage/", rlsample, "_", genome, ".bw"),
    peaks_s3 = paste0("peaks/", rlsample, "_", genome, ".broadPeak"),
    fastq_stats_s3 = paste0("fastq_stats/", rlsample, "_", genome, "__fastq_stats.json"),
    bam_stats_s3 = paste0("bam_stats/", rlsample, "_", genome, "__bam_stats.txt"),
    report_html_s3 = paste0("report/html/", rlsample, "_", genome, ".html"),
    report_rda_s3 = paste0("report/rda/", rlsample, "_", genome, ".rda"),
    rlfs_rda_s3 = paste0("rlfs_rda/", rlsample, "_", genome, ".rlfs.rda")
  ) 
save(rlsamples, file = "misc-data/rlhub/rlsamples/rlsamples.rda", compress = "xz")

## 4. Get the RLFS res
message("- 4. RLFS-res")
dir.create("misc-data/rlhub/rlfsres", showWarnings = FALSE)

# Magic
MANIFEST <- "rlbase-data/rlpipes-out/config.tsv"
MANIFEST_FINAL <- "rlbase-data/rlbase_samples.tsv"
RLFSRDA <- "rlbase-data/misc/rlfsRes.rda"
TODISCARD <- "misc-data/todiscard.rda"

# Load the rlfsRes
load(RLFSRDA)

# Load the models
load("misc-data/model/fftModel.rda")
load("misc-data/model/prepFeatures.rda")

# Get the predictions
rlfsPred <- pbapply::pblapply(rlfsRes, predictCondition, 
                              prepFeatures=prepFeatures, fftModel=fftModel)

# Clean large function
rlfsData <- lapply(rlfsRes, FUN = function(x) {
  x$perTestResults$`regioneR::numOverlaps`$evaluate.function <- NULL
  x
})

# Combine
rlfsres <- lapply(names(rlfsData), function(x) {
  list(
    "rlfsData" = rlfsData[[x]],
    "rlfsPred" = rlfsPred[[x]]
  )
})
names(rlfsres) <- names(rlfsData)

# Save
save(rlfsres, file = "misc-data/rlhub/rlfsres/rlfsres.rda", compress = "xz")

## 5. Models
message("- 5. Models")
dir.create("misc-data/rlhub/models", showWarnings = FALSE)
file.copy("misc-data/model/fftModel.rda", overwrite = TRUE,
          to = "misc-data/rlhub/models/fftModel.rda")
file.copy("misc-data/model/prepFeatures.rda", overwrite = TRUE,
          to = "misc-data/rlhub/models/prepFeatures.rda")

## 6. GSG Correlation
message("- 6. GS Corr")
dir.create("misc-data/rlhub/gsg_correlation/", showWarnings = FALSE)

load("misc-data/gsSignalRLBase.rda")
gsSignalRLBase <- gsSignalRMapDB
save(gsSignalRLBase, file = "misc-data/rlhub/gsg_correlation/gsSignalRLBase.rda", compress = "xz")

## 7. Feature Enrichment Test Results
message("- 7. Feature enrichment")
dir.create("misc-data/rlhub/feature_enrichment/", showWarnings = FALSE)

# Feature enrichment of individual peaks
load("rlbase-data/misc/annotatedPeaks.rda")
feature_enrichment_per_sample <- resAnno
save(feature_enrichment_per_sample, file = "misc-data/rlhub/feature_enrichment/feature_enrichment_per_sample.rda", compress = "xz")

# Feature enrichment of RL-Regions
load("rlbase-data/misc/annotatedPeaks.rlregions.rda")
feature_enrichment_rlregions <- res
save(feature_enrichment_rlregions, file = "misc-data/rlhub/feature_enrichment/feature_enrichment_rlregions.rda", compress = "xz")

## 8. Gene Expression
message("- 8. Gene expression")
# This done within the buildExpression.R script...
file.copy("misc-data/expression/geneexp.rda",
          "misc-data/rlhub/expression/geneexp.rda", overwrite = TRUE)

## 9. RLRegion Counts (CTS, Norm, VST)
# This done within the rlregionCountMat.R script...


