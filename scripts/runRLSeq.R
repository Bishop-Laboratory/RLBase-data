#' This script re-runs RLSeq on all rlsamples
library(pbapply)
library(tidyverse)
library(RLSeq)
pbo <- pboptions(type="txt") 

if (! interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  CORES <- args[1]
} else {
  CORES <- 10
}


dir.create("misc-data/rlranges", showWarnings = FALSE)
dir.create("misc-data/reports", showWarnings = FALSE)

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
    report_html_s3 = paste0("reports/", rlsample, "_", genome, ".html"),
    rlranges_rds_s3 = paste0("rlranges/", rlsample, "_", genome, ".rds"),
    rlfs_rda_s3 = paste0("rlfs_rda/", rlsample, "_", genome, ".rlfs.rda")
  ) 


# Run RLSeq
res2 <- parallel::mclapply(
  seq(rlsamples$rlsample), function(i) {
    
    message(i)
    
    if (TRUE | ! file.exists(paste0(
      "misc-data/rlranges/", rlsamples$rlsample[i],
      "_", rlsamples$genome[i], ".rds"
    ))) {
      
      rlr <- try(
        RLSeq::RLRanges(
          peaks = paste0("rlbase-data/rlpipes-out/peaks/", rlsamples$rlsample[i],
                         "_", rlsamples$genome[i], ".broadPeak"),
          coverage =  paste0("rlbase-data/rlpipes-out/coverage/", rlsamples$rlsample[i],
                             "_", rlsamples$genome[i], ".bw"),
          genome = rlsamples$genome[i],
          mode = rlsamples$mode[i],
          label = rlsamples$label[i],
          sampleName = rlsamples$rlsample[i]
        ), silent = TRUE
      )
      if ("try-error" %in% class(rlr)) {
        return(NULL)
      } 
      
      # Run RLSeq
      rlr2 <- try(RLSeq::RLSeq(rlr, quiet = TRUE), silent = TRUE)
      if ("try-error" %in% class(rlr2)) {
        if (grepl(rlr2[1], pattern = "not available in mask list")) {
          rlr2 <- try(RLSeq::RLSeq(rlr, quiet = TRUE, useMask=FALSE), silent = TRUE)
        } 
        if ("try-error" %in% class(rlr2)) {
          return(NULL)
        } else {
          rlr <- rlr2
        }
      } else {
        rlr <- rlr2
      }
      # Add new prediction...
      PREPMODEL <- "misc-data/model/prepFeatures.rda"
      FFTMODEL <- "misc-data/model/fftModel.rda"
      load(PREPMODEL)
      load(FFTMODEL)
      rlr <- RLSeq:::predictCondition(
        rlr, prepFeatures=prepFeatures,
        fftModel=fftModel
      )
      
      # Save RDS
      saveRDS(
        rlr,
        file = paste0(
          "misc-data/rlranges/", rlr@metadata$sampleName,
          "_", GenomeInfoDb::genome(rlr)[1], ".rds"
        ),
        compress = "xz"
      )
    } else {
      rlr <- readRDS(
        paste0(
          "misc-data/rlranges/", rlsamples$rlsample[i],
          "_", rlsamples$genome[i], ".rds"
        )
      )
    }
    
    
    # Report
    if (FALSE | ! file.exists(paste0(
      "misc-data/reports/", rlr@metadata$sampleName,
      "_", GenomeInfoDb::genome(rlr)[1], ".html"
    )))  {
      try(
        RLSeq::report(
          rlr,
          reportPath = paste0(
            "misc-data/reports/", rlr@metadata$sampleName,
            "_", GenomeInfoDb::genome(rlr)[1], ".html"
          ), quiet = TRUE
        ), silent = TRUE)
    } else {
      return(TRUE)
    }
  }, mc.cores = CORES
)
