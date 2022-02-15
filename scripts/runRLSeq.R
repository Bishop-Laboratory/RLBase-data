#' This script re-runs RLSeq on all rlsamples
library(pbapply)
library(tidyverse)
library(RLSeq)
pbo <- pboptions(type="txt") 

if (! interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  CORES <- args[1]
} else {
  CORES <- 40
}


dir.create("misc-data/rlranges", showWarnings = FALSE)
dir.create("misc-data/reports", showWarnings = FALSE)

# Summarise to RL-sample level and build links
rlsamples <- RLHub::rlbase_samples()

# Run RLSeq
res2 <- parallel::mclapply(
  seq(rlsamples$rlsample), function(i) {
    
    message(i)
    
    if (! file.exists(paste0(
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
    if (! file.exists(paste0(
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
