#' Script for annotating peaks from the RLPipes pipeline

library(RLSeq)
library(tidyverse)

if (! interactive()) {
  # Get args
  args <- commandArgs(trailingOnly = TRUE)
  
  # Get output dir
  PEAK_LOC <- args[1]
  OUTANNO <- args[2]
  CORES <- as.numeric(args[3])
} else {
  CORES <- 44
  PEAK_LOC <- "rlbase-data/rlpipes-out/peaks/"
  OUTANNO <- "rlbase-data/misc/annotatedPeaks.tsv"
}


# Get peaks
peaks <- list.files(PEAK_LOC, pattern = "hg38\\.broadPeak|mm10\\.broadPeak", full.names = TRUE)
peaks <- peaks[file.size(peaks) > 0]
peaks <- tibble(
  experiment = gsub(peaks, pattern = ".+/([ES]{1}RX[0-9]+)_([a-zA-Z0-9]+)\\.broadPeak", replacement = "\\1"),
  genome = gsub(peaks, pattern = ".+/([ES]{1}RX[0-9]+)_([a-zA-Z0-9]+)\\.broadPeak", replacement = "\\2"),
  peakfile = peaks
)

# Run annotations
resAnno <- parallel::mclapply(seq(peaks$experiment), function(i) {
  message(i, " / ", length(peaks$experiment))
  peak <- peaks$peakfile[i] %>%
    regioneR::toGRanges()
  
  anno <- RLSeq::featureEnrich(peak, genome = peaks$genome[i]) %>%
    mutate(experiment = peaks$experiment[i])
}, mc.cores = CORES) %>%
  bind_rows() 
stopifnot(nrow(resAnno) > 10)
save(resAnno,
     file = gsub(OUTANNO, pattern = "\\.tsv", replacement = ".rda"),
     compress = "xz")
write_tsv(resAnno, OUTANNO)
system(paste0("xz -f ", OUTANNO))
