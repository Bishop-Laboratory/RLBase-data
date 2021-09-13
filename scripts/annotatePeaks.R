#' Script for annotating peaks from the RLPipes pipeline

library(RLSeq)
library(tidyverse)
library(pbapply)
pbo <- pboptions(type="txt") 

if (! interactive()) {
  # Get args
  args <- commandArgs(trailingOnly = TRUE)
  
  # Get output dir
  PEAK_LOC <- args[1]
  OUTANNO <- args[2]
  CORES <- as.numeric(args[3])
} else {
  CORES <- 44
  PEAK_LOC <- "../RLBase-data/rlbase-data/rlpipes-out/peaks/"
  OUTANNO <- "../RLBase-data/rlbase-data/misc/annotatedPeaks.tsv"
}

# Get peaks
peaks <- list.files(PEAK_LOC, pattern = "hg38\\.broadPeak|mm10\\.broadPeak", full.names = TRUE)
peaks <- peaks[file.size(peaks) > 0]
peaktbl <- tibble(
  experiment = gsub(peaks, pattern = ".+/([ES]{1}RX[0-9]+)_([a-zA-Z0-9]+)\\.broadPeak", replacement = "\\1"),
  genome = gsub(peaks, pattern = ".+/([ES]{1}RX[0-9]+)_([a-zA-Z0-9]+)\\.broadPeak", replacement = "\\2"),
  peakfile = peaks
)

# Run annotations
annotationLst <- RLSeq::annotations
dir.create("../RLBase-data/tmp/annotatedPeaks")
resAnno <- pblapply(seq(peaktbl$experiment), function(i) {
  peak <- peaktbl$peakfile[i] %>%
    regioneR::toGRanges()
  if (! file.exists(paste0("../RLBase-data/tmp/annotatedPeaks/", peaktbl$experiment[i], ".rda"))) {
    anno <- RLSeq::featureEnrich(peaks = peak, genome = peaktbl$genome[i],
                                 annotations =  annotationLst, cores = CORES) %>%
      mutate(experiment = peaktbl$experiment[i])
    save(anno, file = paste0("../RLBase-data/tmp/annotatedPeaks/", peaktbl$experiment[i], ".rda"))
  } else {
    load(paste0("../RLBase-data/tmp/annotatedPeaks/", peaktbl$experiment[i], ".rda"))
    anno
  }
}) %>% dplyr::bind_rows()
stopifnot(nrow(resAnno) > 10)
save(resAnno,
     file = gsub(OUTANNO, pattern = "\\.tsv", replacement = ".rda"),
     compress = "xz")
write_tsv(resAnno, OUTANNO)
system(paste0("xz -f ", OUTANNO))


