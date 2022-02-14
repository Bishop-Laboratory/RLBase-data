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
cores <- 20

# Load the rlfsRes
load(RLFSRDA)

# Load the difficult-to-classify ones
load(TODISCARD)

# Load all the models
load(PREPMODEL)
load(FFTMODEL)

# Get additional rlfsRes
indir <- "rlbase-data/rlpipes-out/peaks/"
outdir <- "rlbase-data/rlpipes-out/rlfs_rda/"
miscpath <- file.path(dirname(dirname(outdir)), "misc")

peaks <- list.files(indir, full.names = TRUE)
peaks <- peaks[grepl(peaks, pattern = ".+\\.broadPeak$") & file.size(peaks) > 0] %>% unique()
rlfsRes2 <- parallel::mclapply(
  seq(peaks), function(i) {
    message(i, " / ", length(peaks))
    peak <- peaks[i]
    genome <- gsub(peak, pattern = ".*([ES]{1}RX[0-9]+)_([a-zA-Z0-9]+)\\.broadPeak", replacement = "\\2")
    id <- gsub(peak, pattern = ".*([ES]{1}RX[0-9]+)_([a-zA-Z0-9]+)\\.broadPeak", replacement = "\\1")
    out <- file.path(outdir, paste0(id, "_", genome, ".rlfs.rda"))
    
    if (! file.exists(out)) {
      rlr <- RLSeq::RLRanges(peaks = peak, genome = genome)
      res <- RLSeq::analyzeRLFS(rlr, useMask = FALSE)
      res <- RLSeq::rlresult(res, "rlfsRes")
      # Remove large randomization function...
      res$perTestResults$`regioneR::numOverlaps`$randomize.function <- NULL
      res$id <- id
      res$genome <- genome
      
      save(res, file = out, compress = "xz")
      return(res)
    } else {
      message("Already run!")
      load(out)
      return(res)
    }
  }, 
  mc.cores = cores
)

# Get names
names(rlfsRes2) <- sapply(rlfsRes2, function(x) {
  purrr::pluck(x, "id")
})

rlfsRes_full <- rlfsRes2

save(rlfsRes_full, file = file.path(miscpath, "rlfsRes__full.rda"), compress = "xz")

# Get the predictions
pred <- pbapply::pblapply(rlfsRes_full, function(x) {
  predictCondition(
    rlfsRes = x,
    prepFeatures=prepFeatures,
    fftModel=fftModel
  )
})

# Extract label predicted
verd <- sapply(pred, purrr::pluck, "prediction")
verdtb <- tibble(
  experiment = names(verd),
  prediction = verd
)

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
