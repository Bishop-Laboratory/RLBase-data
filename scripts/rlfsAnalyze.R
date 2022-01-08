#' Script which takes a set of peaks and performs RLFS analysis on them
#' Usage: `Rscript rlfsAnalyze.R indir outdir cores`

# Get args
args <- commandArgs(trailingOnly = TRUE)

# Get output dir
if (interactive()) {
  indir <- "rlbase-data/rlpipes-out/peaks/"
  outdir <- "rlbase-data/rlpipes-out/rlfs_rda/"
  cores <- 44
} else {
  indir <- args[1]
  outdir <- args[2]
  cores <- args[3]
}
miscpath <- file.path(dirname(dirname(outdir)), "misc")

message(indir)
message(outdir)
message(cores)

# Make output dir
dir.create(outdir, showWarnings = FALSE)
dir.create(miscpath, showWarnings = FALSE)

# Calculate RLFS enrichment for each peak for which...
# 1. A masked genome is available
# 2. The file size > 0 & for which a .broadPeak file is available
available <- c("hg38", "mm10")
peaks <- list.files(indir, full.names = TRUE, pattern = paste0(available, collapse = "|"))
peaks <- peaks[grepl(peaks, pattern = ".+\\.broadPeak$") & file.size(peaks) > 0]

# Calculate enrichment
rlfsRes <- parallel::mclapply(
  seq(peaks), function(i) {
  
  message(i, " / ", length(peaks))
  
  peak <- peaks[i]
  genome <- gsub(peak, pattern = ".*([ES]{1}RX[0-9]+)_([a-zA-Z0-9]+)\\.broadPeak", replacement = "\\2")
  id <- gsub(peak, pattern = ".*([ES]{1}RX[0-9]+)_([a-zA-Z0-9]+)\\.broadPeak", replacement = "\\1")
  out <- file.path(outdir, paste0(id, "_", genome, ".rlfs.rda"))
  
  if (! file.exists(out)) {
    rlr <- RLSeq::RLRanges(peaks = peak, genome = genome)
    res <- RLSeq::analyzeRLFS(rlr)
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
} , mc.cores = cores)

# Get names
names(rlfsRes) <- sapply(rlfsRes, function(x) {
  purrr::pluck(x, "id")
})

# Compile for saving
save(rlfsRes, file = file.path(miscpath, "rlfsRes.rda"), compress = "xz")

