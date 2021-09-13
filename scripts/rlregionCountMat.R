#' Script for compiling latest RLBase data for RLHub
library(tidyverse)
library(RLSeq)
library(pbapply)
pbo <- pboptions(type="txt") 

### Preliminaries ###

# Directory where RMapDB bigWigs were downloaded to
# This is the result of running in the shell:
# `aws s3 sync s3://rmapdb-data/coverage/ rlbase-data/rlpipes-out/coverage/`
BW_FOLDER <- "rlbase-data/rlpipes-out/coverage/"

# Number of cores to use for parallel operations
CORES_TO_USE <- 6

#####################

# Combined S9.6 and dRNH consensus sites
opts <- c("dRNH", "S96")
rlrs <- lapply(opts, function(opt) {
  message(opt)
  gpat <- "(.+):(.+)\\-(.+):(.+)"
  
  # Wrangle granges from rlregions table
  gr <- paste0("rlbase-data/misc/rlregions_", opt, "_table.tsv") %>%
    read_tsv(show_col_types = FALSE, progress = FALSE) %>%
    select(rlregion, location) %>%
    mutate(seqnames = gsub(location, pattern = gpat, replacement = "\\1"),
           start = as.numeric(gsub(location, pattern = gpat, replacement = "\\2")),
           end =  as.numeric(gsub(location, pattern = gpat, replacement = "\\3")),
           strand =  gsub(location, pattern = gpat, replacement = "\\4")) %>%
    select(-location) %>%
    as.data.frame() %>%
    GenomicRanges::makeGRangesFromDataFrame()
}) 
names(rlrs) <- opts
rlru <- GenomicRanges::union(rlrs$dRNH, rlrs$S96) %>%
  GenomicRanges::reduce() %>%
  unique()

# List the available bw files and wrangle
# Alternatively just download them
bwFiles <- list.files(BW_FOLDER, full.names = TRUE)
bws <- tibble::tibble(
  bw = gsub(bwFiles, pattern = ".+([ES]{1}RX[0-9]+.+)$", replacement = "\\1")
) %>%
  dplyr::mutate(id = gsub(.data$bw, pattern = "(.+)_(.+)\\.bw",
                          replacement = "\\1"),
                genome = gsub(.data$bw, pattern = "(.+)_(.+)\\.bw",
                              replacement = "\\2")) %>%
  dplyr::filter(genome == "hg38")

# Extract the coverage across the gs region windows
# NOTE: this took a long time to finish and segfaulted several times
# This is due to issues in the UCSC bigWig accessors in the rtracklayer 
# package and SSL timeouts. Maybe be easier to just download locally.
resLst <- parallel::mclapply(seq(bws$id)[1], function(rowNow) {
  
  # Get current values
  idNow <- bws$id[rowNow]
  bwFile <- bws$bw[rowNow]
  
  # Read in the bigWig file using these locations
  bwCon <- rtracklayer::BigWigFile(
    paste0(BW_FOLDER, bwFile)
  )
  bw <- rtracklayer::import.bw(con = bwCon, 
                               selection = rlru)
  
  # Wrangle BW into tibble
  bw <- bw %>%
    tibble::as_tibble() %>%
    dplyr::rename(chrom = .data$seqnames) %>%
    dplyr::mutate(chrom = as.character(.data$chrom))
  
  # Summarize across gs intervals with valr
  valr::bed_map(x = hg38wds, y = bw, value = sum(score)) %>%
    dplyr::mutate(location = paste0(
      .data$chrom, "_", .data$start, "_", .data$end
    ),
    id = !! idNow) %>%
    dplyr::select(.data$location, .data$id, .data$value) 
  
}, mc.cores = 44) 

# Wrangle
gsSignalRLBase <- dplyr::bind_rows(resLst) %>%
  tidyr::pivot_wider(id_cols = .data$location,
                     names_from = .data$id,
                     values_from = .data$value)

############################

save(gsSignalRMapDB, file = "misc-data/gsSignalRLBase.rda", compress = "xz")