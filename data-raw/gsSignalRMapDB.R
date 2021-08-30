# Script for generating gsSignalRMapDB
library(magrittr)

### Preliminaries ###

# Minimum allowed size of ranges after liftover
# Number is based on expected ranges size of 100kb
MIN_LIFTED_SIZE <- 10000  

# Directory where RMapDB bigWigs were downloaded to
# This is the result of running in the shell:
# `aws s3 sync s3://rmapdb-data/bigwigs/rseq-coverage-unstranded/ RMapDB_BiwWigs/`
BW_FOLDER <- "data-raw/RMapDB_BigWigs/"

# S3 URI and HTTPS URL where the bigWigs are
BWHTTPS_BASE <- "https://rmapdb-data.s3.us-east-2.amazonaws.com/bigwigs/rseq-coverage-unstranded/"
BWS3_BASE <- "s3://rmapdb-data/bigwigs/rseq-coverage-unstranded/"

# correlation_genes_100kb.bed was graciously provided by Stella Hartono and Fred Chedin
# These intervals correspond to the +/-50kb coordinated of sites profiled by SMRF-Seq
# in 10.1016/j.jmb.2020.02.014 and used in their recent EMBO paper for correlation 
# analysis: https://doi.org/10.15252/embj.2020106394
GSBED_HG19 <- "data-raw/correlation_genes_100kb.bed"

# Number of cores to use for parallel operations
CORES_TO_USE <- 6

#####################

### Perform the analysis ###

# Get the GS regions as a granges
hg19gs <- rtracklayer::import(GSBED_HG19)

# They were lifted from hg19 to hg38
chain <- RSeqR:::getChain(genomeFrom = "hg19", genomeTo = "hg38")

# Lift over
ranges <- GenomicRanges::reduce(hg19gs)
names(ranges) <- seq(GenomicRanges::start(ranges))
lifted <- unlist(rtracklayer::liftOver(ranges, chain = chain)) 

# Clean up the ranges 
liftedRed <- GenomicRanges::reduce(lifted)
hg38gs <- liftedRed[GenomicRanges::width(liftedRed) > MIN_LIFTED_SIZE,]

# Windows of 1kb
hg38Tbl <- hg38gs %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  dplyr::rename(chrom = .data$seqnames)
hg38wds <- valr::bed_makewindows(hg38Tbl, win_size = 1000) %>%
  dplyr::mutate(chrom = as.character(.data$chrom))

# Make into GRanges
positions <- hg38wds %>%
  GenomicRanges::makeGRangesFromDataFrame() %>%
  rtracklayer::BigWigSelection()

# List the available bw files and wrangle
# Alternatively just download them
# `aws s3 sync s3://rmapdb-data/bigwigs/rseq-coverage-unstranded/ RMapDB_BiwWigs/`
bwFiles <- system(paste0("aws s3 ls ", BWS3_BASE), intern = TRUE)
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
resLst <- lapply(seq(bws$id), function(rowNow) {
  
  # Get current values
  idNow <- bws$id[rowNow]
  bwFile <- bws$bw[rowNow]
  
  # Print
  message(idNow, " - ", rowNow, "/", length(bws$id))
  
  # Read in the bigWig file using these locations
  bwCon <- rtracklayer::BigWigFile(
    paste0(BWHTTPS_BASE, bwFile)
  )
  bw <- rtracklayer::import.bw(con = bwCon, 
                               selection = positions)
  
  # Wrangle BW into tibble
  bw <- bw %>%
    as.data.frame() %>%
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
  
}) 

# Wrangle
gsSignalRMapDB <- dplyr::bind_rows(resLst) %>%
  tidyr::pivot_wider(id_cols = .data$location,
                     names_from = .data$id,
                     values_from = .data$value)
  
############################

# Save it
usethis::use_data(gsSignalRMapDB, overwrite = TRUE, compress = "xz")
