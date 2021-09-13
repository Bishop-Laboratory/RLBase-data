#' Script for compiling latest RLBase data for RLHub
library(tidyverse)
library(Rsubread)
library(pbapply)
pbo <- pboptions(type="txt") 

### Preliminaries ###

# Directory where RMapDB bigWigs were downloaded to
# This is the result of running in the shell:
# `aws s3 sync s3://rmapdb-data/coverage/ rlbase-data/rlpipes-out/coverage/`
BW_FOLDER <- "rlbase-data/rlpipes-out/coverage/"

# Number of cores to use for parallel operations
args <- commandArgs(trailingOnly = TRUE)
CORES <- args[1]

#####################

# Combined S9.6 and dRNH consensus sites
opts <- c("dRNH", "S96")
rlrs <- lapply(opts, function(opt) {
  message(opt)
  gpat <- "(.+):(.+)\\-(.+):(.+)"
  
  # Wrangle granges from rlregions table
  gr <- paste0("rlregions_", opt, "_table.tsv.xz") %>%
    read_tsv(show_col_types = FALSE, progress = FALSE) %>%
    filter(! is_repeat) %>%
    select(rlregion, location) %>%
    mutate(seqnames = gsub(location, pattern = gpat, replacement = "\\1"),
           start = as.numeric(gsub(location, pattern = gpat, replacement = "\\2")),
           end =  as.numeric(gsub(location, pattern = gpat, replacement = "\\3")),
           strand =  gsub(location, pattern = gpat, replacement = "\\4")) %>%
    select(-location) %>%
    as.data.frame() %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}) 
names(rlrs) <- opts
combo <- c(rlrs$dRNH, rlrs$S96)
types <- gsub(combo$rlregion, pattern = "(.+)_.+", replacement = "\\1")
rvmp <- combo %>%
  GenomicRanges::reduce(with.revmap=TRUE) 
labs <- parallel::mclapply(seq(rvmp$revmap), function(i) {
  x <- rvmp$revmap[[i]]
  paste0(unique(types[x]), collapse = " ")
}, mc.cores = CORES)
rlregions_union <- rvmp %>%
  unique() %>%
  tibble::as_tibble() %>%
  mutate(rlrID = paste0("OL_RLRegion_", seq(seqnames)),
         sources = unlist(!! labs),
         sources = case_when(
           sources == "S96 dRNH" ~ "dRNH S96",
           TRUE ~ sources
         )) %>%
  select(-revmap) %>%
  relocate(rlrID, .before = seqnames) 
save(rlregions_union, file = "rlbase-data/misc/rlregions_union.rda", compress = "xz")
save(rlregions_union, file = "misc-data/rlhub/rlregions/rlregions_union.rda", compress = "xz")


# Get seq info
samps <- read_tsv("rlbase-data/rlbase_samples.tsv") %>%
  filter(genome == "hg38", group == "rl") %>%
  select(experiment, condType, strand_specific, paired_end) %>%
  unique() %>%
  mutate(bam = paste0("rlbase-data/rlpipes-out/bam/", experiment, "/", experiment, "_hg38.bam"),
         bam_avail = file.exists(bam)) %>%
  filter(bam_avail)

# Use RSubRead to get the counts within each
fcRes <- featureCounts(
  samps$bam, 
  annot.ext =  rename(rlregions_union, Chr = seqnames, GeneID = rlrID,
                      Start = start, End = end, Strand = strand),
  allowMultiOverlap = TRUE,  # Reads may span multiple rlregions easily
  minMQS = 10,
  isPairedEnd = samps$paired_end
)
save(fcRes, file = "tmp/fcRes.rda")

