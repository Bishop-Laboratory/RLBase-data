library(tidyverse)

# Input
if (! interactive()) {
  nps <- snakemake@input[['peaks']]
  manifests <- snakemake@input[['manifests']]
  output <- snakemake@output
  threads <- snakemake@threads
  patpeak <- "rlregions/rlregions_(.+)\\.narrowPeak"
  patmani <- "rlregions/hg38_manifest_(.+)\\.csv"
} else {
  nps <- list(
    "rlbase-data/rlregions/rlregions_dRNH.narrowPeak",
    "rlbase-data/rlregions/rlregions_S96.narrowPeak"
  )
  manifests <- list(
    "rlbase-data/rlregions/hg38_manifest_dRNH.csv",
    "rlbase-data/rlregions/hg38_manifest_S96.csv"
  )
  threads <- 44
  output <- "rlbase-data/rlregions/rlregions_All.narrowPeak"
  patpeak <- "rlbase-data/rlregions/rlregions_(.+)\\.narrowPeak"
  patmani <- "rlbase-data/rlregions/rlregion_manifest_(.+)\\.csv"
}

## Combined S9.6 and dRNH consensus sites ##

# Extract type (RNH or S96)
names(nps) <- sapply(nps, function(x) {
  gsub(x, pattern = patpeak, replacement = "\\1")
})

# Get number of samples from manifests
names(manifests) <- sapply(manifests, function(x) {
  gsub(x, pattern = patmani, replacement = "\\1")
})
numsamps <- sapply(manifests, FUN = function(x) {
  read_csv(x, col_names = FALSE, show_col_types = FALSE) %>% nrow()
})

# Load peaks
rlrs <- lapply(names(nps), function(opt) {
  # Wrangle granges from rlregions table
  nps[[opt]] %>%
    read_tsv(show_col_types = FALSE, progress = FALSE, 
             skip = 1, col_names = c(
               "chrom", "start", "end", "name", "score", "strand", 
               "signalVal", "pVal", "qVal", "peak"
             )) %>% 
    mutate(type = {{ opt }},
           score = score / (10 * numsamps[[opt]])) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}) 
names(rlrs) <- names(nps)

# Combine into a singal GRanges
combo <- c(rlrs$dRNH, rlrs$S96)

# Store type column
types <- combo$type
scores <- combo$score

# Reduce with revmap
rvmp <- combo %>%
  GenomicRanges::reduce(with.revmap=TRUE) 

# Match revmap to original types
revMatch <- parallel::mclapply(seq(rvmp$revmap), function(i) {
  x <- rvmp$revmap[[i]]
  list( 
    "type" = paste0(unique(types[x]), collapse = " "), 
    "score" = mean(scores[x]))
}, mc.cores = threads)
labs <- sapply(revMatch, function(x) {x$type})
scores <-  round(sapply(revMatch, function(x) {x$score}) * 100, digits = 2)


# Wrangle
rlregions_union <- rvmp %>%
  unique() %>%
  tibble::as_tibble() %>%
  mutate(name = paste0("RLRegion_", seq(seqnames)),
         sources = unlist({{ labs }}),
         score = {{ scores }},
         strand = ".",
         sources = case_when(
           sources == "S96 dRNH" ~ "dRNH S96",
           TRUE ~ sources
         )) %>%
  select(-revmap) %>%
  relocate(name, .before = width) %>%
  relocate(score, .before = strand) %>%
  select(-width)

readr::write_tsv(rlregions_union, file = output[[1]], col_names = FALSE)

message("Done")
