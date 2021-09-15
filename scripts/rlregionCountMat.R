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
if (! interactive()) {
  CORES <- args[1]
} else {
  CORES <- 44
}

#####################

# Combined S9.6 and dRNH consensus sites
opts <- c("dRNH", "S96")
rlrs <- lapply(opts, function(opt) {
  message(opt)
  gpat <- "(.+):(.+)\\-(.+):(.+)"
  
  # Wrangle granges from rlregions table
  gr <- paste0("../RLBase-data/rlbase-data/misc/rlregions_", opt, "_table.tsv") %>%
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
samps <- read_tsv("../RLBase-data/rlbase-data/rlbase_samples.tsv") %>%
  filter(genome == "hg38", group == "rl") %>%
  select(experiment, condType, strand_specific, paired_end, mode, verdict, discarded, numPeaks) %>%
  unique() %>%
  mutate(bam = paste0("../RLBase-data/rlbase-data/rlpipes-out/bam/", experiment, "/", experiment, "_hg38.bam"),
         bam_avail = file.exists(bam)) %>%
  filter(bam_avail)


message("Starting Feature Counts")
samps$paired_end[samps$experiment == "SRX9684573"] <- FALSE  # Correct a mistaken one
# https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
chunks <- chunk2(x = seq(nrow(samps)), n = 3)
fcResLst <- lapply(names(chunks), function(x) {
  chunk <- chunks[[x]]
  if (! file.exists(paste0("../RLBase-data/tmp/", x, "_fcRes.rda"))) {
    fcRes <- featureCounts(
      samps$bam[chunk], 
      annot.ext =  rename(rlregions_union, Chr = seqnames, GeneID = rlrID,
                          Start = start, End = end, Strand = strand),
      allowMultiOverlap = TRUE,  # Reads may span multiple rlregions easily
      minMQS = 10,
      isPairedEnd = samps$paired_end[chunk],
      nthreads=CORES
    )
    save(fcRes, file = paste0("../RLBase-data/tmp/", x, "_fcRes.rda"))
  } else {
    message("Already done")
    load(paste0("../RLBase-data/tmp/", x, "_fcRes.rda"))
    fcRes
  }
})

# Combine the results from rsubread
cts <- lapply(
  fcResLst, 
  function(x) {
    x$counts
  }
) %>% do.call("cbind", .) 
rownames(cts)
colnames(cts) <- gsub(colnames(cts), pattern = "_hg38.bam", replacement = "")

# Get log2 TPM
# From https://support.bioconductor.org/p/91218/
rllengths <- rlregions_union$width
x <- cts/rllengths
tpm <- log2(t( t(x) * 1e6 / colSums(x) ) + 1)

# Wrangle into DESeq2 dataset and get vst
dds <- DESeq2::DESeqDataSetFromMatrix(
  cts + 1, samps, design = ~1
)
vsd <- DESeq2::vst(dds)

ctsLst <- list(
  "cts" = cts,
  "vst" = vsd@assays@data@listData[[1]],
  "tpm" = tpm
) 
rlregions_counts <- lapply(seq(ctsLst), function(i) {
  ctsNow <- ctsLst[[i]]
  typeCts <- names(ctsLst)[i]
  ctsNow %>%
    as.data.frame() %>%
    rownames_to_column(var = "rlregion") %>%
    pivot_longer(
      values_to = typeCts, names_to = "experiment",
      cols = -1
    )
}) %>% purrr::reduce(inner_join, by = c("rlregion", "experiment"))

rlregions_counts <- SummarizedExperiment(
  assays = ctsLst, colData = samps, rowData = rowData(dds)
) 

save(rlregions_counts, 
     file = "../RLBase-data/misc-data/rlhub/rlregions/rlregions_counts.rda", 
     compress = "xz")





