#' Script for compiling latest RLBase data for RLHub
library(tidyverse)
library(Rsubread)
library(pbapply)
pbo <- pboptions(type="txt") 

### Preliminaries ###

dir.create("misc-data/rlhub", showWarnings = FALSE)

# Directory where RMapDB bigWigs were downloaded to
# This is the result of running in the shell:
# `aws s3 sync s3://rmapdb-data/coverage/ rlbase-data/rlpipes-out/coverage/`
BW_FOLDER <- "rlbase-data/rlpipes-out/coverage/"

# Number of cores to use for parallel operations
args <- commandArgs(trailingOnly = TRUE)
if (! interactive()) {
  CORES <- args[1]
  rerun <- args[2] == 1
} else {
  CORES <- 44
  rerun <- FALSE
}

#####################

# Load the features
union_regions <- read_tsv("rlbase-data/misc/rlregions_All_table.tsv")
pattern <- "(.+):(.+)\\-(.+):(.+)"
cts_regions <- union_regions %>%
  dplyr::select(GeneID=rlregion, location) %>%
  mutate(
    Chr=gsub(location, pattern = pattern, replacement = "\\1"),
    Start=as.numeric(gsub(location, pattern = pattern, replacement = "\\2")),
    End=as.numeric(gsub(location, pattern = pattern, replacement = "\\3")),
    Strand=gsub(location, pattern = pattern, replacement = "\\4")
  ) %>%
  dplyr::select(-location)
  

# Get seq info
samps <- read_tsv("../RLBase-data/rlbase-data/rlbase_samples.tsv", show_col_types = FALSE) %>%
  dplyr::filter(genome == "hg38", group == "rl") %>%
  dplyr::select(experiment, label, strand_specific, paired_end, mode, prediction, discarded, numPeaks) %>%
  unique() %>%
  mutate(bam = paste0("../RLBase-data/rlbase-data/rlpipes-out/bam/", experiment, "/", experiment, "_hg38.bam"),
         bam_avail = file.exists(bam)) %>%
  dplyr::filter(bam_avail)


message("Starting Feature Counts")
samps$paired_end[samps$experiment == "SRX9684573"] <- FALSE  # Correct a mistaken one
# https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
chunks <- chunk2(x = seq(nrow(samps)), n = 3)
fcResLst <- lapply(names(chunks), function(x) {
  if (! file.exists( paste0("../RLBase-data/tmp/", x, "_fcRes.rda")) | rerun) {
    chunk <- chunks[[x]]
    fcRes <- featureCounts(
      samps$bam[chunk], 
      annot.ext =  cts_regions,
      allowMultiOverlap = TRUE,  # Reads may span multiple rlregions easily
      minMQS = 10,
      isPairedEnd = samps$paired_end[chunk],
      nthreads=CORES
    )
    save(fcRes, file = paste0("../RLBase-data/tmp/", x, "_fcRes.rda"))
  } else {
    load(paste0("../RLBase-data/tmp/", x, "_fcRes.rda"))
  }
  return(fcRes)
})

# Combine the results from rsubread
message("Collating results..")
cts <- lapply(
  fcResLst, 
  function(x) {
    x$counts
  }
) %>% do.call("cbind", .) 
colnames(cts) <- gsub(colnames(cts), pattern = "_hg38.bam", replacement = "")

# Get log2 TPM
# From https://support.bioconductor.org/p/91218/
rllengths <- cts_regions$End - cts_regions$Start
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

rlregions_counts <- SummarizedExperiment::SummarizedExperiment(
  assays = ctsLst, colData = samps, rowData = SummarizedExperiment::rowData(dds)
) 

save(rlregions_counts, 
     file = "../RLBase-data/misc-data/rlhub/rlregions/rlregions_counts.rda", 
     compress = "xz")





