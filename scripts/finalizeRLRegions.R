#' This script will wrangle the results of the RLRegions pipeline to create final tables

library(tidyverse)
library(ChIPpeakAnno)
library(pbapply)
pbo <- pboptions(type="txt") 

#' Helper function for converting CSV to GR
csvToGR <- function(csv) {
  readr::read_csv(csv) %>%
    dplyr::mutate(
      seqnames = gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\1"),
      start = gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\2"),
      end = gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\3"),
      strand = gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\4"),
      # from https://stackoverflow.com/questions/45146688/execute-dplyr-operation-only-if-column-exists
      score = if ("confidence_level" %in% colnames(.)) confidence_level else 0
    ) %>%
    dplyr::filter(! grepl(pattern = "LRG.+", x = id)) %>%
    dplyr::select(
      seqnames,
      start,
      end,
      name=id,
      score,
      strand
    ) %>%
    dplyr::distinct(seqnames, start, end, .keep_all = TRUE) %>%
    as.data.frame() %>%
    # Resolves issue of having NAs in the start/end column
    dplyr::mutate(
      start = as.numeric(as.character(start)),
      end = as.numeric(as.character(end))
    ) %>%
    na.omit() %>%
    ChIPpeakAnno::toGRanges() 
}

#' Makes the "Rloops" table
makeRLoops <- function(ORIG_RL, MAX_RL_SIZE, RLOOPS) {
  
  cond <- gsub(RLOOPS, pattern = ".+rlregions_([a-zA-Z0-9]+)\\.csv", replacement = "\\1")
  
  readRL <- function(x, max_rl_size) {
    readr::read_tsv(x, col_names = c(
      "seqnames", "start", "end", "name", "score", "strand", "signalValue", "pval", "qval", "peak"
    ),skip = 1) %>%
      dplyr::filter(end - start < max_rl_size) %>%
      dplyr::mutate(group = {{ cond }}) %>%
      as.data.frame() %>%
      regioneR::toGRanges() %>%
      return()
  }
  
  # Get rloops
  rloops <- ORIG_RL %>%
    readRL(MAX_RL_SIZE)
  
  # Check RLFS, chrom_sizes, and mask
  RLFS <- RLSeq:::getRLFSAnno("hg38")
  
  # Prevent stranded assignment
  GenomicRanges::strand(RLFS) <- "*"
  RLFS <- GenomicRanges::reduce(RLFS)
  
  # Olap
  rloopsRLFS <- list("RLoops" = rloops, "RLFS"=RLFS)
  rlrlfsol <- findOverlapsOfPeaks(rloopsRLFS)
  olNames <- unlist(rlrlfsol$peaklist$`RLoops///RLFS`$peakNames) %>%
    gsub(pattern = ".+__(RLoops\\..+)", replacement = "\\1") 
  rloops <- rlrlfsol$all.peaks$RLoops
  mcols(rloops)$names <- names(rloops)
  
  rloops %>%
    as_tibble() %>%
    rownames_to_column(var = "id") %>%
    mutate(id = paste0(group, "_RL", id),
           type = "Unknown",
           location = paste0(seqnames, ":", start, "-", end, ":", strand),
           confidence_level = score,
           is_rlfs = case_when(
             names %in% !! olNames ~ TRUE,
             TRUE ~ FALSE
           ),
           origName = name) %>%
    dplyr::select(id, type, location, confidence_level, is_rlfs, origName) %>%
    write_csv(RLOOPS)
  return(NULL)
}

#' Makes the "RLoopSignal" table for RMapDB
makeRLoopSignal <- function(RLOOPS, RLOOP_SIGNAL_MAT, RLOOP_OLAP) {
  
  # Get RLoop to peak mapping
  rloops <- read_tsv(RLOOPS) %>%
    dplyr::rename(rloop_id = id)
  rlol <- read_tsv(RLOOP_OLAP, col_names = c("origName", "olPeak", "signalVal", "qVal", "blank")) %>%
    mutate(rmap_sample_id = gsub(olPeak, pattern = "^([ES]{1}RX[0-9]+)_.+", replacement = "\\1")) %>%
    select(-olPeak, -blank) %>%
    full_join(rloops, by = "origName")
  rlol <- rlol %>%
    select(signalVal, qVal, rloop_id, rmap_sample_id) %>%
    group_by(rloop_id, rmap_sample_id) %>%
    summarise(
      numOlap = n(),
      qVal = max(qVal),
      signalVal = max(signalVal)
    )
  
  
  # Get RLoops
  rloops <- RLOOPS %>% csvToGR() 
  GenomicRanges::mcols(rloops)$keepName <- GenomicRanges::mcols(rloops) %>% rownames()
  rloops <- rloops %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::rename(names = keepName) %>%
    dplyr::mutate(loc = paste0(seqnames, "_", start, "_", end)) %>%
    dplyr::select(-c("seqnames", "start", "end", "width", "strand", "score"))
  
  # Read in the TSV
  rlsm <- readr::read_tsv(RLOOP_SIGNAL_MAT) %>% 
    dplyr::mutate(loc = paste0(`#'chr'`, "_", `'start'` + 1, "_", `'end'`)) %>%
    dplyr::select(-c(`#'chr'`, `'start'`,`'end'`)) %>%
    dplyr::right_join(rloops, by = "loc") %>%
    dplyr::select(-loc) %>%
    dplyr::rename(rloop_id=names) %>%
    tidyr::pivot_longer(cols = ! rloop_id,
                        names_to = "rmap_sample_id",
                        values_to = "norm_counts") %>%
    dplyr::mutate(rmap_sample_id = gsub(rmap_sample_id,
                                        pattern = "'([ES]+RX[0-9]+)_.+",
                                        replacement = "\\1"))
  
  # Combine the rlsm with the rlol
  full_join(rlol, rlsm, by = c("rloop_id", "rmap_sample_id")) %>%
    replace_na(replace = list(numOlap = 0, qVal = 0, signalVal = 0, norm_counts = 0)) %>%
    ungroup()
  
}

# Input
if (! interactive()) {
  hg38rl <- snakemake@input[['rl']]
  rlregions <- snakemake@input[["np"]]
  rlregout <- snakemake@output[['tsv']]
  rlregoutbed <- snakemake@output[['bed']]
  sigtsv <- snakemake@output[['sigtsv']]
  manifest <- snakemake@output[['manifest']]
  
} else {
  hg38rl <- "rlbase-data/rlregions/hg38_manifest_All.csv"
  rlregions <- "rlbase-data/rlregions/rlregions_All.narrowPeak"
  rlregout <- "rlbase-data/rlregions/rlregions_All.csv"
  rlregoutbed <- "rlbase-data/rlregions/rlregions_All.bed"
  sigtsv <- "rlbase-data/rlregions/rlregions_All_signal.tsv"
  manifest <- "rlbase-data/rlbase_manifest_final.tsv"
}

peaks <- read_csv(hg38rl)
manifest <- read_tsv(manifest)
MAX_RL_SIZE <- 50000
makeRLoops(ORIG_RL = rlregions, MAX_RL_SIZE, RLOOPS = rlregout)

# Need to create BED file for deeptools to use
rlgr <- rlregout %>% csvToGR()
rlgr %>% rtracklayer::export.bed(con = rlregoutbed)

# Get overlap for all peaks and the R-loops
message("- Obtaining peak list")
rltbl <- read_csv(rlregout) %>%
  mutate(
    seqnames = gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\1"),
    start = as.numeric(gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\2")),
    end = as.numeric(gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\3")),
    strand = gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\4")
  ) %>%
  select(
    chrom=seqnames, start, end, strand, rlregion=id, confidence_level, is_rlfs
  )
pklst <- pblapply(
    peaks[,1,drop=TRUE],
    function(peak) {
      read_tsv(peak, col_names = c(
        "chrom", "start", "end", "name", "score", "strand", "signalVal", "pVal", "qVal"
      ), show_col_types = FALSE, progress = FALSE)
    }
  )
names(pklst) <- gsub(peaks[,1,drop=TRUE], pattern = ".+/peaks/([ES]{1}RX[0-9]+)_.+", replacement = "\\1")
message("- Intersecting peaks with RLRegions")
isct <- valr::bed_intersect(rltbl, pklst[1:length(pklst)])

# Get the average signal per sample-rlregion combination
message("- Averaging signal within sample-rlregion combinations")
isctSR <- isct %>%
  group_by(rlregion.x, .source) %>%
  mutate(
    signalVal.y = mean(signalVal.y),
    pVal.y = mean(pVal.y),
    qVal.y = mean(qVal.y)
  ) %>%
  distinct(rlregion.x, .source, .keep_all=TRUE) %>%
  ungroup()

# Get average and median pval, signalval, and qval per rlregion
message("- getting average signal per RLRegion")
isctSRa <- isctSR %>%
  group_by(rlregion.x) %>%
  mutate(avgSignalVal = mean(signalVal.y),
         medSignalVal = median(signalVal.y),
         avgPVal = mean(pVal.y),
         medPVal = median(pVal.y),
         avgQVal = mean(qVal.y),
         medQVal = median(qVal.y)) %>%
  ungroup()

# Intersect sample metadata and RL regions
rlregionsOl <- select(isctSRa, rlregion=rlregion.x, sample=.source, 
                      is_rlfs=is_rlfs.x, confidence_level=confidence_level.x, 
                      signalVal=signalVal.y, pVal=pVal.y, qVal=qVal.y,
                      contains("med"), contains("avg"))
message("- Adding in metadata")
manifestToUse <- manifest %>%
  # filter(condType == "POS") %>%
  select(
    sample=experiment, 
    mode, tissue, ip_type, 
    study, condType, verdict, numPeaks
  ) 

rlregions <- rlregionsOl %>%
  inner_join(
      manifestToUse, by = "sample"
  ) %>%
  distinct(rlregion, sample, .keep_all = TRUE) %>%
  group_by(rlregion) %>%
  mutate(
    samples = paste0(sample, collapse = "\n"),
    nSamples = length(samples),
    studies = paste0(unique(study), collapse = "\n"),
    nStudies = length(unique(study)),
    modes = paste0(unique(mode), collapse = "\n"),
    nModes = length(unique(mode)),
    tissues = paste0(unique(tissue), collapse = "\n"),
    nTissues = length(unique(tissue)),
    ip_types = paste0(unique(ip_type), collapse = "\n"),
    nIPTypes = length(unique(ip_type)),
    pct_case = 100*(sum(verdict == "Case") / length(samples)),
    avgNumPeaks = mean(numPeaks),
    medNumPeaks = median(numPeaks)
  ) %>%
  select(-sample, -signalVal, -pVal, -qVal, -tissue, -verdict, 
         -study, -mode, -ip_type, -condType) %>%
  distinct(rlregion, .keep_all = TRUE)
rlregions



# Wrangle intersect with metadata
select(isct, rlregion=name.x, experiment=.source, seqnames=chrom,
       start=start.x, end=end.x, strand=strand.x, confidenceScore=score.x) 

# Make the RLoop Signal File
message("Making RL signal output")
makeRLoopSignal(RLOOPS = RLOOPSXZ, RLOOP_SIGNAL_MAT = RLOOP_SIGNAL_MAT, RLOOP_OLAP = RLOOP_OLAP) %>%
  readr::write_csv(RLOOP_SIGNAL_OUT)

# Compress
message("Compress")







