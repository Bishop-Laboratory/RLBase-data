#' Makes the "Rloops" table
makeRLoops <- function(ORIG_RL, MAX_RL_SIZE, RLOOPS) {
  require(ChIPpeakAnno)
  require(tidyverse)
  
  urlExists <- function(url) {
    identical(
      httr::status_code(
        # Checks HEAD only due to size constraints
        httr::HEAD(
          url
        )
      ), 200L  # Checks if response is ok
    )
  }
  
  readRL <- function(x, max_rl_size) {
    readr::read_tsv(x, col_names = c(
      "seqnames", "start", "end", "name", "score", "strand", "signalValue", "pval", "qval", "peak"
    ),skip = 1) %>%
      dplyr::filter(end - start < max_rl_size) %>%
      as.data.frame() %>%
      regioneR::toGRanges() %>%
      return()
  }
  checkRLFSAnno <- function(genome) {
    return(
      urlExists(
        paste0(
          "https://rmapdb-data.s3.us-east-2.amazonaws.com/rlfs-beds/", 
          genome, ".rlfs.bed"
        )
      )
    )
  }
  
  getRLFSAnno <- function(genome) {
    
    # Check if annotations available first
    if (! checkRLFSAnno(genome)) {
      stop("No RLFS annotations available for ", genome)
    }
    
    # Return as a GRanges object
    return(
      regioneR::toGRanges(
        as.data.frame(
          suppressMessages(readr::read_tsv(
            paste0(
              "https://rmapdb-data.s3.us-east-2.amazonaws.com/rlfs-beds/", 
              genome, ".rlfs.bed"
            ),
            col_names = FALSE))
        )
      )
    )
  }
  
  # Get rloops
  rloops <- ORIG_RL %>%
    readRL(MAX_RL_SIZE)
  
  # Check RLFS, chrom_sizes, and mask
  RLFS <- getRLFSAnno("hg38")
  
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
    mutate(id = paste0("RL", id),
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
  system(paste0("xz -f ", RLOOPS))
  return(NULL)
}

#' Makes the "RLoopSignal" table for RMapDB
#'
#' @importFrom magrittr %>%
#' @import rlang
makeRLoopSignal <- function(RLOOPS, RLOOP_SIGNAL_MAT, RLOOP_OLAP) {
  
  # Get RLoop to peak mapping
  rloops <- read_csv(RLOOPS) %>%
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


if (sys.nframe() == 0L | TRUE) {
  require(rlang)
  require(magrittr)
  source("analyses/Prepare-RMapDB-Tables/GeneRLoopOverlap_GenFeatRLoopOverlap/gene_rl_overlap__gf_rl_overlap.R")
  
  # TODO: Needs to become part of the snake pipe
  RLOOPS <- "analyses/Prepare-RMapDB-Tables/rloops.csv"
  RLOOPSXZ <- paste0(RLOOPS, ".xz")
  RMAPSAMPS <- "analyses/Prepare-RMapDB-Tables/rmap_samples.csv.xz"
  RLOOP_SIGNAL_MAT <- "analyses/Prepare-RMapDB-Tables/RLoops/data/rloop_signal.tsv.xz"
  RLOOP_SIGNAL_OUT <- "analyses/Prepare-RMapDB-Tables/rloop_signal.csv"
  RLOOP_OLAP <- "analyses/Prepare-RMapDB-Tables/RLoops/data/rloops/strategy.a__10bp__peaks.countMat.tsv.xz"
  ORIG_RL <- "analyses/Prepare-RMapDB-Tables/RLoops/data/rloops/strategy.a__10bp__peaks.narrowPeak"
  MAX_RL_SIZE <- 50000
  
  makeRLoops(ORIG_RL, MAX_RL_SIZE, RLOOPS)
  
  # Need to create BED file for deeptools to use
  RLOOPS_BED <- "analyses/Prepare-RMapDB-Tables/RLoops/data/rloops.bed"
  RLOOPSXZ %>% csvToGR() %>% rtracklayer::export.bed(con = RLOOPS_BED)
  system(paste0("xz -f ", RLOOPS_BED))
  
  ### IMPORTANT: NEED TO RUN DEEPTOOLS HERE -- USE SNAKEMAKE ###
  
  # Make the RLoop Signal File
  message("Making RL signal output")
  makeRLoopSignal(RLOOPS = RLOOPSXZ, RLOOP_SIGNAL_MAT = RLOOP_SIGNAL_MAT, RLOOP_OLAP = RLOOP_OLAP) %>%
    readr::write_csv(RLOOP_SIGNAL_OUT)
  
  # Compress
  message("Compress")
  system(paste0("xz -f ", RLOOP_SIGNAL_OUT))  
  
}

