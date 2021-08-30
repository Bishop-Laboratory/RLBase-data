#' Makes the overlap tables for RMapDB
#'
#' @importFrom magrittr %>%
#' @import rlang
makeOverlapTable <- function(
  peakLst
) {
  # TODO: Should we consider a separate category for promoter overlap?
  ChIPpeakAnno::findOverlapsOfPeaks(peakLst, ignore.strand = FALSE) %>%
    purrr::pluck("overlappingPeaks") %>%
    purrr::pluck(1) %>%
    dplyr::select(
      !! names(peakLst)[1] := peaks1,
      !! names(peakLst)[2] := peaks2,
    )
}

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


# From https://stackoverflow.com/questions/47932246/rscript-detect-if-r-script-is-being-called-sourced-from-another-script/47932989
if (sys.nframe() == 0L) {
  
  RLOOPS <- "analyses/Prepare-RMapDB-Tables/rloops.csv.xz"
  GENES <- "analyses/Prepare-RMapDB-Tables/genes.csv.xz" 
  GENFEATS <- "analyses/Prepare-RMapDB-Tables/genomic_features.csv.xz"
  GENE_RLOOP_OVERLAP <- "analyses/Prepare-RMapDB-Tables/gene_rl_overlap.csv"
  GENFEAT_RLOOP_OVERLAP <- "analyses/Prepare-RMapDB-Tables/gf_rl_overlap.csv"
  
  require(rlang)
  require(magrittr)

  # Process GENE_RLOOP_OVERLAP
  message("Starting Genes")
  peakLst <- list(
    rloop_id = csvToGR(RLOOPS),
    gene_id = csvToGR(GENES)
  ) %>%
    makeOverlapTable() %>%
    readr::write_csv(file = GENE_RLOOP_OVERLAP)
  
  # Process GENFEAT_RLOOP_OVERLAP
  message("Starting GenFeat")
  peakLst <- list(
    rloop_id = csvToGR(RLOOPS),
    feature_id = csvToGR(GENFEATS)
  ) %>%
    makeOverlapTable() %>%
    readr::write_csv(file = GENFEAT_RLOOP_OVERLAP) 
  
  # xzip Compress
  message("Compress")
  system(paste0("xz -f ", GENE_RLOOP_OVERLAP))  
  system(paste0("xz -f ", GENFEAT_RLOOP_OVERLAP))
  
}

