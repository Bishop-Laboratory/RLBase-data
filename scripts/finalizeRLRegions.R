#' This script will wrangle the results of the RLRegions pipeline to create final tables

library(tidyverse)
library(pbapply)
pbo <- pboptions(type="txt") 

#' Helper function for converting CSV to GR
csvToGR <- function(csv) {
  readr::read_csv(csv, show_col_types = FALSE) %>%
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
    GenomicRanges::makeGRangesFromDataFrame() 
}

#' Makes the "Rloops" table
makeRLoops <- function(ORIG_RL, MAX_RL_SIZE, RLOOPS) {
  
  cond <- gsub(RLOOPS, pattern = ".*rlregions_([a-zA-Z0-9]+)\\.csv", replacement = "\\1")
  
  readRL <- function(x, max_rl_size) {
    readr::read_tsv(x, col_names = c(
      "chrom", "start", "end", "name", "score", "strand", "signalValue", "pval", "qval", "peak"
    ),skip = 1, show_col_types = FALSE) %>%
      dplyr::filter(end - start < max_rl_size) %>%
      dplyr::mutate(group = {{ cond }}) %>%
      return()
  }
  
  # Get rloops
  rloops <- ORIG_RL %>%
    readRL(MAX_RL_SIZE) %>%
    rownames_to_column(var = "id") %>%
    mutate(id = paste0(group, "_RL", id),
           type = "Unknown",
           location = paste0(chrom, ":", start, "-", end, ":", strand),
           confidence_level = score,
           origName = name) 
  
    
  # Check RLFS, chrom_sizes, and mask
  RLFS <- rtracklayer::import.bed(
    file.path(
      RLSeq:::RLBASE_URL,
      "rlfs-beds",
      paste0("hg38", ".rlfs.bed")
    )
  )
  
  # Prevent stranded assignment
  GenomicRanges::strand(RLFS) <- "*"
  RLFS <- GenomicRanges::reduce(RLFS)
  
  # Make this a tbl
  rlfs <- RLFS %>%
    as_tibble() %>%
    dplyr::rename(chrom=seqnames)
  
  # Intersect
  olr <- valr::bed_intersect(rloops, rlfs) %>%
    pull(id.x) %>% 
    unique()
  
  source <- gsub(RLOOPS, pattern = ".*rlregions_(.+)\\.csv", replacement = "\\1")
  if (source == "All") {
    sources <- rloops$signalValue
  } else {
    sources <- source
  }
  rlregout <- rloops %>%
    mutate(is_rlfs = id %in% !! olr,
           source = {{ sources }}) %>%
    dplyr::select(id, type, location, confidence_level, is_rlfs, origName, source, type)
  
  write_csv(rlregout, RLOOPS)
  return(NULL)
}

# Input
if (! interactive()) {
  hg38rl <- snakemake@input[['rl']]
  rlregionsIn <- snakemake@input[["np"]]
  rlregout <- snakemake@output[['csv']]
  rlregtbl <- snakemake@output[['tbl']]
  rlregoutbed <- snakemake@output[['bed']]
  sigtsv <- snakemake@output[['sigtsv']]
  manifest <- snakemake@params[['manifest']]
  blacklist <- snakemake@params[['blacklist']]
} else {
  hg38rl <- "rlregions/hg38_manifest_dRNH.csv"
  rlregionsIn <- "rlregions/rlregions_dRNH.narrowPeak"
  rlregout <- "rlregions/rlregions_dRNH.csv"
  rlregtbl <- "misc/rlregions_dRNH_table.tsv"
  rlregoutbed <- "rlregions/rlregions_dRNH.bed"
  sigtsv <- "rlregions/rlregions_dRNH_signal.tsv"
  manifest <- "rlbase_samples.tsv"
  blacklist <- "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
  setwd("rlbase-data/")
}

message(hg38rl)
message(rlregionsIn)
message(rlregout)
message(rlregtbl)
message(rlregoutbed)
message(sigtsv)
message(manifest)
message(blacklist)
message(getwd())

peaks <- read_csv(hg38rl, show_col_types = FALSE)
manifest <- read_tsv(manifest, show_col_types = FALSE)
MAX_RL_SIZE <- 50000

message("- Overlapping with RLFS")
a_ <- makeRLoops(ORIG_RL = rlregionsIn, MAX_RL_SIZE, RLOOPS = rlregout)

# Need to create BED file for deeptools to use
rlgr <- rlregout %>% csvToGR()
rlgr %>% rtracklayer::export.bed(con = rlregoutbed)

# Get overlap for all peaks and the R-loops
message("- Obtaining peak list")
rltbl <- read_csv(rlregout, show_col_types = FALSE) %>%
  mutate(
    seqnames = gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\1"),
    start = as.numeric(gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\2")),
    end = as.numeric(gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\3")),
    strand = gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\4")
  ) %>%
  dplyr::select(
    chrom=seqnames, start, end, strand, rlregion=id, confidence_level, is_rlfs, source
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
pklst <- pklst[sapply(pklst, function(x) {nrow(x) > 0})]
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

# Get the locations
loc <- isctSRa %>%
  mutate(location = paste0(chrom, ":", start.x, "-", end.x, ":", strand.x)) %>%
  dplyr::select(
    rlregion=rlregion.x,
    location
  ) %>% distinct(rlregion, location)

# Intersect sample metadata and RL regions
rlregionsOl <- dplyr::select(isctSRa, rlregion=rlregion.x, sample=.source, 
                             is_rlfs=is_rlfs.x, 
                             source = source.x,
                             confidence_level=confidence_level.x, 
                             signalVal=signalVal.y,
                             pVal=pVal.y, qVal=qVal.y,
                      contains("med"), contains("avg")) %>%
  left_join(loc) %>%
  relocate(location, .after = rlregion)
message("- Adding in metadata")
manifestToUse <- manifest %>%
  dplyr::filter(label == "POS") %>%
  dplyr::select(
    sample=experiment, 
    mode, tissue, ip_type, 
    study, label, prediction, numPeaks
  ) 
rlregions <- rlregionsOl %>%
  inner_join(
      manifestToUse, by = "sample"
  ) %>%
  distinct(rlregion, sample, .keep_all = TRUE) %>%
  group_by(rlregion) %>%
  mutate(
    samples = paste0(sample, collapse = ","),
    nSamples = length(samples),
    studies = paste0(unique(study), collapse = ","),
    nStudies = length(unique(study)),
    modes = paste0(unique(mode), collapse = ","),
    nModes = length(unique(mode)),
    tissues = paste0(unique(tissue), collapse = ","),
    nTissues = length(unique(tissue)),
    ip_types = paste0(unique(ip_type), collapse = ","),
    nIPTypes = length(unique(ip_type)),
    pct_case = 100*(sum(prediction == "POS") / length(samples)),
    avgNumPeaks = mean(numPeaks),
    medNumPeaks = median(numPeaks)
  ) %>%
  dplyr::select(-sample, -signalVal, -pVal, -qVal, -tissue, -prediction, 
         -study, -mode, -ip_type, -label) %>%
  distinct(rlregion, .keep_all = TRUE)
rlregions <- rlregions %>% 
  arrange(desc(nStudies), desc(avgQVal)) 

## Add in genes ##
message("- Adding in genes")

# Get latest gene symbols
genes <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86 %>%
  AnnotationDbi::select(
    keys = AnnotationDbi::keys(.), 
    columns = c("SYMBOL", "GENEBIOTYPE", 
                "SEQNAME", "GENESEQSTART", 
                "GENESEQEND", "SEQSTRAND")
  ) %>% 
  dplyr::select(
    geneid = GENEID,
    symbol = SYMBOL,
    chrom = SEQNAME,
    start = GENESEQSTART,
    end = GENESEQEND,
    strand = SEQSTRAND
  ) %>%
  mutate(strand = case_when(strand == 1 ~ "+",
                            strand == -1 ~ "-",
                            TRUE ~ "*"),
         chrom = paste0("chr", chrom))

# Convert rlregions for compatibility with valr
rlregionsSplit <- rlregions %>%
  mutate(
    chrom = gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\1"),
    start = as.numeric(gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\2")),
    end = as.numeric(gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\3")),
    strand = gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\4")
  ) %>%
  dplyr::select(
    chrom, start, end, strand, rlregion
  )

# Intersect
geneol <- valr::bed_intersect(rlregionsSplit, genes) %>%
  dplyr::select(rlregion=rlregion.x, geneIDs=geneid.y, allGenes=symbol.y)

# Get commonly-annotated genes from pathway databases
mainGeneVec <- pblapply(
    c("C2", "C5", "H", "C8"), 
    msigdbr::msigdbr,
    species = "Homo sapiens"
  ) %>%
  bind_rows() %>%
  pull(gene_symbol) %>%
  unique()

# Main vs supplemental genes
geneol <- geneol %>%
  mutate(
    mainGenes = ifelse(allGenes %in% !! mainGeneVec, allGenes, "")
  ) %>%
  group_by(
    rlregion
  ) %>%
  mutate(
    geneIDs = paste0(geneIDs, collapse = ","),
    allGenes = paste0(allGenes, collapse = ","),
    mainGenes = paste0(mainGenes[mainGenes != ""], collapse = ",")
  ) %>%
  ungroup() %>%
  distinct(rlregion, .keep_all = TRUE)

# Add back to main table
if ("mainGenes" %in% colnames(rlregions)) {
  rlregions <- select(rlregions, -contains("gene"))
}
rlregions <- left_join(
  rlregions, geneol, by = "rlregion"
) %>%
  mutate(across(contains("gene"), function(x) {
    ifelse(is.na(x), "", x)
  })) 

## Get the repeat regions & blacklists ##

# Get blacklist 
blklst <- read_tsv(blacklist, col_names = c("chrom", "start", "end"), show_col_types = FALSE) %>%
  mutate(strand = "*",
         type = "blacklist")
# Get repeats
annos <- RLHub::annots_full_hg38()
rmsk <- annos[grepl(names(annos), pattern = "Repeat_Masker.+|Centromeres+")] %>%
  bind_rows(.id = "type")
rmsk <- rmsk %>%
  mutate(type = gsub(type, pattern = ".+__(.+)", replacement = "\\1")) %>%
  filter(type %in% c("Centromeres", "Retroposon", "rRNA", "Satellite", 
                      "scRNA", "snRNA", "tRNA", "RC", "srpRNA"))
# Bind together
repeats <- bind_rows(blklst, rmsk)

# Add within 1MB of centromeres to censor list
repeats <- repeats %>%
  mutate(
    start = case_when(
      type == "Centromeres" ~ start - 0.5E6,
      TRUE ~ start
    ),
    end = case_when(
      type == "Centromeres" ~ end + 0.5E6,
      TRUE ~ end
    )
  ) 

# Intersect
repol <- valr::bed_intersect(rlregionsSplit, repeats) %>%
  dplyr::select(rlregion=rlregion.x) %>%
  pull(rlregion)

# Modify the rlregions table
common_chrom <- paste0("chr", c(1:22, "X", "Y", "M"))
rlregions <- rlregions %>%
  mutate(
    is_repeat = rlregion %in% repol,
    chrom = gsub(location, pattern = "(chr[A-Z0-9]+):.+", replacement = "\\1")
  ) %>%
  dplyr::filter(
    chrom %in% {{ common_chrom }}
  )

## Save outputs ##
write_tsv(rlregions, file = rlregtbl)
write_tsv(rlregionsOl, file = sigtsv)

message("Done")
