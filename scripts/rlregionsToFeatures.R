#' This script aligns RLRegions and genomic features
opts <- c("All", "dRNH", "S96")
library(tidyverse)
library(pbapply)
pbo <- pboptions(type="txt") 

# Get annots
load("misc-data/annotations/annotations_hg38.rda")
annotationsTbl <- bind_rows(annotations)

a_ <- pblapply(opts, function(opt) {
  message(opt)
  
  # Get rlr
  fl <- paste0("rlbase-data/misc/rlregions_", opt, "_table.tsv")
  locpattern <- "(.+):(.+)\\-(.+):(.+)"
  rlregions <- read_tsv(fl, show_col_types = FALSE, progress = FALSE)
  rlregions <- rlregions %>% 
    mutate(chrom = gsub(location, pattern = locpattern, replacement = "\\1"),
           start = as.numeric(gsub(location, pattern = locpattern, replacement = "\\2")),
           end = as.numeric(gsub(location, pattern = locpattern, replacement = "\\3")),
           strand = gsub(location, pattern = locpattern, replacement = "\\4")) %>%
    dplyr::select(chrom, start, end, strand, name=rlregion)
  
  # Overlap
  intsct <- valr::bed_intersect(rlregions, annotationsTbl)
  intsct %>% 
    dplyr::select(
      rlregion=name.x, annotation=name.y
    ) %>% 
    write_csv(paste0("rlbase-data/misc/rlregions_", opt, "_annotations.csv"), progress = FALSE)
})

message("Done")

