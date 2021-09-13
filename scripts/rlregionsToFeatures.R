#' This script aligns RLRegions and genomic features
opts <- c("All", "dRNH", "S96")
library(tidyverse)
library(pbapply)
pbo <- pboptions(type="txt") 

# Get annots
annot_pat <- "([a-zA-Z0-9]+)/([a-zA-Z0-9_ \\-]+)\\.csv\\.gz"
annots <- tibble(Key = list.files("misc-data/annotations/", recursive = TRUE)) %>%
  dplyr::filter(grepl(Key, pattern = annot_pat)) %>%
  dplyr::mutate(genome = gsub(Key, pattern = annot_pat, replacement = "\\1"),
                db = gsub(Key, pattern = annot_pat, replacement = "\\2")) %>%
  filter(genome == "hg38")

# Download into R environment as a list
annotations <- pbapply::pblapply(
  seq(nrow(annots)), 
  function(i) {
    fl <- annotGen$Key[i]
    readr::read_csv(
      paste0("misc-data/annotations/", fl), 
      show_col_types = FALSE, 
      progress = FALSE
    ) 
  }
) %>% bind_rows()

lapply(opts, function(opt) {
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
  intsct <- valr::bed_intersect(rlregions, annotations)
  intsct %>% 
    dplyr::select(
      rlregion=name.x, annotation=name.y
    ) %>% 
    write_csv(paste0("rlbase-data/misc/rlregions_", opt, "_annotations.csv"), progress = FALSE)
})



