library(tidyverse)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(parallel)

# Tables to query along with column name in each table that contains "type"
# TODO: Acc U2 snRNP sites and Poly-A Signal (PAS) sites
# TODO: Add G4Q sites 
# TODO: Add SkewR-predicted regions (or G/C-skew regions)
# See this ref https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC6125637/
tablesToUse <- data.frame(
  row.names = c("table", "typeCol"),
  "CpG_Islands"= c("cpgIslandExt", "typeNow"),
  "Centromeres" = c("centromeres", "typeNow"),
  "Encode_CREs" = c("encodeCcreCombined", "ucscLabel"),
  "knownGene_RNAs" = c("knownGene", "transcriptType"),
  "Microsatellite" = c("microsat", "typeNow"),
  "Repeat_Masker" = c("rmsk", "repClass"),
  "Splice_Events" = c("knownAlt", "name"),
  "snoRNA_miRNA_scaRNA" = c("wgRna", "type"),
  "tRNAs" = c("tRNAs", "typeNow")
) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  as_tibble()

# TODO: Remove Problematic annotations
# problems <- c("Mt_rRNA", "Mt_tRNA")

# Get the genes in a list
genomes <- list("hg38" = TxDb.Hsapiens.UCSC.hg38.knownGene,
                "mm10" = TxDb.Mmusculus.UCSC.mm10.knownGene)

# For each genome, retrieve all annotations
annotationLst <- lapply(names(genomes), function(genome) {
  
  message(genome)
  
  # Get the tables from UCSC ** LONG TUNNING
  message("Getting tables from UCSC...")
  dir.create("tmp", showWarnings = FALSE)
  TMPFILE <- paste0("tmp/tabLst_", genome, ".rda")
  if (! file.exists(TMPFILE)) {
    tabLst <- mclapply(
      pull(tablesToUse, table), function(tabNow) {
        message(tabNow)
        ucscTableQuery(genome, table = tabNow) %>%
          getTable()
      }, mc.cores = length(pull(tablesToUse, table)) 
    )
    names(tabLst) <- pull(tablesToUse, group)
    save(tabLst, file = TMPFILE)
  } else {
    load(TMPFILE)
  }
  
  # Remove any which are not available
  message("Wrangling...")
  errors <- sapply(tabLst, function(x) {
    class(x) == "try-error"
  })
  if (any(errors)) {
    tabLst <- tabLst[-which(errors)]
  }
  
  # Add the type information
  tabLst2 <- lapply(names(tabLst), function(nm) {
    x <- tabLst[[nm]]
    x$typeNow <- nm
    x
  })
  names(tabLst2) <- names(tabLst)
  
  # Perform the group assignment and clean column names
  tabDF <- lapply(names(tabLst2), function(nm) {
    message(nm)
    x <- tabLst2[[nm]]
    tosplit <- tablesToUse %>%
      dplyr::filter(group == !! nm) %>%
      pull(typeCol)
    
    if (! "strand" %in% colnames(x)) {
      x$strand <- "."
    }
    
    x <- x %>%
      as_tibble() %>%
      dplyr::mutate(strand = case_when(strand == "." ~ "*",
                                       TRUE ~ strand)) %>%
      {
        if (! tosplit %in% colnames(.)) { 
          dplyr::mutate(., !!quo_name(tosplit) := typeNow) 
        } else {
          .
        }
      } %>%
      dplyr::select(contains(c("chrom", "geno", "tx")) & ! contains(c("Left", "Starts", "Url")), 
                    strand,
                    group = !! tosplit,
                    DB = typeNow) 
    
    if ("genoName" %in% colnames(x)) {
      x <- x %>% dplyr::rename(
        chrom=genoName, chromStart = genoStart, chromEnd = genoEnd
      )
    } 
    
    if ("txStart" %in% colnames(x)) {
      x <- x %>% dplyr::rename(
        chromStart = txStart, chromEnd = txEnd
      )
    } 
    x
  }) %>% bind_rows()
  
  # Annotations to keep  
  keepers <- table(tabDF$group) %>%
    as.data.frame() %>%
    filter(Freq > 1000 | Var1 %in% c("tRNA", "Centromeres", "srpRNA", "snoRNA",
                                     "rRNA", "scRNA", "origin_of_replication",
                                     "non_stop_decay", "Mt_rRNA", "Mt_tRNA",
                                     "RC", "ribozyme", "scaRNA"),
           ! Var1 %in% c("retained_intron", "processed_transcript", 
                         "knownGene_RNAs", "Unknown", "Other",
                         "misc_RNA", "TEC", "enhancer")) %>%
    pull(Var1)
  
  # Get the final annotations
  tabDF2 <- tabDF %>%
    dplyr::filter(group %in% keepers) %>%
    dplyr::rename(seqnames = chrom, start = chromStart, end = chromEnd)
  
  # Get the TxDb annotations
  message("Getting TxDb annotations...")
  txdb <- genomes[[genome]]
  # TODO: Clean this part up
  fiveUTR <- GenomicFeatures::fiveUTRsByTranscript(txdb) %>%
    unlist() %>%
    as.data.frame(row.names = seq(names(.))) %>%
    mutate(group = "fiveUTR") %>%
    dplyr::select(seqnames, start, end , strand, group)
  threeUTR <- GenomicFeatures::threeUTRsByTranscript(txdb) %>%
    unlist() %>%
    as.data.frame(row.names = seq(names(.))) %>%
    mutate(group = "threeUTR") %>%
    dplyr::select(seqnames, start, end, strand, group)
  exon <- GenomicFeatures::tidyExons(txdb) %>%
    as.data.frame() %>%
    mutate(group = "Exon") %>%
    dplyr::select(seqnames, start, end, strand, group) %>%
    distinct()
  intron <- GenomicFeatures::tidyIntrons(txdb) %>%
    as.data.frame() %>%
    mutate(group = "Intron") %>%
    dplyr::select(seqnames, start, end, strand, group) %>%
    distinct()
  TSS <- GenomicFeatures::tidyTranscripts(txdb) %>%
    as.data.frame() %>%
    mutate(group = "TSS",
           end = start + 1) %>%
    dplyr::select(seqnames, start, end, strand, group) %>%
    distinct()
  TTS <- GenomicFeatures::tidyTranscripts(txdb) %>%
    as.data.frame() %>%
    mutate(group = "TTS",
           start = end - 1) %>%
    dplyr::select(seqnames, start, end, strand, group) %>%
    distinct()
  Intergenic <- GenomicFeatures::tidyTranscripts(txdb) %>%
    gaps() %>%
    as.data.frame() %>%
    mutate(group = "Intergenic") %>%
    dplyr::select(seqnames, start, end, strand, group) %>%
    distinct()
  genomicAnno <- bind_rows(list(
    TSS, fiveUTR, exon, intron, threeUTR, TTS
  ))
  
  # Add to the annotation table
  message("Compiling...")
  bind_rows(tabDF2, genomicAnno) %>%
    dplyr::rename(type = group) %>%
    dplyr::select(-DB) %>%
    dplyr::group_by(type) %>%
    {setNames(group_split(.), group_keys(.)[[1]])} 
  
})
names(annotationLst) <- names(genomes)  # Add names back

# Save
usethis::use_data(annotationLst,
                  overwrite = TRUE, 
                  compress = "xz",
                  internal = FALSE)
