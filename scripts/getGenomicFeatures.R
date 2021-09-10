library(tidyverse)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(parallel)
library(pbapply)
pbo <- pboptions(type="txt") 

dir.create("misc-data/annotations", showWarnings = FALSE)

### Get genomic features from UCSC ###

# Get the additional tables from UCSC
# "typeNow" indicates no column to stratify by
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
  "tRNAs" = c("tRNAs", "typeNow"),
  "CNA" = c("coriellDelDup", "CN_State"),
  "PolyA" = c("wgEncodeGencodePolyaV38", "name2")
) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  as_tibble()

# Only human for now
# Get the genes in a list
# TODO: Need to remove LRGs from all ENsDB
# TODO: Need to choose genes from MSigDB based off the ens gene ID, not symbols...
# TODO: Need to use the GTF files instead of the EnsDB packages
genomes <- list("hg38" = GenomicFeatures::makeTxDbFromUCSC(genome = "hg38"),
                "mm10" = GenomicFeatures::makeTxDbFromUCSC(genome = "mm10"))

# For each genome, retrieve all annotations
annotationLst <- lapply(names(genomes), function(genome) {
  
  message(genome)
  
  # Get the tables from UCSC ** LONG TUNNING
  message("Getting tables from UCSC...")
  dir.create("tmp", showWarnings = FALSE)
  TMPFILE <- paste0("tmp/tabLst_", genome, ".rda")
  # Get the tables from UCSC ** LONG TUNNING
  tabLst <- pblapply(
    pull(tablesToUse, table), function(tabNow, genome) {
      library(magrittr)
      message(tabNow)
      try(
        rtracklayer::ucscTableQuery(x=genome, table = tabNow) %>%
          rtracklayer::getTable()
      ) 
    },
    genome=genome
  )
  
  # Remove any which are not available
  message("Wrangling...")
  errors <- sapply(tabLst, function(x) {
    class(x) == "try-error"
  })
  tabNow <- tablesToUse
  if (any(errors)) {
    tabLst <- tabLst[-which(errors)]
    tabNow <- tabNow[-which(errors),]
  }
  names(tabLst) <- pull(tabNow, group)
  
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
    tosplit <- tabNow %>%
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
                    db = typeNow) 
    
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
    filter(! Var1 %in% c("retained_intron", "processed_transcript", 
                         "knownGene_RNAs", "Unknown", "Other",
                         "LINE?", "DNA?", "RC?", "LTR?", "SINE?",
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
    dplyr::select(seqnames, start, end , strand, group) %>% 
    mutate(db="Transcript Features")
  threeUTR <- GenomicFeatures::threeUTRsByTranscript(txdb) %>%
    unlist() %>%
    as.data.frame(row.names = seq(names(.))) %>%
    mutate(group = "threeUTR") %>%
    dplyr::select(seqnames, start, end, strand, group) %>% 
    mutate(db="Transcript Features")
  exon <- GenomicFeatures::tidyExons(txdb) %>%
    as.data.frame() %>%
    mutate(group = "Exon") %>%
    dplyr::select(seqnames, start, end, strand, group) %>% 
    mutate(db="Transcript Features") %>%
    distinct()
  intron <- GenomicFeatures::tidyIntrons(txdb) %>%
    as.data.frame() %>%
    mutate(group = "Intron") %>%
    dplyr::select(seqnames, start, end, strand, group) %>% 
    mutate(db="Transcript Features") %>%
    distinct()
  TSS <- GenomicFeatures::tidyTranscripts(txdb) %>%
    as.data.frame() %>%
    mutate(group = "TSS",
           end = start + 1) %>%
    dplyr::select(seqnames, start, end, strand, group) %>% 
    mutate(db="Transcript Features") %>%
    distinct()
  TTS <- GenomicFeatures::tidyTranscripts(txdb) %>%
    as.data.frame() %>%
    mutate(group = "TTS",
           start = end - 1) %>%
    dplyr::select(seqnames, start, end, strand, group) %>% 
    mutate(db="Transcript Features") %>%
    distinct()
  Intergenic <- GenomicFeatures::tidyTranscripts(txdb) %>%
    gaps() %>%
    as.data.frame() %>%
    mutate(group = "Intergenic") %>%
    dplyr::select(seqnames, start, end, strand, group) %>% 
    mutate(db="Transcript Features") %>%
    distinct()
  genomicAnno <- bind_rows(list(
    TSS, fiveUTR, exon, intron, threeUTR, TTS
  ))
  
  # Add to the annotation table
  message("Compiling...")
  bind_rows(tabDF2, genomicAnno) %>%
    dplyr::rename(type = group) %>%
    dplyr::group_by(db) %>%
    {setNames(group_split(.), group_keys(.)[[1]])} 
  
})
names(annotationLst) <- names(genomes)  # Add names back

### Get DNase and TFBS ###

# DNase
dnase <- read_tsv("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/wgEncodeRegDnaseClustered.txt.gz",
                  col_names = c("x", "chrom", "start", "end")) %>%
  mutate(type = "DNaseHS", db=type, strand="*",) %>%
  dplyr::select(chrom, start, end, strand, type, db)

# TFBS
tfbs <- read_tsv("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/encRegTfbsClustered.txt.gz",
                  col_names = c("x", "chrom", "start", "end", "type")) %>%
  mutate(db="encodeTFBS", strand="*",) %>%
  dplyr::select(chrom, start, end, strand, type, db)
  
# 

### Get G4Q DNA ###

# Data was obtained from: https://doi.org/10.6084/m9.figshare.c.3498270.v1
# It contains predicted G4Q sites in bed files
download.file(
  "https://figshare.com/ndownloader/files/6432597",
  destfile = "misc-data/G4Q.tar.gz"
)
untar("misc-data/G4Q.tar.gz", exdir = "misc-data/")
file.rename("misc-data/DATA/", "misc-data/G4Q")

## Parse G4Q pred sites ##
# The cutoffs for labeling were determined adhoc from examintation of quantiles
# They attempt to bin the data such that interesting biology will be retained
# but dimensionality will be reduced.
# Refer to the paper for the meaning of the terms
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0165101
g4qpred <- read_tsv("misc-data/G4Q/hg38_QFS_ALL.txt", col_names = c("chrom", "loc", "G4Pred", "seq", "strand")) %>%
  mutate(
    start = as.numeric(gsub(loc, pattern = "(.+)-(.+)", replacement = "\\1")),
    end = as.numeric(gsub(loc, pattern = "(.+)-(.+)", replacement = "\\2")),
    tractlen = as.numeric(gsub(G4Pred, pattern = "(.+):(.+):(.+)", replacement = "\\1")),
    numloc = as.numeric(gsub(G4Pred, pattern = "(.+):(.+):(.+)", replacement = "\\2")),
    g4num = as.numeric(gsub(G4Pred, pattern = "(.+):(.+):(.+)", replacement = "\\3")),
    tractlen = case_when(
      tractlen > 12 ~ "13+",
      tractlen > 5 ~ "6-12",
      TRUE ~ "4-5"
    ),
    tractlen = paste0("tl:", tractlen),
    numloc = case_when(
      numloc > 9 ~ "10+",
      numloc > 2 ~ "3-9",
      TRUE ~ "1-2"
    ),
    numloc = paste0("nl:", numloc),
    g4num = case_when(
      g4num > 4 ~ "5+",
      g4num > 1 ~ "2-4",
      TRUE ~ "1"
    ),
    g4num = paste0("gn:", g4num),
    G4Pred = paste0("G4Pred__", tractlen, "_", numloc, "_", g4num)
  ) %>%
  dplyr::select(chrom, start, end, strand, G4Pred)

## Get the experimentally-obtained G4Q datasets
plus <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874_Na_K_PDS_plus_hits_intersect.bed.gz"
minus <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874_Na_K_PDS_minus_hits_intersect.bed.gz"
g4qexpp <- read_tsv(plus, col_names = c("chrom", "start", "end")) %>% mutate(strand = "+")
g4qexpm <- read_tsv(minus, col_names = c("chrom", "start", "end")) %>% mutate(strand = "-")
g4qexp <- bind_rows(g4qexpp, g4qexpm) %>% mutate(type = "G4Q_PDS_NaK_GSE63874")

## Save annotations


### Get Transcription Factors Meta-Clusters from GTRD ###

# Wrangle factor names into uniprot IDs
factorsToGet <- c("NFE2L2", "BRCA1", "BRCA2", "STAG2", "STAG1", "RAD23A",
                  "SMC1A", "SLC3A2", "SLC7A11", "EZH2", "HDAC1", "HDAC2", "HDAC3",
                  "TEAD1", "CTCF", "SRSF2", "SRSF1", "U2AF1", "SF3B1", "DHX9")



### Get annotations from other sources ###

# 2. Histone modifications
# Selected from the literature:
# https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC6259673/
# https://www.frontiersin.org/articles/10.3389/fgene.2020.00043/full
# https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC5123964/
histList <- list(
  "H3K27ac" = c(
    "https://www.encodeproject.org/files/ENCFF805FFP/@@download/ENCFF805FFP.bed.gz",
    "https://www.encodeproject.org/files/ENCFF518VCX/@@download/ENCFF518VCX.bed.gz",
    "https://www.encodeproject.org/files/ENCFF926NKP/@@download/ENCFF926NKP.bed.gz",
    "https://www.encodeproject.org/files/ENCFF530TUK/@@download/ENCFF530TUK.bed.gz",
    "https://www.encodeproject.org/files/ENCFF023LTU/@@download/ENCFF023LTU.bed.gz",
    "https://www.encodeproject.org/files/ENCFF035VXU/@@download/ENCFF035VXU.bed.gz",
    "https://www.encodeproject.org/files/ENCFF317QGQ/@@download/ENCFF317QGQ.bed.gz",
    "https://www.encodeproject.org/files/ENCFF553MQN/@@download/ENCFF553MQN.bed.gz",
    "https://www.encodeproject.org/files/ENCFF922LHZ/@@download/ENCFF922LHZ.bed.gz",
    "https://www.encodeproject.org/files/ENCFF805GNH/@@download/ENCFF805GNH.bed.gz",
    "https://www.encodeproject.org/files/ENCFF437OIO/@@download/ENCFF437OIO.bed.gz",
    "https://www.encodeproject.org/files/ENCFF757HGV/@@download/ENCFF757HGV.bed.gz",
    "https://www.encodeproject.org/files/ENCFF405SDB/@@download/ENCFF405SDB.bed.gz",
    "https://www.encodeproject.org/files/ENCFF899XEF/@@download/ENCFF899XEF.bed.gz",
    "https://www.encodeproject.org/files/ENCFF831XSS/@@download/ENCFF831XSS.bed.gz",
    "https://www.encodeproject.org/files/ENCFF371BKL/@@download/ENCFF371BKL.bed.gz",
    "https://www.encodeproject.org/files/ENCFF340KSH/@@download/ENCFF340KSH.bed.gz",
    "https://www.encodeproject.org/files/ENCFF097LTN/@@download/ENCFF097LTN.bed.gz"
  ),
  "H3K9ac" = c(
    "https://www.encodeproject.org/files/ENCFF436XTS/@@download/ENCFF436XTS.bed.gz",
    "https://www.encodeproject.org/files/ENCFF069KAG/@@download/ENCFF069KAG.bed.gz",
    "https://www.encodeproject.org/files/ENCFF306KRD/@@download/ENCFF306KRD.bed.gz",
    "https://www.encodeproject.org/files/ENCFF105NDA/@@download/ENCFF105NDA.bed.gz",
    "https://www.encodeproject.org/files/ENCFF350BNS/@@download/ENCFF350BNS.bed.gz",
    "https://www.encodeproject.org/files/ENCFF044QMC/@@download/ENCFF044QMC.bed.gz",
    "https://www.encodeproject.org/files/ENCFF361HMC/@@download/ENCFF361HMC.bed.gz",
    "https://www.encodeproject.org/files/ENCFF555OKO/@@download/ENCFF555OKO.bed.gz",
    "https://www.encodeproject.org/files/ENCFF671RKZ/@@download/ENCFF671RKZ.bed.gz",
    "https://www.encodeproject.org/files/ENCFF698RYW/@@download/ENCFF698RYW.bed.gz",
    "https://www.encodeproject.org/files/ENCFF021PYM/@@download/ENCFF021PYM.bed.gz",
    "https://www.encodeproject.org/files/ENCFF348DEB/@@download/ENCFF348DEB.bed.gz",
    "https://www.encodeproject.org/files/ENCFF694ORI/@@download/ENCFF694ORI.bed.gz",
    "https://www.encodeproject.org/files/ENCFF146OHQ/@@download/ENCFF146OHQ.bed.gz"
  ),
  "H3K4me1" = c(
    "https://www.encodeproject.org/files/ENCFF540NGG/@@download/ENCFF540NGG.bed.gz",
    "https://www.encodeproject.org/files/ENCFF750HRF/@@download/ENCFF750HRF.bed.gz",
    "https://www.encodeproject.org/files/ENCFF238YJA/@@download/ENCFF238YJA.bed.gz",
    "https://www.encodeproject.org/files/ENCFF861QCD/@@download/ENCFF861QCD.bed.gz",
    "https://www.encodeproject.org/files/ENCFF209XDL/@@download/ENCFF209XDL.bed.gz",
    "https://www.encodeproject.org/files/ENCFF699ABD/@@download/ENCFF699ABD.bed.gz",
    "https://www.encodeproject.org/files/ENCFF572NWJ/@@download/ENCFF572NWJ.bed.gz",
    "https://www.encodeproject.org/files/ENCFF035USI/@@download/ENCFF035USI.bed.gz",
    "https://www.encodeproject.org/files/ENCFF666PTL/@@download/ENCFF666PTL.bed.gz",
    "https://www.encodeproject.org/files/ENCFF321BVG/@@download/ENCFF321BVG.bed.gz",
    "https://www.encodeproject.org/files/ENCFF095KMQ/@@download/ENCFF095KMQ.bed.gz",
    "https://www.encodeproject.org/files/ENCFF148TUG/@@download/ENCFF148TUG.bed.gz"
  ),
  "H3K4me2" = c(
    "https://www.encodeproject.org/files/ENCFF749KLQ/@@download/ENCFF749KLQ.bed.gz",
    "https://www.encodeproject.org/files/ENCFF946FQI/@@download/ENCFF946FQI.bed.gz",
    "https://www.encodeproject.org/files/ENCFF438DUR/@@download/ENCFF438DUR.bed.gz",
    "https://www.encodeproject.org/files/ENCFF101KOJ/@@download/ENCFF101KOJ.bed.gz",
    "https://www.encodeproject.org/files/ENCFF270NMV/@@download/ENCFF270NMV.bed.gz",
    "https://www.encodeproject.org/files/ENCFF283LNH/@@download/ENCFF283LNH.bed.gz",
    "https://www.encodeproject.org/files/ENCFF836LZM/@@download/ENCFF836LZM.bed.gz",
    "https://www.encodeproject.org/files/ENCFF429OQI/@@download/ENCFF429OQI.bed.gz",
    "https://www.encodeproject.org/files/ENCFF262SBC/@@download/ENCFF262SBC.bed.gz",
    "https://www.encodeproject.org/files/ENCFF926IMD/@@download/ENCFF926IMD.bed.gz",
    "https://www.encodeproject.org/files/ENCFF188VRU/@@download/ENCFF188VRU.bed.gz",
    "https://www.encodeproject.org/files/ENCFF904TGH/@@download/ENCFF904TGH.bed.gz",
    "https://www.encodeproject.org/files/ENCFF714YFS/@@download/ENCFF714YFS.bed.gz"
  ),
  "H3K4me3" = c(
    "https://www.encodeproject.org/files/ENCFF954FXD/@@download/ENCFF954FXD.bed.gz",
    "https://www.encodeproject.org/files/ENCFF521CXH/@@download/ENCFF521CXH.bed.gz",
    "https://www.encodeproject.org/files/ENCFF378AAO/@@download/ENCFF378AAO.bed.gz",
    "https://www.encodeproject.org/files/ENCFF706WUF/@@download/ENCFF706WUF.bed.gz",
    "https://www.encodeproject.org/files/ENCFF549DKP/@@download/ENCFF549DKP.bed.gz",
    "https://www.encodeproject.org/files/ENCFF187LLD/@@download/ENCFF187LLD.bed.gz",
    "https://www.encodeproject.org/files/ENCFF304CMG/@@download/ENCFF304CMG.bed.gz",
    "https://www.encodeproject.org/files/ENCFF711MPL/@@download/ENCFF711MPL.bed.gz",
    "https://www.encodeproject.org/files/ENCFF268RXB/@@download/ENCFF268RXB.bed.gz",
    "https://www.encodeproject.org/files/ENCFF965NTW/@@download/ENCFF965NTW.bed.gz",
    "https://www.encodeproject.org/files/ENCFF904UJC/@@download/ENCFF904UJC.bed.gz",
    "https://www.encodeproject.org/files/ENCFF704LTU/@@download/ENCFF704LTU.bed.gz",
    "https://www.encodeproject.org/files/ENCFF242DXJ/@@download/ENCFF242DXJ.bed.gz"
  ),
  "H3K9me2" = c(
    "https://www.encodeproject.org/files/ENCFF367ELM/@@download/ENCFF367ELM.bed.gz",
    "https://www.encodeproject.org/files/ENCFF577ZNV/@@download/ENCFF577ZNV.bed.gz",
    "https://www.encodeproject.org/files/ENCFF632KWA/@@download/ENCFF632KWA.bed.gz",
    "https://www.encodeproject.org/files/ENCFF300WHS/@@download/ENCFF300WHS.bed.gz",
    "https://www.encodeproject.org/files/ENCFF447MUW/@@download/ENCFF447MUW.bed.gz",
    "https://www.encodeproject.org/files/ENCFF069OCN/@@download/ENCFF069OCN.bed.gz",
    "https://www.encodeproject.org/files/ENCFF050ALA/@@download/ENCFF050ALA.bed.gz",
    "https://www.encodeproject.org/files/ENCFF538BHD/@@download/ENCFF538BHD.bed.gz"
  ),
  "H3K27me3" = c(
    "https://www.encodeproject.org/files/ENCFF684ZZH/@@download/ENCFF684ZZH.bed.gz",
    "https://www.encodeproject.org/files/ENCFF584RYA/@@download/ENCFF584RYA.bed.gz",
    "https://www.encodeproject.org/files/ENCFF468ZHH/@@download/ENCFF468ZHH.bed.gz",
    "https://www.encodeproject.org/files/ENCFF801AHF/@@download/ENCFF801AHF.bed.gz",
    "https://www.encodeproject.org/files/ENCFF294LZM/@@download/ENCFF294LZM.bed.gz",
    "https://www.encodeproject.org/files/ENCFF959VSM/@@download/ENCFF959VSM.bed.gz",
    "https://www.encodeproject.org/files/ENCFF394AZZ/@@download/ENCFF394AZZ.bed.gz",
    "https://www.encodeproject.org/files/ENCFF428MSN/@@download/ENCFF428MSN.bed.gz",
    "https://www.encodeproject.org/files/ENCFF296RYM/@@download/ENCFF296RYM.bed.gz",
    "https://www.encodeproject.org/files/ENCFF403ODS/@@download/ENCFF403ODS.bed.gz",
    "https://www.encodeproject.org/files/ENCFF360IXV/@@download/ENCFF360IXV.bed.gz",
    "https://www.encodeproject.org/files/ENCFF841QVP/@@download/ENCFF841QVP.bed.gz"
  ),
  "H3K36me3" = c(
    "https://www.encodeproject.org/files/ENCFF739SBL/@@download/ENCFF739SBL.bed.gz",
    "https://www.encodeproject.org/files/ENCFF423AWX/@@download/ENCFF423AWX.bed.gz",
    "https://www.encodeproject.org/files/ENCFF457NNG/@@download/ENCFF457NNG.bed.gz",
    "https://www.encodeproject.org/files/ENCFF475QVQ/@@download/ENCFF475QVQ.bed.gz",
    "https://www.encodeproject.org/files/ENCFF195FSD/@@download/ENCFF195FSD.bed.gz",
    "https://www.encodeproject.org/files/ENCFF233BNK/@@download/ENCFF233BNK.bed.gz",
    "https://www.encodeproject.org/files/ENCFF479TQU/@@download/ENCFF479TQU.bed.gz",
    "https://www.encodeproject.org/files/ENCFF924TKB/@@download/ENCFF924TKB.bed.gz",
    "https://www.encodeproject.org/files/ENCFF243SKO/@@download/ENCFF243SKO.bed.gz",
    "https://www.encodeproject.org/files/ENCFF584FSY/@@download/ENCFF584FSY.bed.gz",
    "https://www.encodeproject.org/files/ENCFF257ZFX/@@download/ENCFF257ZFX.bed.gz",
    "https://www.encodeproject.org/files/ENCFF813VFV/@@download/ENCFF813VFV.bed.gz",
    "https://www.encodeproject.org/files/ENCFF432EMI/@@download/ENCFF432EMI.bed.gz",
    "https://www.encodeproject.org/files/ENCFF127LEC/@@download/ENCFF127LEC.bed.gz",
    "https://www.encodeproject.org/files/ENCFF655QDU/@@download/ENCFF655QDU.bed.gz",
    "https://www.encodeproject.org/files/ENCFF054TNJ/@@download/ENCFF054TNJ.bed.gz"
  ),
  "H3F3A" = c(
    "https://www.encodeproject.org/files/ENCFF516NJV/@@download/ENCFF516NJV.bed.gz",
    "https://www.encodeproject.org/files/ENCFF864SBW/@@download/ENCFF864SBW.bed.gz",
    "https://www.encodeproject.org/files/ENCFF276HPR/@@download/ENCFF276HPR.bed.gz",
    "https://www.encodeproject.org/files/ENCFF052DKK/@@download/ENCFF052DKK.bed.gz",
    "https://www.encodeproject.org/files/ENCFF376NEQ/@@download/ENCFF376NEQ.bed.gz",
    "https://www.encodeproject.org/files/ENCFF421FUZ/@@download/ENCFF421FUZ.bed.gz",
    "https://www.encodeproject.org/files/ENCFF595UXF/@@download/ENCFF595UXF.bed.gz",
    "https://www.encodeproject.org/files/ENCFF926ZHI/@@download/ENCFF926ZHI.bed.gz",
    "https://www.encodeproject.org/files/ENCFF103RHB/@@download/ENCFF103RHB.bed.gz",
    "https://www.encodeproject.org/files/ENCFF549AYZ/@@download/ENCFF549AYZ.bed.gz",
    "https://www.encodeproject.org/files/ENCFF509WEL/@@download/ENCFF509WEL.bed.gz"
  )
)

# Pull the data and get consensus peaks
pblapply(seq(histList), function(i) {
  mark <- names(histList)[i]
  fls <- histList[[i]]
  message(mark)
  
  downnames <- pbsapply(fls, function(fl) {
    downname <- paste0("tmp/", basename(fl))
    if (! file.exists(gsub(downname, pattern =  "(.+)\\.gz", replacement =  "\\1"))) {
      download.file(fl, destfile = downname, quiet = TRUE)
      system(paste0("gunzip ", downname))
    } else {
      message("Already downloaded!")
    }
    gsub(downname, pattern =  "(.+)\\.gz", replacement =  "\\1")
  }, cl = makeCluster(length(fls)))
  if (! file.exists(paste0("misc-data/annotations/", mark, "_optimal.bed"))) {
    system(paste0("chipr -i ", paste0(downnames, collapse = " "), " -m 2 -o misc-data/annotations/", mark), wait = FALSE)
  } else {
    message("Already run!")
  }
})



### Get SkewR annotations ###

# Get genes
EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86 %>%
  AnnotationDbi::select(keys =  AnnotationDbi::keys(.),
                        columns = c("GENEID", "SEQNAME", "GENESEQSTART",
                                    "GENESEQEND", "SEQSTRAND")) %>%
  mutate(
    strand = case_when(
      SEQSTRAND == 1 ~ "+",
      SEQSTRAND == -1 ~ "-",
      TRUE ~ "*"
    ),
    seqnames = paste0("chr", SEQNAME)
  ) %>%
  select(seqnames, start = GENESEQSTART, end = GENESEQEND, name = GENEID, strand) %>%
  GenomicRanges::makeGRangesFromDataFrame() %>%
  rtracklayer::export.bed(con = "~/.rseq_genomes/hg38/hg38.ensGene.bed")

# Get CpG islands
annotationLst$hg38$CpG_Islands %>%
  as.data.frame() %>%
  GenomicRanges::makeGRangesFromDataFrame() %>%
  rtracklayer::export.bed(con = "~/.rseq_genomes/hg38/hg38.cpg.bed")


### Save annotations ###
# annotationLst <- RLSeq::annotationLst


