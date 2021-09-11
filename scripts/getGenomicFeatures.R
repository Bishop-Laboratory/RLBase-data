library(tidyverse)
library(rtracklayer)
library(parallel)
library(pbapply)
library(GenomicFeatures)
pbo <- pboptions(type="txt") 

args <- commandArgs(trailing = TRUE)
CORES <- args[1]

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

# Get the genes in a list
message("- Make TxDBs")
if (! file.exists("tmp/genomes_txdb.rda")) {
  genomes <- list("hg38" = GenomicFeatures::makeTxDbFromGFF(file = "~/.rlpipes_genomes/hg38/hg38.ensGene.gtf"),
                  "mm10" = GenomicFeatures::makeTxDbFromGFF(file = "~/.rlpipes_genomes/mm10/mm10.ensGene.gtf"))
  saveDb(genomes$hg38, file = "tmp/hg38.db")
  saveDb(genomes$mm10, file = "tmp/mm10.db")
} else {
  hg38 <- loadDb(file = "tmp/hg38.db" )
  mm10 <- loadDb(file = "tmp/mm10.db")
  genomes <- list(
    "hg38" = hg38,
    "mm10" = mm10
  )
}

# For each genome, retrieve all annotations
message("- Building UCSC tables")
if (! file.exists("tmp/annotationLst_partI.rda")) {
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
  save(annotationLst, file = "tmp/annotationLst_partI.rda")
} else {
  load("tmp/annotationLst_partI.rda")
}


### Get DNase and TFBS ###
# DNase
dnase <- read_tsv("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/wgEncodeRegDnaseClustered.txt.gz",
                  col_names = c("x", "chrom", "start", "end"), show_col_types = FALSE, progress = FALSE) %>%
  mutate(type = "DNaseHS", db=type, strand="*",) %>%
  dplyr::select(chrom, start, end, strand, type, db)

# TFBS
tfbs <- read_tsv("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/encRegTfbsClustered.txt.gz",
                  col_names = c("x", "chrom", "start", "end", "type"), show_col_types = FALSE, progress = FALSE) %>%
  mutate(db="encodeTFBS", strand="*",) %>%
  dplyr::select(chrom, start, end, strand, type, db)

### Get G4Q DNA ###

message("- Getting G4Q DBs")

# Data was obtained from: https://doi.org/10.6084/m9.figshare.c.3498270.v1
# It contains predicted G4Q sites in bed files
dir.create("misc-data/G4Q", showWarnings = FALSE)
download.file(
  "https://figshare.com/ndownloader/files/6432597",
  destfile = "misc-data/G4Q/G4Q.tar.gz", quiet = TRUE
)
untar("misc-data/G4Q/G4Q.tar.gz", exdir = "misc-data/G4Q/")
file.rename("misc-data/G4Q/DATA/", "misc-data/G4Q/pred")

## Parse G4Q pred sites ##
# The cutoffs for labeling were determined adhoc from examintation of quantiles
# They attempt to bin the data such that interesting biology will be retained
# but dimensionality will be reduced.
# Refer to the paper for the meaning of the terms
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0165101
g4qpred <- read_tsv("misc-data/G4Q/pred/hg38_QFS_ALL.txt", col_names = c("chrom", "loc", "G4Pred", "seq", "strand"), show_col_types = FALSE) %>%
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
    type = paste0("G4Pred__", tractlen, "_", numloc, "_", g4num),
    db="G4Qpred"
  ) %>%
  dplyr::select(chrom, start, end, strand, type, db)

## Get the experimentally-obtained G4Q datasets
plus <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874_Na_K_PDS_plus_hits_intersect.bed.gz"
minus <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874_Na_K_PDS_minus_hits_intersect.bed.gz"
g4qexpp <- read_tsv(plus, col_names = c("chrom", "start", "end"), show_col_types = FALSE) %>% mutate(strand = "+", db="G4Qexp")
g4qexpm <- read_tsv(minus, col_names = c("chrom", "start", "end"), show_col_types = FALSE) %>% mutate(strand = "-", db="G4Qexp")
g4qexp <- bind_rows(g4qexpp, g4qexpm) %>% mutate(type = "G4Q_PDS_NaK_GSE63874")

### Get histone marks ###
message("- Getting histone marks")
# 2. Histone modifications
# Selected from the literature:
# https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC6259673/
# https://www.frontiersin.org/articles/10.3389/fgene.2020.00043/full
# https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC5123964/
histManifest <- read_csv("misc-data/histone_encode_manifest.csv", show_col_types = FALSE)
histList <- histManifest %>%
  dplyr::group_by(mark) %>%
  {setNames(group_split(.), group_keys(.)[[1]])}

# Pull the data and get consensus peaks
dir.create("misc-data/chipr_out/", showWarnings = FALSE)
pblapply(seq(histList), function(i, histList) {
  mark <- names(histList)[i]
  fls <- histList[[i]]$links
  message(mark)

  downnames <- pbapply::pbsapply(fls, function(fl) {
    downname <- paste0("tmp/", basename(fl))
    if (! file.exists(gsub(downname, pattern =  "(.+)\\.gz", replacement =  "\\1"))) {
      download.file(fl, destfile = downname, quiet = TRUE)
      system(paste0("gunzip ", downname))
    } else {
      message("Already downloaded!")
    }
    gsub(downname, pattern =  "(.+)\\.gz", replacement =  "\\1")
  }, cl = parallel::makeCluster(length(fls)))
  if (! file.exists(paste0("misc-data/chipr_out/", mark, "_optimal.bed"))) {
    system(paste0("chipr -i ", paste0(downnames, collapse = " "), " -m 2 -o misc-data/chipr_out/", mark), wait = TRUE)
  } else {
    message("Already run!")
  }
}, cl = parallel::makeCluster(length(histList)), histList=histList)

# Wrangle
histones <- pblapply(names(histList), function(mark) {
  read_tsv(paste0("misc-data/chipr_out/", mark, "_optimal.bed"), 
           show_col_types = FALSE, progress = FALSE,
           col_names = c("chrom", "start", "end", "name", "score",
                         "strand", "signalVal", "pVal", "qVal")) %>%
    mutate(
      type = !! mark,
      db = "Encode Histone"
    ) %>% dplyr::select(chrom, start, end, strand, type, db)
}) %>% bind_rows()

### Get SkewR annotations ###
message("- Running SkewR")
# hg38 #

# Genes
genesgr <- GenomicFeatures::genes(genomes$hg38)
rtracklayer::export.bed(genesgr, con = "~/.rlpipes_genomes/hg38/hg38.ensGene.bed")
toSymbhg38 <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86 %>%
  AnnotationDbi::select(keys =  AnnotationDbi::keys(.),
                        columns = "SYMBOL")
geneshg38 <- genesgr %>%
  as_tibble() %>%
  inner_join(toSymbhg38, by = c("gene_id" = "GENEID")) %>%
  pull(SYMBOL)

# CpG
annotationLst$hg38$CpG_Islands %>%
  GenomicRanges::makeGRangesFromDataFrame() %>%
  rtracklayer::export.bed(con = "~/.rlpipes_genomes/hg38/hg38.cpg.bed")

# mm10 #

# Genes
genesgr <- GenomicFeatures::genes(genomes$mm10)
rtracklayer::export.bed(genesgr, con = "~/.rlpipes_genomes/mm10/mm10.ensGene.bed")
toSymbmm10 <- EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79 %>%
  AnnotationDbi::select(keys =  AnnotationDbi::keys(.),
                        columns = "SYMBOL")
genesmm10 <- genesgr %>%
  as_tibble() %>%
  inner_join(toSymbmm10, by = c("gene_id" = "GENEID")) %>%
  pull(SYMBOL)

# CpG
annotationLst$mm10$CpG_Islands %>%
  GenomicRanges::makeGRangesFromDataFrame() %>%
  rtracklayer::export.bed(con = "~/.rlpipes_genomes/mm10/mm10.cpg.bed")

# Run SkewR 

sapply(names(genomes), function(genome) {
  message(genome)
  if (! file.exists(paste0("misc-data/SkewR/skewr_res_", genome, "/weak.txt"))) {
    system(paste0(
      "cd misc-data/SkewR/ && perl bin/RunGC-SKEW.pl -s ~/.rlpipes_genomes/",
      genome, "/", genome, ".fa -m model/GC_SKEW_20k.hmm -g ~/.rlpipes_genomes/", 
      genome, "/", genome, ".ensGene.bed -b ~/.rlpipes_genomes/",
      genome, "/", genome, ".cpg.bed -o skewr_res_", genome, " -z ", CORES))
  } else {
    message("Already run...")
  }
})

# Collate Skewr results
hg38_skewr <- read_tsv("misc-data/SkewR/skewr_res_hg38/probBedFile.bed", 
         col_names = c("chrom", "start", "end", "type", "score", "strand"), show_col_types = FALSE) %>%
  mutate(db = "skewr") %>%
  dplyr::select(chrom, start, end, strand, type, db)
mm10_skewr <- read_tsv("misc-data/SkewR/skewr_res_mm10/probBedFile.bed", show_col_types = FALSE,
                       col_names = c("chrom", "start", "end", "type", "score", "strand")) %>%
  mutate(db = "skewr") %>%
  dplyr::select(chrom, start, end, strand, type, db)

### Get RBPs ###
message("- Getting RBPs")
## Get RBP predictions from Ornament ##
# http://rnabiology.ircm.qc.ca/oRNAment/dashboard/
# HG38 only for now

dir.create("misc-data/hs_rbp_ornament", showWarnings = FALSE)
if (! file.exists("misc-data/hs_rbp_ornament/A1CF.bed.gz")) {
  download.file("http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/Homo_sapiens_oRNAment.bed.tar.gz", 
                destfile = "Homo_sapiens_oRNAment.bed.tar.gz")
  system("tar -xvzf misc-data/Homo_sapiens_oRNAment.bed.tar.gz -C misc-data/hs_rbp_ornament/ --strip-components=1")
}
fls <- list.files("misc-data/hs_rbp_ornament/", pattern = "\\.bed\\.gz$", full.names = TRUE)
rbps <- gsub(fls, pattern = ".+//([a-zA-Z0-9]+)\\.bed\\.gz", replacement = "\\1") %>% unique()
names(fls) <- rbps
tocheck <- rbps[which(! rbps %in% unique(tfbs$type))]
toProcess <- tocheck[! grepl(tocheck, pattern = "[a-z]+|^ENS|CONSTR|TAIR|CADAN|DRAFT|^TVAG|^ANIA|^RO3|^MAL|^NCU|STARP|^CG[0-9]+") &
        tocheck %in% geneshg38]
fls <- fls[toProcess]
RBPpred <- pblapply(seq(fls), function(i) {
  fl <- fls[i]
  rbp <- names(fl)
  read_tsv(fl, col_names = c("chrom", "start", "end", "type"), show_col_types = FALSE, progress = FALSE) %>%
    mutate(strand = "*", db = "RBPs_pred")
})
names(RBPpred) <- names(fls)

## Get eCLIP from encode (RNA-binding) ##

# Helper function for converting encode experiments to file URLs for peaks
encExp2Peaks <- function(acc, getIDR=FALSE) {
  resp <- httr::GET(url = paste0("https://www.encodeproject.org/search/?type=File&dataset=/experiments/",
                                 acc,"/&file_type=bed+narrowPeak&format=json&frame=object&limit=all"), httr::accept_json()) %>%
    httr::content(as="parsed")
  res <- lapply(seq(resp$`@graph`), function(i) {
    entry <- resp$`@graph`[[i]]
    tibble(
      accfile = entry$accession,
      sup = ! is.null(entry$superseded_by %>% unlist()),
      type = entry$file_type,
      assembly = entry$assembly,
      date = as.Date(entry$date_created),
      url = entry$s3_uri,
      output_type = entry$output_type,
      biol_reps = entry$biological_replicates %>% unlist() %>% paste0(collapse = "_"),
      acc = acc
    )
  }) %>% bind_rows() %>%
    mutate(best_samps = case_when(
      "1_2" %in% biol_reps ~ "1_2",
      TRUE ~ "1"
    )) %>%
    dplyr::filter(! sup,
                  assembly == "GRCh38",
                  biol_reps == best_samps) %>%
    mutate(dateMax = max(date)) %>%
    dplyr::filter(date == dateMax) %>%
    mutate(url = paste0("https://www.encodeproject.org/files/", accfile, "/@@download/", accfile, ".bed.gz"))
  if (nrow(res) > 1 & getIDR) {
    res <- dplyr::filter(res, output_type == "IDR thresholded peaks")
  } else {
    res
  }
}

# From https://www.nature.com/articles/s41586-020-2077-3
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2077-3/MediaObjects/41586_2020_2077_MOESM6_ESM.xlsx", 
              destfile = "misc-data/eclip_exps.xlsx")
eclips <- readxl::read_excel("misc-data/eclip_exps.xlsx", skip = 1) %>%
  dplyr::select(tissue=`Cell Line`,
                GENEID=`Official ENSGID`,
                encacc = `ENCODE Experiment Accession`) %>%
  mutate(GENEID = gsub(GENEID, pattern = "\\..+", replacement = "")) %>%
  inner_join(toSymbhg38)

forDownload <- pblapply(eclips$encacc, encExp2Peaks) %>% bind_rows()
manifestLst <- forDownload %>%
  inner_join(eclips, by = c("acc" = "encacc")) %>%
  write_csv("misc-data/eclips_rbps_manifest.csv") %>%
  group_by(SYMBOL) %>%
  {setNames(group_split(.), group_keys(.)[[1]])}
eclipres <- pblapply(seq(manifestLst), function(i) {
  manifest <- manifestLst[[i]]
  cols <- c("chrom", "start", "end", "name", "score", "strand", "pVal", "qVal")
  pklst <- lapply(seq(manifest$url), function(j) {
    pks <- read_tsv(manifest$url[j], show_col_types=FALSE, col_names=cols, progress = FALSE)
    pks %>%
      # Score is high (1000) vs low (200) confidence cutoff from peak caller
      dplyr::filter(score == 1000) 
  })
  names(pklst) <- paste0(manifest$SYMBOL, "_", manifest$tissue, "_", manifest$biol_reps)
  if (length(pklst) > 1) {
    ints <- valr::bed_intersect(pklst[[1]], pklst[2:length(pklst)])
    ints %>%
      mutate(start = ifelse(start.x > start.y, start.x, start.y),
             end = ifelse(end.x > end.y, end.x, end.y),
             strand = strand.x, 
             type = manifest$SYMBOL[1],
             db = "RBP_eCLiP") %>%
      dplyr::select(chrom, start, end, strand, type, db)
  } else {
    pklst %>%
      pluck(1) %>%
      mutate(
        type = manifest$SYMBOL[1],
        db = "RBP_eCLiP"
      ) %>%
      dplyr::select(chrom, start, end, strand, type, db)
  }
  
}) %>%
  bind_rows()

## Get ChIP from encode (Chromatin-binding) ##
# From https://www.nature.com/articles/s41586-020-2077-3
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2077-3/MediaObjects/41586_2020_2077_MOESM9_ESM.xlsx", 
              destfile = "misc-data/chiprbp_exps.xlsx")
chips <- readxl::read_excel("misc-data/chiprbp_exps.xlsx", skip = 0) %>%
  dplyr::select(tissue=Cell,
                SYMBOL=Target,
                encacc = `Accession Number`) 
forDownload <- pblapply(chips$encacc, encExp2Peaks, getIDR=TRUE) %>% bind_rows()
manifestLst <- forDownload %>%
  inner_join(chips, by = c("acc" = "encacc")) %>%
  write_csv("misc-data/chips_rbps_manifest.csv") %>%
  group_by(SYMBOL) %>%
  {setNames(group_split(.), group_keys(.)[[1]])}
chipres <- pblapply(seq(manifestLst), function(i) {
  manifest <- manifestLst[[i]]
  cols <- c("chrom", "start", "end", "name", "score", "strand", "pVal", "qVal")
  pklst <- lapply(seq(manifest$url), function(j) {
    read_tsv(manifest$url[j], show_col_types=FALSE, col_names=cols, progress = FALSE)
  })
  names(pklst) <- paste0(manifest$SYMBOL, "_", manifest$tissue, "_", manifest$biol_reps)
  if (length(pklst) > 1) {
    ints <- valr::bed_intersect(pklst[[1]], pklst[2:length(pklst)])
    ints %>%
      mutate(start = ifelse(start.x > start.y, start.x, start.y),
             end = ifelse(end.x > end.y, end.x, end.y),
             strand = strand.x, 
             type = manifest$SYMBOL[1],
             db = "RBP_ChIP") %>%
      dplyr::select(chrom, start, end, strand, type, db)
  } else {
    pklst %>%
      pluck(1) %>%
      mutate(
        type = manifest$SYMBOL[1],
        db = "RBP_ChIP"
      ) %>%
      dplyr::select(chrom, start, end, strand, type, db)
  }
}) %>%
  bind_rows()

### Cohesin ###
message("- Getting cohesin")
coh <- list("SA1"="STAG1", "SA2"="STAG2")
pblapply(seq(coh), function(i) {
  prot <- names(coh)[i]
  gen <- coh[[i]]
  fls <- list.files("misc-data/cohesin_peaks/", pattern = prot, full.names = TRUE)
  if (! file.exists(paste0("misc-data/chipr_out/", gen, "_optimal.bed"))) {
    system(paste0("chipr -i ", paste0(fls, collapse = " "), " -m 2 -o misc-data/chipr_out/", gen))
  }
})

# Wrangle
cohesins <- pblapply(coh, function(mark) {
  read_tsv(paste0("misc-data/chipr_out/", mark, "_optimal.bed"), 
           show_col_types = FALSE, progress = FALSE,
           col_names = c("chrom", "start", "end", "name", "score",
                         "strand", "signalVal", "pVal", "qVal")) %>%
    mutate(
      type = !! mark,
      db = "Cohesin"
    ) %>% dplyr::select(chrom, start, end, strand, type, db)
}) %>% bind_rows()

## Save annotations
message("- Saving")
# Add to annotations
annotationLst$hg38 <- c(as.list(annotationLst$hg38), list(
  "DNaseHS" = dnase,
  "encodeTFBS" = tfbs,
  "G4Qexp" = g4qexp,
  "G4Qpred" = g4qpred,
  "RBP_eCLiP" = eclipres,
  "RBP_ChIP" = chipres,
  "Encode Histone" = histones,
  "SkewR" = hg38_skewr,
  "Cohesin" = cohesins
))
annotationLst$mm10 <- c(as.list(annotationLst$mm10),
                        list(
                          "SkewR"=mm10_skewr
                        ))

# Write csv files
dir.create('misc-data/annotations/', showWarnings = FALSE)
a_ <- lapply(names(genomes), function(genome) {
  message(genome)
  annoOut <- paste0("misc-data/annotations/", genome)
  dir.create(annoOut, showWarnings = FALSE)
  annoLstGen <- annotationLst[[genome]]
  pbsapply(names(annoLstGen), function(x) {
    outname <- paste0(annoOut, "/", x, ".csv")
    write_csv(annoLstGen[[x]], outname, progress = FALSE)
  })
  return(NULL)
})

message("Done")
