## Script to set up the genome browser tracks ##
# Guide is here: https://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html#UseOneFile

# Libraries
library(tidyverse)
library(magrittr)

### MAGIC ###
GEN_BROWSE_HTML_DIR <- "misc-data/RLBase_TrackHub/"
baseURL <- RLSeq:::RLBASE_URL
#############

## Make output dir
dir.create(GEN_BROWSE_HTML_DIR, showWarnings = FALSE)

# Make HTML files
load("misc-data/rlhub/rlsamples/rlsamples.rda")
samples <- rlsamples %>%
  filter(genome == "hg38") %>%
  mutate(
    bigDataUrl = file.path(RLSeq:::RLBASE_URL, coverage_s3),
    type = "bigWig",
    visibility = "hide",
    longLabel = ifelse(genotype != "WT",
                       paste0(tissue, " - ", mode, " (", genotype, ", ", other, ") [", rlsample, "]"), 
                       paste0(tissue, " - ", mode, " (", other, ") [", rlsample, "]")),
    track = rlsample,
    shortLabel = paste0(tissue, " (", condition, ") [", mode, "]"),
    line = paste0("<h2>Description</h2>\n<p>Sample ID: ", rlsample,
                  "</p><p>Cell/tissue: ", tissue,
                  "</p><p>Genotype: ", genotype,
                  "</p><p>Other_id: ", other,
                  "</p><p>Condition: ", condition,
                  "</p><p>Label: ", label,
                  "</p><p>Prediction: ", prediction, 
                  "</p><p>PMID: ", PMID,
                  "</P>"),
    html = file.path(baseURL, "misc", "RLBase_TrackHub", paste0(rlsample, ".bw.html")),
    autoScale = "on",
    smoothingWindow = "8",
    maxHeightPixels = "100:40:8"
  ) %T>% {
    group_by(., rlsample) %>%
      {setNames(group_split(.), group_keys(.)[[1]])} %>%  # Split tibble into list by group with names
      lapply(
        function(x) {
          write_lines(x$line, file = paste0(GEN_BROWSE_HTML_DIR, x$rlsample, ".bw.html"))
        }
      )
  } 


# oneFile
modecols <- RLSeq:::auxdata$mode_cols
modergb <- as.data.frame(t(col2rgb(modecols$col)))
samples <- samples %>%
  mutate(modeNow = case_when(mode %in% modecols$mode ~ mode,
                             TRUE ~ "misc"))
samples <- modergb %>%
  mutate(modeNow = modecols$mode,
         color = paste0(red, ",", green, ",", blue)) %>%
  select(modeNow, color) %>%
  inner_join(samples, by = "modeNow")

# oneFile -- one per mode type...
group_by(samples, modeNow) %>%
  {setNames(group_split(.), group_keys(.)[[1]])} %>%  # Split tibble into list by group with names
  pbapply::pblapply(., function(x) {
    x$priority <- c(seq(x$rlsample))
    linesNow <- sapply(seq(x$rlsample), function(i) {
      paste0("\ntrack ", x$track[i],
             "\ntype ", x$type[i],
             "\nshortLabel ", x$shortLabel[i],
             "\nlongLabel ", x$longLabel[i],
             "\nvisibility ", x$visibility[i],
             "\nbigDataUrl ", x$bigDataUrl[i],
             "\nhtml ", x$html[i],
             "\ncolor ", x$color[i], 
             "\npriority ", x$priority[i],
             "\nautoScale ", x$autoScale[i],
             "\nsmoothingWindow ", x$smoothingWindow[i],
             "\nmaxHeightPixels ", x$maxHeightPixels[i], "\n")
    })
    linesNow <- paste0(linesNow, collapse = "")
    linesNow <- paste0("hub RLBase (", x$modeNow[1], ")",
                       "\nshortLabel R-loop database (mode: ", x$modeNow[1], ")",
                       "\nlongLabel A web database for R-loops and R-loop mapping experiments (mode: ", x$modeNow[1], ")",
                       "\nuseOneFile on\nemail millerh1@uthscsa.edu\n\ngenome hg38\n", linesNow)
    writeLines(linesNow, con = file.path(GEN_BROWSE_HTML_DIR, paste0("oneFile__", x$modeNow[1],".txt")))
} )

# Put the bigWigs (probably should move this part)
# aws.s3::put_object(file = "rlbase-data/rlregions/rlregions_S96.bw", object = "misc/rlregions_S96.bw", bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)
# aws.s3::put_object(file = "rlbase-data/rlregions/rlregions_dRNH.bw", object = "misc/rlregions_dRNH.bw", bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)

sapply(c("S96", "dRNH"), function(x) {
  download.file("https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes", destfile = "tmp/chrom.sizes")
  infile <- paste0("rlbase-data/rlregions/rlregions_", x, ".narrowPeak")
  read_tsv(infile, skip = 1, col_names = FALSE) %>%
    mutate(X4 = paste0(x, "__", row_number())) %>% 
    write_tsv(infile, col_names = FALSE)
  infile2 <- paste0("rlbase-data/rlregions/rlregions_", x, ".sort.narrowPeak")
  outfile <- paste0("rlbase-data/rlregions/rlregions_", x, ".bb")
  system(paste0("tail -n +2 ", infile, " | awk '$5=$5/10' | sort -k1,1 -k2,2n > ", infile2))
  system(paste0("~/miniconda3/condabin/conda run -n rlbaseData bedToBigBed -type=bed6+4 -as=tmp/bigNarrowPeak.as ", infile2, " tmp/chrom.sizes ", outfile))
})
aws.s3::put_object(file = "rlbase-data/rlregions/rlregions_S96.bb", object = "misc/rlregions_S96.bb", bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)
aws.s3::put_object(file = "rlbase-data/rlregions/rlregions_dRNH.bb", object = "misc/rlregions_dRNH.bb", bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)

# Add in the consensus signal
linesNow <- sapply(c("S96", "dRNH"), function(x) {
  write_lines(paste0("<h2>Description</h2>\n<p>Name: ", x, " (signal)</p>"), file = paste0(GEN_BROWSE_HTML_DIR, x, ".bw.html"))
  paste0("\ntrack ", paste0(x, "_signal"),
         "\ntype ", "bigWig",
         "\nshortLabel ", paste0("Consensus signal (", x, ")"),
         "\nlongLabel ",paste0("Consensus R-loop signal (", x, ")"),
         "\nvisibility ", "show",
         "\nbigDataUrl ", file.path(baseURL, "misc", paste0("rlregions_", x, ".bw")),
         "\nhtml ", file.path(baseURL, "misc", "RLBase_TrackHub", paste0(x, ".bw.html")),
         "\ncolor ", ifelse(x == "S96", "209,86,77", "83,157,166"), 
         "\npriority ", sample(1:2, 1),
         "\nautoScale ", "on",
         "\nsmoothingWindow ", "8",
         "\nmaxHeightPixels ", "100:40:8", "\n")
})
linesNow <- paste0(linesNow, collapse = "")
linesNow2 <- sapply(c("S96", "dRNH"), function(x) {
  write_lines(paste0("<h2>Description</h2>\n<p>Name: ", x, " (peaks)</p>"), file = paste0(GEN_BROWSE_HTML_DIR, x, ".bb.html"))
  paste0("\ntrack ", paste0(x, "_peaks"),
         "\ntype ", "bigNarrowPeak",
         "\nshortLabel ", paste0("Consensus peaks (", x, ")"),
         "\nlongLabel ",paste0("Consensus R-loop peaks (", x, ")"),
         "\nvisibility ", "show",
         "\nbigDataUrl ", file.path(baseURL, "misc", paste0("rlregions_", x, ".bb")),
         "\nhtml ", file.path(baseURL, "misc", "RLBase_TrackHub", paste0(x, ".bb.html")),
         "\ncolor ", ifelse(x == "S96", "209,86,77", "83,157,166"), 
         "\npriority ", sample(1:2, 1),
         "\nautoScale ", "on",
         "\nsmoothingWindow ", "8",
         "\nmaxHeightPixels ", "100:40:8", "\n")
})
linesNow <- paste0(linesNow, paste0(linesNow2, collapse = ""), collapse = "")
write_lines(paste0("<h2>Description</h2>\n<p>Name: RL Regions (RLBase V1)</p>"), file = paste0(GEN_BROWSE_HTML_DIR, "rlregions_table.bb.html"))
load("misc-data/rlhub/rlregions/rlregions_table.rda")
locpat <- "(.+):(.+)-(.+):(.+)"
crp <- RColorBrewer::brewer.pal(name = "Greens", n = 8)
x <- seq(1000)
n <- 8
rgbTbl <- split(x, sort(x%%n)) %>% 
  lapply(as_tibble) %>% 
  bind_rows(.id = "id") %>%
  mutate(id = as.numeric(id) + 1) %>%
  inner_join(tibble(
    "col" = crp,
    "id" = seq(crp)
  )) %>%
  mutate(itemRgb = map_chr(col, function(x) {
    paste0(col2rgb(x)[,1], collapse = ",")
  })) %>% select(score=value, itemRgb)


dd <- rlregion_table %>%
  mutate(
    name = gsub(rlregion, pattern = "All_", replacement = ""),
    chrom = gsub(location, pattern = locpat, replacement = "\\1"),
    start = as.numeric(gsub(location, pattern = locpat, replacement = "\\2")),
    end = as.numeric(gsub(location, pattern = locpat, replacement = "\\3")),
    strand = gsub(location, pattern = locpat, replacement = "\\4"),
    confidence_score=ifelse(confidence_score > 2.0, 2.0, confidence_score),
    score = floor(1000*((confidence_score-min(confidence_score))/(max(confidence_score)-min(confidence_score)))),
    thickStart = start,
    thickEnd = end
  ) %>%
  left_join(rgbTbl, by = "score") %>%
  filter(end > start & ! is.na(itemRgb)) %>%
  select(chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb) %>%
  write_tsv("misc-data/rlhub/rlregions/rlregions_table.unsorted.bed", col_names = FALSE)
infile <- "misc-data/rlhub/rlregions/rlregions_table.unsorted.bed"
infile2 <- "misc-data/rlhub/rlregions/rlregions_table.bed"
outfile <- "misc-data/rlhub/rlregions/rlregions_table.bb"
system(paste0("sort -k1,1 -k2,2n ", infile, " > ", infile2))
system(paste0("~/miniconda3/condabin/conda run -n rlbaseData bedToBigBed -type=bed9 ", infile2, " tmp/chrom.sizes ", outfile))
aws.s3::put_object(file = outfile, object = "misc/rlregions_table.bb", bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)
linesNow3 <- paste0("\ntrack RLRegions",
       "\ntype ", "bigBed 9",
       "\nshortLabel RL Regions",
       "\nlongLabel RL Regions (RLBase v1)",
       "\nvisibility show",
       "\nbigDataUrl ", file.path(baseURL, "misc", "rlregions_table.bb"),
       "\nhtml ", file.path(baseURL, "misc", "RLBase_TrackHub", "rlregions_table.bb.html"),
       "\nitemRgb on",
       "\npriority 1",
       "\nmaxHeightPixels ", "100:40:8", "\n")
linesNow <- paste0(linesNow, paste0(linesNow3, collapse = ""), collapse = "")

# Add the RLFS
write_lines(paste0("<h2>Description</h2>\n<p>Name: R-loop forming sequences (QmRLFS-finder.py)</p>"), file = paste0(GEN_BROWSE_HTML_DIR, "rlfs_hg38.bb.html"))
rlfs <- rtracklayer::import.bed(file.path(RLSeq:::RLBASE_URL, "rlfs-beds", "hg38.rlfs.bed"))
names(rlfs) <- paste0("RLFS__", seq(rlfs$name))
rlfs$name <- names(rlfs)
rlfs$score <- 1
si <- GenomeInfoDb::getChromInfoFromUCSC("hg38", as.Seqinfo = TRUE)
sl <- seqlevels(si)
seqlevels(rlfs) <- sl
GenomeInfoDb::seqinfo(rlfs) <- si
rtracklayer::export.bb(rlfs, con = "misc-data/rlhub/rlregions/rlfs_hg38.bb")
aws.s3::put_object(file = "misc-data/rlhub/rlregions/rlfs_hg38.bb", object = "misc/rlfs_hg38.bb", bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)
linesNow4 <- paste0("\ntrack RLFS",
                    "\ntype ", "bigBed",
                    "\nshortLabel RLFS",
                    "\nlongLabel R-loop forming sequences (QmRLFS-Finder.py)",
                    "\nvisibility show",
                    "\nbigDataUrl ", file.path(baseURL, "misc", "rlfs_hg38.bb"),
                    "\nhtml ", file.path(baseURL, "misc", "RLBase_TrackHub", "rlfs_hg38.bb.html"),
                    "\ncolor 115,50,168", 
                    "\npriority 1",
                    "\nautoScale ", "on",
                    "\nsmoothingWindow ", "8",
                    "\nmaxHeightPixels ", "100:40:8", "\n")
linesNow <- paste0(linesNow, paste0(linesNow4, collapse = ""), collapse = "")
linesNow <- paste0("hub RLBase (RL Regions)",
                   "\nshortLabel R-loop database (RL Regions)",
                   "\nlongLabel A web database for R-loops and R-loop mapping experiments (RL Regions)",
                   "\nuseOneFile on\nemail millerh1@uthscsa.edu\n\ngenome hg38\n", linesNow)
writeLines(linesNow, con = file.path(GEN_BROWSE_HTML_DIR, paste0("oneFile__RLRegions.txt")))




# https://genome.ucsc.edu/s/millerh1%40livemail.uthscsa.edu/RLBase
