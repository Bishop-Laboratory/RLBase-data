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
aws.s3::put_object(file = "rlbase-data/rlregions/rlregions_S96.bw", object = "misc/rlregions_S96.bw", bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)
aws.s3::put_object(file = "rlbase-data/rlregions/rlregions_dRNH.bw", object = "misc/rlregions_dRNH.bw", bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)

# Add in the consensus signal
linesNow <- sapply(c("S96", "dRNH"), function(x) {
  write_lines(paste0("<h2>Description</h2>\n<p>Name: ", x, "</p>"), file = paste0(GEN_BROWSE_HTML_DIR, x, ".bw.html"))
  paste0("\ntrack ", x,
         "\ntype ", "bigWig",
         "\nshortLabel ", paste0("RLRegions (", x, ")"),
         "\nlongLabel ",paste0("RLRegions consensus signal (", x, ")"),
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
linesNow <- paste0("hub RLBase (RL Regions)",
                   "\nshortLabel R-loop database (RL Regions)",
                   "\nlongLabel A web database for R-loops and R-loop mapping experiments (RL Regions)",
                   "\nuseOneFile on\nemail millerh1@uthscsa.edu\n\ngenome hg38\n", linesNow)
writeLines(linesNow, con = file.path(GEN_BROWSE_HTML_DIR, paste0("oneFile__RLRegions.txt")))

# https://genome.ucsc.edu/s/millerh1%40livemail.uthscsa.edu/RLBase
