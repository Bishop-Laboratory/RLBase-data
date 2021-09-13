#' Script for compiling latest RLBase data for RLHub
dir.create("../RLBase-data/misc-data/rlhub/", showWarnings = FALSE)
library(tidyverse)
library(RLSeq)
library(pbapply)
pbo <- pboptions(type="txt") 

AWS_HTTPS_URL <- "https://rlbase-data.s3.amazonaws.com/"

## 1. Curated Genomic Annotations
dir.create("../RLBase-data/misc-data/rlhub/annotations", showWarnings = FALSE)

# Obtain the annotation paths from local
annot_pat <- "([a-zA-Z0-9]+)/([a-zA-Z0-9_ \\-]+)\\.csv\\.gz"
annots <- tibble(Key = list.files("../RLBase-data/misc-data/annotations/", recursive = TRUE)) %>%
  dplyr::filter(grepl(Key, pattern = annot_pat)) %>%
  dplyr::mutate(genome = gsub(Key, pattern = annot_pat, replacement = "\\1"),
                db = gsub(Key, pattern = annot_pat, replacement = "\\2")) 

# Download into R environment as a list
annotationLst <- lapply(
  unique(annots$genome),
  function(genome) {
    message(genome)
    annotGen <- dplyr::filter(annots, genome == {{ genome }})
    annotations <- pbapply::pblapply(
      seq(nrow(annotGen)), 
      function(i) {
        fl <- annotGen$Key[i]
        readr::read_csv(
          paste0("../RLBase-data/misc-data/annotations/", fl), 
          show_col_types = FALSE, 
          progress = FALSE
        )
      }
    )
    names(annotations) <- annotGen$db
    annotations
  }
) 
names(annotationLst) <- unique(annots$genome)
annotationLst2 <- lapply(names(annotationLst), function(genome) {
  annotationGen <- annotationLst[[genome]]
  # Flatten annotations into tbl list
  subpat <- "(.+)__(.+)__(.+)"
  message(" - Preparing annotations...")
  annots <- lapply(annotationGen, FUN = function(x) {
    if ("tbl" %in% class(x)) {
      x
    } else {
      dplyr::bind_rows(x)
    }
  }) %>% 
    dplyr::bind_rows() %>%
    dplyr::mutate(
      comb = gsub(.data$name, pattern = subpat, 
                  replacement = "\\1__\\2", perl=TRUE)
    ) %>% dplyr::group_by(.data$comb) %>%
    {setNames(dplyr::group_split(.), dplyr::group_keys(.)[[1]])}
})
names(annotationLst2) <- names(annotationLst)
annotations <- annotationLst2
annotations$hg38 <- annotationLst2$hg38[! names(annotationLst2$hg38) %in% c("DNaseHS__DNaseHS")]
save(annotations, file = "../RLBase-data/misc-data/rlhub/annotations/annotations.rda", compress="xz")

## 2. RL Regions

dir.create("misc-data/rlhub/rlregions/", showWarnings = FALSE)
opts <- c("All", "dRNH", "S96")
lapply(opts, function(opt) {
  # Get the rlregions table
  fl <- paste0("rlbase-data/misc/rlregions_", opt, "_table.tsv")
  outtbl <- paste0("misc-data/rlhub/rlregions/", opt, "_table.rda")
  rlregions_table <- read_tsv(fl, show_col_types=FALSE, progress=FALSE)
  expcor <- paste0("rlbase-data/misc/", opt, "_rlExpCorr.tsv")
  exp <- read_tsv(expcor)
  left_join(rlregions_table, exp) %>%
    dplyr::select(-chrom) %>%
    save(compress = "xz", file = outtbl)
  
  # Get the annotated rlregions table
  fl <- paste0("rlbase-data/misc/rlregions_", opt, "_annotations.csv")
  outtbl <- paste0("misc-data/rlhub/rlregions/", opt, "_annotations.rda")
  rlregions_table <- read_csv(fl, show_col_types = FALSE, progress = FALSE)
  save(rlregions_table, compress = "xz", file = outtbl)
  
})

## 3. RL Samples

dir.create("misc-data/rlhub/rlsamples", showWarnings = FALSE)
# Get the condition names for matching exp to rl
expToCond <- read_csv("misc-data/rl_to_exp.csv", show_col_types = FALSE) %>%
  right_join(read_tsv("rlbase-data/rlbase_samples.tsv", show_col_types = FALSE), 
             by = c("rl" = "experiment")) %>%
  mutate(toChoose = gsub(toChoose, pattern = "_", replacement =  "|"))
matchConds <- pbsapply(seq(nrow(expToCond)), function(i) {
  if (! is.na(expToCond$toChoose[i])) {
    cols <- paste0(colnames(expToCond)[grep(colnames(expToCond[i,]),
                                            pattern = expToCond$toChoose[i])])
    dplyr::select(expToCond[i,], all_of(cols[order(cols)])) %>% 
      pivot_longer(everything()) %>% 
      summarise(paste0(value, collapse = "_")) %>% 
      pull(1)
  } else {
    NA
  }
})
expToCond$matchCond <- matchConds

# Summarise to RL-sample level and build links
rlsamples <- expToCond %>%
  dplyr::select(-toChoose, rlsample = rl, expsamples=exp) %>%
  group_by(rlsample) %>%
  mutate(expsamples = paste0(expsamples, collapse=",")) %>%
  ungroup() %>%
  dplyr::filter(group == "rl") %>%
  distinct() %>%
  dplyr::rename(exp_matchCond = matchCond) %>%
  relocate(expsamples, .after = numPeaks) %>%
  mutate(
    coverage_s3 = paste0("coverage/", rlsample, "_", genome, ".bw"),
    peaks_s3 = paste0("peaks/", rlsample, "_", genome, ".broadPeak"),
    fastq_stats_s3 = paste0("fastq_stats/", rlsample, "_", genome, "__fastq_stats.json"),
    bam_stats_s3 = paste0("bam_stats/", rlsample, "_", genome, "__bam_stats.txt"),
    report_html_s3 = paste0("report/html/", rlsample, "_", genome, ".html"),
    report_rda_s3 = paste0("report/rda/", rlsample, "_", genome, ".rda"),
    rlfs_rda_s3 = paste0("rlfs_rda/", rlsample, "_", genome, ".rlfs.rda")
  ) 
save(rlsamples, file = "misc-data/rlhub/rlsamples/rlsamples.rda", compress = "xz")

## 4. Get the RLFS res

dir.create("misc-data/rlhub/rlfsres", showWarnings = FALSE)

# Magic
MANIFEST <- "rlbase-data/rlpipes-out/config.tsv"
MANIFEST_FINAL <- "rlbase-data/rlbase_samples.tsv"
RLFSRDA <- "rlbase-data/misc/rlfsRes.rda"
TODISCARD <- "misc-data/todiscard.rda"

# Load the rlfsRes
load(RLFSRDA)

# Get the predictions
rlfsPred <- pbapply::pblapply(rlfsRes, predictCondition)

# Clean large function
rlfsData <- lapply(rlfsRes, FUN = function(x) {
  x$perTestResults$`regioneR::numOverlaps`$evaluate.function <- NULL
  x
})

# Combine
rlfsres <- lapply(names(rlfsData), function(x) {
  list(
    "rlfsData" = rlfsData[[x]],
    "rlfsPred" = rlfsPred[[x]]
  )
})
names(rlfsres) <- names(rlfsData)

# Save
save(rlfsResults, file = "misc-data/rlhub/rlfsres/rlfsres.rda", compress = "xz")

## 5. Models

dir.create("misc-data/rlhub/models", showWarnings = FALSE)
file.copy("misc-data/model/fftModel.rda", to = "misc-data/rlhub/models/fftModel.rda")
file.copy("misc-data/model/prepFeatures.rda", to = "misc-data/rlhub/models/prepFeatures.rda")

## 6. GSG Correlation

dir.create("misc-data/rlhub/gsg_correlation/", showWarnings = FALSE)

file.copy("misc-data/gsSignalRLBase.rda", "misc-data/rlhub/gsg_correlation/gsSignalRLBase.rda")

# Get the corr
load("misc-data/rlhub/gsg_correlation/gsSignalRLBase.rda")
corrPearson <- gsSignalRLBase %>%
  column_to_rownames("location") %>%
  as.matrix() %>%
  cor()
corrSpearman <- gsSignalRLBase %>%
  column_to_rownames("location") %>%
  as.matrix() %>%
  cor(method = "spearman")
corrKendall <- gsSignalRLBase %>%
  column_to_rownames("location") %>%
  as.matrix() %>%
  cor(method = "kendall")





