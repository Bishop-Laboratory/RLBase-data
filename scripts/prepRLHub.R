#' Script for compiling latest RLBase data for RLHub
dir.create("../RLBase-data/misc-data/rlhub/", showWarnings = FALSE)
library(tidyverse)
library(RLSeq)
library(pbapply)
pbo <- pboptions(type="txt") 

AWS_HTTPS_URL <- "https://rlbase-data.s3.amazonaws.com/"

## 1. Curated Genomic Annotations
message("- 1. Annotations")
dir.create("../RLBase-data/misc-data/rlhub/annotations", showWarnings = FALSE)

# Get the basic version of annotations
a_ <- pblapply(
  c("hg38", "mm10"), 
  function(genome) {
    
    # Wrangle file paths
    fl <- paste0("misc-data/annotations/annotations_", genome, ".rda")
    annots_prime <- paste0("misc-data/rlhub/annotations/annotations_primary_", genome, ".rda")
    annots_all <- paste0("misc-data/rlhub/annotations/annotations_full_", genome, ".rda")
    
    # Decide which databases to keep in primary annotations
    db_to_keep <- c("CpG_Islands", "Encode_CREs", "G4Qpred__G4Pred", 
                    "knownGene_RNAs", "PolyA", "Repeat_Masker", "skewr",
                    "snoRNA_miRNA_scaRNA", "Transcript_Features", "tRNAs")
    pat <- "(.+)__(.+)"
    load(fl)
    keep <- which(gsub(names(annotations), pattern = pat, replacement = "\\1") %in% db_to_keep)
    
    # Clean up some names
    names(annotations) <- gsub(names(annotations), pattern = "G4Qpred__G4Pred", replacement = "G4Qpred")
    
    # Load annotations
    message("loading")
    annotations <- lapply(
      annotations, 
      function(x) {
        mutate(x,
               id = as.numeric(gsub(
                 x$name, pattern = ".+__.+__([0-9]+)$", replacement = "\\1"
               ))) %>%
          dplyr::select(-name)
      }
    )
    
    # Subset annotations to get primary
    annotations_primary <- annotations[keep]
    annotations_all <- annotations
    
    # Save 
    message("Save primary")
    save(annotations_primary, file = annots_prime, compress = "xz")
    message("Save all")
    save(annotations_all, file = annots_all, compress = "xz")
    
  }
)
## 2. RL Regions
message("- 2. RL-Regions")
dir.create("misc-data/rlhub/rlregions/", showWarnings = FALSE)
opt <- "All"
# Get the rlregions table
fl <- paste0("rlbase-data/misc/rlregions_", opt, "_table.tsv")
# outtbl <- paste0("misc-data/rlhub/rlregions/", opt, "_table.rda")
outtbl <- paste0("misc-data/rlhub/rlregions/rlregions_table.rda")
rlregions_table <- read_tsv(fl, show_col_types=FALSE, progress=FALSE)
expcor <- paste0("rlbase-data/misc/", opt, "_rlExpCorr.tsv")
if (opt != "dRNH") {
  exp <- read_tsv(expcor, show_col_types = FALSE, progress = FALSE)
  rlregions_table <- left_join(rlregions_table, exp) %>%
    dplyr::select(-chrom) 
  save(rlregions_table, compress = "xz", file = outtbl)
} else {
  rlregions_table <- rlregions_table %>%
    mutate(
      corrR = NA,
      corrPVal = NA,
      corrPAdj = NA
    ) 
  save(rlregions_table, compress = "xz", file = outtbl)
  
}

# Get the annotated rlregions table
# outtbl <- paste0("misc-data/rlhub/rlregions/", opt, "_annotations.rda")
outtbl <- paste0("misc-data/rlhub/rlregions/rlregions_annotations.rda")
fl <- paste0("rlbase-data/misc/rlregions_", opt, "_annotations.csv")
rlregions_annotated <- read_csv(fl, show_col_types = FALSE, progress = FALSE)
rlregions_annotated$annotation <- gsub(rlregions_annotated$annotation, pattern = "G4Qpred__G4Pred", replacement = "G4Qpred")
save(rlregions_annotated, compress = "xz", file = outtbl)

## 3. RL Samples
message("- 3. RL-Samples")
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
    report_html_s3 = paste0("reports/", rlsample, "_", genome, ".html"),
    rlranges_rds_s3 = paste0("rlranges/", rlsample, "_", genome, ".rds"),
    rlfs_rda_s3 = paste0("rlfs_rda/", rlsample, "_", genome, ".rlfs.rda")
  ) 

# Clean up predictions (some ended up as NAs...)

rlres <- parallel::mclapply(seq(rlsamples$rlsample), function(i) {
  samp <- rlsamples$rlsample[i]
  gen <- rlsamples$genome[i]
  if (file.exists(file.path(
    "misc-data", "rlranges", paste0(samp, "_", gen, ".rds")
  ))) {
    rlr <- readRDS(file.path(
      "misc-data", "rlranges", paste0(samp, "_", gen, ".rds")
    ))
    rlr
  } else {
    NA
  }
}, mc.cores = 44)
names(rlres) <- rlsamples$rlsample

# Get prediction
rlsamples$prediction <- pbsapply(
  rlres, function(rlr) {
    if (! is.na(rlr)) {
      rlr@metadata$results@predictRes$prediction
    } else {
      NA
    }
  }
)

# Save
save(rlsamples, file = "misc-data/rlhub/rlsamples/rlsamples.rda", compress = "xz")

## 4. Get the RLFS res
message("- 4. RLFS-res")
dir.create("misc-data/rlhub/rlfsres", showWarnings = FALSE)

# Also fix the rlfs_rda
rlfsres <- pbsapply(seq(rlres), function(i) {
  rlr <- rlres[[i]]
  if (! is.na(rlr)) {
    lst <- list(
      rlfsPred=rlr@metadata$results@predictRes,
      rlfsData=rlr@metadata$results@rlfsRes
    )
    lst$rlfsData$perTestResults$`regioneR::numOverlaps`$evaluate.function <- NULL
    lst$rlfsData$perTestResults$`regioneR::numOverlaps`$randomize.function <- NULL
    lst
  } else {
    NA
  }
})
names(rlfsres) <- names(rlres)

# Save
save(rlfsres, file = "misc-data/rlhub/rlfsres/rlfsres.rda", compress = "xz")

## 5. Models
message("- 5. Models")
dir.create("misc-data/rlhub/models", showWarnings = FALSE)
file.copy("misc-data/model/fftModel.rda", overwrite = TRUE,
          to = "misc-data/rlhub/models/fftModel.rda")
file.copy("misc-data/model/prepFeatures.rda", overwrite = TRUE,
          to = "misc-data/rlhub/models/prepFeatures.rda")

## 6. GSG Correlation
message("- 6. GS Corr")
dir.create("misc-data/rlhub/gsg_correlation/", showWarnings = FALSE)

load("misc-data/gsSignalRLBase.rda")
save(gsSignalRLBase, file = "misc-data/rlhub/gsg_correlation/gsSignalRLBase.rda", compress = "xz")

## 7. Feature Enrichment Test Results
message("- 7. Feature enrichment")
dir.create("misc-data/rlhub/feature_enrichment/", showWarnings = FALSE)

# Feature enrichment of individual peaks
load("rlbase-data/misc/annotatedPeaks.rda")
feature_enrichment_per_sample <- resAnno
feature_enrichment_per_sample$db <- gsub(feature_enrichment_per_sample$db, pattern = "G4Qpred__G4Pred", replacement = "G4Qpred")
save(feature_enrichment_per_sample, file = "misc-data/rlhub/feature_enrichment/feature_enrichment_per_samples.rda", compress = "xz")

# Feature enrichment of RL-Regions
load("rlbase-data/misc/annotatedPeaks.rlregions.rda")
feature_enrichment_rlregions <- res
feature_enrichment_rlregions$db <- gsub(feature_enrichment_rlregions$db, pattern = "G4Qpred__G4Pred", replacement = "G4Qpred")
save(feature_enrichment_rlregions, file = "misc-data/rlhub/feature_enrichment/feature_enrichment_rlregions.rda", compress = "xz")

## 8. Gene Expression
message("- 8. Gene expression")
# This done within the buildExpression.R script...
file.copy("misc-data/expression/geneexp.rda",
          "misc-data/rlhub/expression/geneexp.rda", overwrite = TRUE)

## 9. RLRegion Counts (CTS, Norm, VST)
# This done within the rlregionCountMat.R script...


## 10. R-Loop binding proteins
dir.create("misc-data/rlhub/rlbps", showWarnings = FALSE)
# Update gene symbols and include Entrez IDs
rlbp <- readr::read_csv("misc-data/R_Loop_Binding_Proteins/Results/COMPARISONS/Combined.Proteins.csv")
symbs <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db),
  columns = c("ALIAS", "SYMBOL")
)
rlbpukn <- rlbp[! rlbp$geneName %in% symbs$SYMBOL,]
rlbpsym <- rlbp[rlbp$geneName %in% symbs$SYMBOL,]
rlbps <- inner_join(
  rlbpukn, symbs, by = c("geneName" = "ALIAS")
) %>% select(-geneName, -ENTREZID, geneName = SYMBOL) %>%
  bind_rows(rlbpsym) %>%
  distinct(geneName, .keep_all = TRUE) %>%
  relocate(geneName) %>%
  arrange(desc(combinedScore))

# Save
save(rlbps, file = "misc-data/rlhub/rlbps/rlbps.rda", compress = "xz")

## 11. RL and Expression correlation
#' Analyze the expression-r-loop correlation
load("misc-data/rlhub/rlsamples/rlsamples.rda")
load("misc-data/rlhub/rlregions/rlregions_table.rda")
load("misc-data/rlhub/expression/geneexp.rda")
tpm <- gene_exp@assays@data$tpm
rlregions <- rlregions_table
## Convert tpm gene to rl ##
rlregion_to_gene <- rlregions %>%
  dplyr::select(rlregion, geneIDs) %>%
  dplyr::filter(! is.na(geneIDs)) %>%
  mutate(gene = map(geneIDs, function(x) {unlist(strsplit(x, split = ","))})) %>%
  dplyr::select(-geneIDs) %>%
  unnest(gene)
# TODO: Should this be a sum?
tpmexp <- tpm %>%
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -"gene") %>%
  left_join(rlregion_to_gene, by = "gene") %>%
  group_by(name, rlregion) %>%
  summarise(value = mean(value)) %>%
  dplyr::rename(exp = value)
condmap <- rlsamples %>%
  dplyr::filter(prediction == "POS") %>%
  dplyr::select(rlsample, expsamples, exp_matchCond) %>%
  mutate(expsamples = map(expsamples, function(x) {unlist(strsplit(x, split = ","))})) %>%
  unnest(expsamples)
exp2cond <- unique(dplyr::select(condmap, -rlsample))
rl2cond <- unique(dplyr::select(condmap, -expsamples))
tpmexp2 <- inner_join(exp2cond, tpmexp,  by = c("expsamples" = "name")) %>%
  group_by(rlregion, exp_matchCond) %>%
  summarise(exp = mean(exp))
# Combine with R-loop signal
load("misc-data/rlhub/rlregions/rlregions_counts.rda")
tpmrl <- rlcounts@assays@data$tpm %>%
  as.data.frame() %>% 
  rownames_to_column(var = "rlregion") %>%
  pivot_longer(cols = -"rlregion") %>%
  dplyr::rename(rl = value)
# Match up to shared conds
tpmrl2 <- inner_join(rl2cond, tpmrl,  by = c("rlsample" = "name")) %>%
  group_by(rlregion, exp_matchCond) %>%
  summarise(rl = mean(rl))
# Combined expression and r-loop
tpm_join <- inner_join(ungroup(tpmrl2), ungroup(tpmexp2), by = c("rlregion", "exp_matchCond"))
# Get corr
corr_estimate <- tpm_join %>%
  group_by(rlregion) %>%
  group_split() %>%
  pbapply::pblapply(function(x) {
    ct <- suppressWarnings(cor.test(x$rl, x$exp, method="spearman"))
    tibble(
      rlregion = x$rlregion[1],
      pval = ct$`p.value`,
      estimate=ct$estimate
    )
  }) %>% bind_rows()
corr_estimate <- corr_estimate %>%
  mutate(corrPAdj = p.adjust(pval)) %>%
  dplyr::rename(corrPVal = pval,
                corrR = estimate)
rlregions <- rlregions %>% dplyr::select(
  -corrPAdj,
  -corrPVal,
  -corrR
)
rlregion_table <- left_join(rlregions, corr_estimate, by = "rlregion")
save(rlregion_table, file = "misc-data/rlhub/rlregions/rlregions_table.rda", compress = "xz")
tpm_rl_exp <- tpm_join
save(tpm_rl_exp, file = "misc-data/rlhub/rlregions/tpm_rl_exp.rda", compress = "xz")

# rlr <- rlregion_table %>%
#   select(corrR, allGenes) %>%
#   mutate(allGenes = map(allGenes, function(x) {unlist(strsplit(x, split = ","))})) %>%
#   unnest(allGenes) %>%
#   group_by(allGenes) %>%
#   summarize(
#     corrR = corrR[which.max(abs(corrR))]
#   ) %>%
#   write_csv("corGSEA.csv")
