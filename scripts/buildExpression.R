#' Script for wrangling the gene expression datasets 
#' and for matching the gene expression to r-loop mapping
# Wrangle to get the table

library(tidyverse)
library(pbapply)
pbo <- pboptions(type="txt") 

args <- commandArgs(trailingOnly = TRUE)

if (! interactive()) {
  SALMON_OUT <- args[1]
  GENE_EXP_TABLE <- args[2]
  rlbase_samples <- args[3]
} else {
  SALMON_OUT <- "rlbase-data/rlpipes-out/quant/"
  GENE_EXP_TABLE <- "rlbase-data/misc/gene_expression.csv"
  rlbase_samples <- "rlbase-data/rlbase_samples.tsv"
}

rlbase_samps <- read_tsv(rlbase_samples, show_col_types = FALSE)

## Match up expression and r-loops ##

# All allowed combinations to use for matching -- in order of priority
combs <- list(
  c("study", "tissue", "genotype", "other"),
  c("study", "tissue", "other"),
  c("study", "tissue", "genotype"),
  c("study", "tissue"),
  c("tissue", "genotype")
)

# Matching R-loop to expression
message("- Matching R-loop samples to expression samples via biological condition")
exp <- rlbase_samps %>% filter(group == "exp")
rl <- rlbase_samps %>% filter(group == "rl")
matches <- lapply(combs, function(comb) {
  inner_join(
    select(rl, experiment, all_of(combs[[1]])),
    select(exp, experiment, all_of(combs[[1]])),
    by = comb
  ) %>%
    select(rl=experiment.x, !!quo_name(paste0(comb, collapse = "_")) := experiment.y)
})
matchesJoined <- purrr::reduce(matches, full_join)
expChosen <- matchesJoined %>%
  mutate(
    toChoose = case_when(
      ! is.na(study_tissue_genotype_other) ~ "study_tissue_genotype_other",
      ! is.na(study_tissue_other) ~ "study_tissue_other",
      ! is.na(study_tissue_genotype) ~ "study_tissue_genotype",
      ! is.na(study_tissue) ~ "study_tissue",
      ! is.na(tissue_genotype) ~ "tissue_genotype",
      TRUE ~ "NA"
    )
  ) %>%
  pivot_longer(
    cols = contains("_")
  ) %>%
  filter(toChoose == name) %>%
  select(-name, exp=value) %>%
  distinct(rl, toChoose, exp)

# Save the matches
write_csv(expChosen, file = "misc-data/rl_to_exp.csv")

## Compile quant ##
message("- Compiling quant")
quants <- list.files(SALMON_OUT, recursive = TRUE,
                     pattern = "quant.sf", full.names = TRUE)
quanttbl <- tibble(
  sample = gsub(quants, pattern = ".+//([ES]{1}RX[0-9]+)_([a-zA-Z0-9]+)/quant.sf", replacement = "\\1"),
  genome = gsub(quants, pattern = ".+//([ES]{1}RX[0-9]+)_([a-zA-Z0-9]+)/quant.sf", replacement = "\\2"),
  quants
) %>%
  filter(genome == "hg38")
quanttbl <- left_join(quanttbl, exp, by = c("sample" = "experiment"))

tx2gene <- AnnotationDbi::select(
  EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, 
  keys = AnnotationDbi::keys(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86),
  columns = c("TXNAME", "GENEID")
)
fls <- quanttbl$quants
names(fls) <- quanttbl$sample

# Get TXI
txi <- tximport::tximport(files = fls, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Get variance-stabilizing transform
dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = quanttbl, design = ~1)
dds <- dds[,colSums(dds@assays@data$counts) > 5E6]
vsd <- DESeq2::vst(dds)
vsd <- vsd@assays@data[[1]]
vsd <- vsd %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  pivot_longer(cols = -gene_id) %>%
  dplyr::rename(exp_sample_id = name,
                vst = value)
# Get TPM
tpm <- txi$abundance %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  pivot_longer(cols = -gene_id) %>%
  dplyr::rename(exp_sample_id = name,
                log2tpm = value) %>%
  mutate(log2tpm = log2(log2tpm + 1))
# Get Counts
counts <- txi$counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  pivot_longer(cols = -gene_id) %>%
  dplyr::rename(exp_sample_id = name,
                counts = value)
# Compile
gene_exp <- purrr::reduce(list(counts, tpm, vsd), inner_join, by = c("gene_id", "exp_sample_id"))
write_csv(gene_exp, GENE_EXP_TABLE)

## Calculate gene expression correlation with R-loop levels ##
message("- Obtaining biological condition names")

# Get the condition names
expToCond <- expChosen %>%
  inner_join(rl, by = c("rl" = "experiment")) %>%
  mutate(toChoose = gsub(toChoose, pattern = "_", replacement =  "|"))
matchConds <- pbsapply(seq(nrow(expToCond)), function(i) {
  cols <- paste0(colnames(expToCond)[grep(colnames(expToCond[i,]),
                                          pattern = expToCond$toChoose[i])])
  select(expToCond[i,], all_of(cols[order(cols)])) %>% 
    pivot_longer(everything()) %>% 
    summarise(paste0(value, collapse = "_")) %>% 
    pull(1)
})
expToCond$matchCond <- matchConds

# Get the files
rlregionfiles <- list.files("rlbase-data/misc/", pattern = "rlregions.*\\.tsv", full.names = TRUE)
names(rlregionfiles) <- gsub(rlregionfiles, pattern = ".+//rlregions_([a-zA-Z0-9]+)_table.tsv$", replacement = "\\1")

# For each, get the condition -> expression mapping
message("- Matching expression to biological conditions")
rlCondExp <- pblapply(seq(rlregionfiles), function(i) {
  message(i)
  rlrfile <- rlregionfiles[i]
  cond <- names(rlrfile)
  rlrs <- read_tsv(rlrfile, show_col_types = FALSE)
  rlrs %>%
    filter(! is.na(geneIDs))  %>%
    mutate(
      samples = strsplit(samples, split = ","),
      geneIDs = strsplit(geneIDs, split = ",")
    ) %>%
    unnest(cols = geneIDs) %>%
    unnest(cols = samples) %>%
    inner_join(
      select(expToCond, rl, exp, matchCond),
      by = c("samples" = "rl")
    ) %>%
    select(rlregion, exp, matchCond, geneIDs) %>%
    inner_join(gene_exp, by = c("exp" = "exp_sample_id", "geneIDs" = "gene_id")) %>%
    group_by(rlregion, matchCond) %>%
    summarise(
      maxTPM = max(log2tpm),
      maxVST = max(vst)
    )
})
names(rlCondExp) <- names(rlregionfiles)

# For each, get the condition -> rlsignal mapping & combine
message("- Matching R-loops to expression")
rlsigfiles <- list.files("rlbase-data/rlregions/", pattern = "rlregions.*_signal\\.tsv", full.names = TRUE)
names(rlsigfiles) <- gsub(rlsigfiles, pattern = ".+//rlregions_([a-zA-Z0-9]+)_signal.tsv$", replacement = "\\1")
rlCondSignals <- pblapply(seq(rlsigfiles), function(i) {
  message(i)
  rlrfile <- rlsigfiles[i]
  cond <- names(rlrfile)
  rlrs <- read_tsv(rlrfile, show_col_types = FALSE)
  rlrs %>%
    select(rlregion, sample, signalVal, pVal, qVal) %>%
    inner_join(
      select(expToCond, sample=rl, matchCond)
    ) %>%
    group_by(rlregion, matchCond) %>%
    summarise(
      maxSignalVal = max(signalVal),
      maxPVal = max(pVal),
      maxQVal = max(qVal)
    ) %>%
    inner_join(
      rlCondExp[[cond]], by = c("rlregion", "matchCond")
    )
})
names(rlCondSignals) <- names(rlsigfiles)

# Get the correlations
rlCorrs <- pblapply(
  rlCondSignals[-2], function(x) {
    tokeep <- x %>% count() %>% filter(n > 2) %>% pull(rlregion)
    corres <- x %>%
      filter(rlregion %in% tokeep) %>%
      summarise(
        corrR = cor.test(x = maxQVal, y = maxVST, method = "spearman")$estimate,
        corrPVal = cor.test(x = maxQVal, y = maxVST, method = "spearman")$p.value
      ) 
    corres$corrPAdj <- p.adjust(corres$corrPVal)
    return(corres)
  }
)

# Save the correlations
lapply(names(rlCorrs), function(cond) {
  tblcor <- rlCorrs[[cond]]
  write_tsv(tblcor, file = paste0("rlbase-data/misc/", cond, "_rlExpCorr.tsv"))
})

message("Done")















