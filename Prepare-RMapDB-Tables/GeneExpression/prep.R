source("../RSeq/src/scripts/utils.R")
load("../RSeq/src/data/available_genomes.rda")

library(tidyverse)

GENE_EXP_SAMPLES <- "analyses/Prepare-RMapDB-Tables/gene_exp_samples.csv.xz"
RUN_SHEET <- "analyses/Prepare-RMapDB-Tables/GeneExpression/run_sheet.csv"
RUN_SHEET_SMALL <- "analyses/Prepare-RMapDB-Tables/GeneExpression/run_sheet.small.csv"

# Get gene exp samps
gene_exp_samps <- read_csv(GENE_EXP_SAMPLES)

# Get public run info
pri_ges <- get_public_run_info(gene_exp_samps$gene_exp_sample_id)

# Prep for running gene exp
pri_ges %>%
  filter(sra_experiment %in% gene_exp_samps$gene_exp_sample_id) %>%
  select(sample=sra_experiment, experiment, genome, paired_end) %>%
  write_csv(RUN_SHEET)

pri_ges %>%
  filter(sra_experiment %in% gene_exp_samps$gene_exp_sample_id) %>%
  select(sample=sra_experiment, experiment, genome, paired_end) %>%
  slice(1:3) %>%
  write_csv(RUN_SHEET_SMALL)


