# Wrangle to get the table
library(tidyverse)
GENE_EXP_SAMPLES <- "analyses/Prepare-RMapDB-Tables/gene_exp_samples.csv.xz"
SALMON_OUT <- "analyses/Prepare-RMapDB-Tables/GeneExpression/quant"
GENE_EXP_TABLE <- "analyses/Prepare-RMapDB-Tables/gene_expression.csv"
GENOME <- "hg38"

gene_exp_samps <- read_csv(GENE_EXP_SAMPLES) %>%
  mutate(quant_file = file.path(SALMON_OUT,
                                paste0(gene_exp_sample_id,
                                       "_",
                                       GENOME), 
                                "quant.sf"))

library(EnsDb.Hsapiens.v86)
library(tximport)

tx2gene <- AnnotationDbi::select(
  EnsDb.Hsapiens.v86, keys = keys(EnsDb.Hsapiens.v86), columns = c("TXNAME", "GENEID")
)
fls <- gene_exp_samps$quant_file
names(fls) <- gene_exp_samps$gene_exp_sample_id

txi <- tximport(files = fls, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)


tpm <- txi$abundance %>%
  apply(MARGIN = 1:2, function(x) log2(x + 1)) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  pivot_longer(cols = contains("SRX")) %>%
  dplyr::rename(exp_sample_id = name,
         log2tpm = value)


vsd <- txi$counts %>% 
  round() %>%
  DESeq2::vst() %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  pivot_longer(cols = contains("SRX")) %>%
  dplyr::rename(exp_sample_id = name,
                vst = value)


counts <- txi$counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  pivot_longer(cols = contains("SRX")) %>%
  dplyr::rename(exp_sample_id = name,
                counts = value)


gene_exp <- purrr::reduce(list(counts, tpm, vsd), inner_join, by = c("gene_id", "exp_sample_id"))
write_csv(gene_exp, GENE_EXP_TABLE)
system(paste0("xz -f ", GENE_EXP_TABLE))

