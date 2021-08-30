# Script for creating rmapSamps

# This manifest was hand-curated and RSeq was used to generate peaks and
# coverage tracks for all the samples in it
require(rlang)
require(magrittr)
library(dplyr)

# TODO: Needs to become part of the snake pipe
RMAP_ORIG_MANIFEST <- "data-raw/R_loop_map_accessions_11132020.xlsx"
PRI <- "data-raw/rlorig_pri.rda"
RLSEQ_TYPES <- c("DRIP", "bisDRIP", "DRIPc", "DRIVE",
                 "MapR", "qDRIP", "R-ChIP", "RDIP",
                 "RNH-CnR", "RR-ChIP", "S1-DRIP", "sDRIP",
                 "SMRF", "ssDRIP")
NORMAL_LIKE <- c("D209N", "D210N", "delta-HC", "FLAG", "S96", "bisFoot")
INPUT_LIKE <- c("IgG", "Input")
RNH_LIKE <- c("RNH", "ACTD", "WKKD")
AVAILABLE_GENOMES <- "https://gist.github.com/millerh1/c3cf57c7715c7456c88f25c1526bb6fb/raw/c67abe8162030572351b332f670af45a103dbd75/available_genomes.tsv"
PUBLIC_RUN_INFO_GIST <- "https://gist.github.com/millerh1/1f4c32b9823567e446ccc504096243f5/raw/89375203bfd843b43ae6784124e5b227f5cd5e10/getPublicSampleInfo.R"


# Get data and source code
source(PUBLIC_RUN_INFO_GIST)
available_genomes <- readr::read_tsv(AVAILABLE_GENOMES)

# Wrangle tibble
rlorig <- readxl::read_excel(RMAP_ORIG_MANIFEST) %>%
  dplyr::select(-.data$Run, -.data$SpikeIn, -.data$AddInfo, -.data$Issues, -.data$Reference) %>%
  dplyr::distinct(.data$GSM, .keep_all = TRUE) %>%
  dplyr::filter(! grepl(x = .data$GSM, pattern = "_|\\."))

# Get the SRX for the input controls to match with normal-like samples
rlorigCont <- rlorig %>%
  filter(ControlSample != "NA")
controls <- rlorigCont %>%
  pull(ControlSample) %>%
  unique()

# Get public run info from SRA
if (! file.exists(PRI)) {
  rlorig_pri <- get_public_run_info(rlorig$GSM)
  rlorig_priInput <- get_public_run_info(controls)
  priFinal <- dplyr::bind_rows(rlorig_pri, rlorig_priInput)
  save(priFinal, file = PRI)
} else {
  load(PRI)
}

# Fix problematic samples that have the wrong genome annotation due to spikeins
HG38_SPIKES <- c("SRX7583981", "SRX7583982", "SRX7583983", "SRX7583984")
priFinal <- priFinal %>% 
  mutate(genome = case_when(sra_experiment %in% !! HG38_SPIKES ~ "hg38",
                            TRUE ~ genome))

# Add SRX for control sampels
rlcomb <- rlorig %>%
  left_join(select(priFinal, accessions_original, sra_experiment), 
            by = c("ControlSample"="accessions_original")) %>%
  mutate(ControlSample = sra_experiment) %>%
  select(-sra_experiment) %>%
  dplyr::left_join(y = priFinal, by = c("GSM" = "accessions_original")) %>%
  unique()

# Wrangle futher...
rmapSamps <- rlcomb %>%
  dplyr::filter(Type %in% !! RLSEQ_TYPES) %>%
  dplyr::rename(study=condition) %>%
  dplyr::mutate(condition = case_when(
    Condition == "RNASEH1" ~ "RNH",
    TRUE ~ Condition
  )) %>%
  dplyr::mutate(group = case_when(
    condition %in% !! NORMAL_LIKE ~ "Normal-like",
    condition %in% !! INPUT_LIKE ~ "Input-like",
    condition %in% !! RNH_LIKE ~ "RNH-like"
  )) %>%
  dplyr::select(id=sra_experiment, study, mode=Type,
                lab=Group, genome, tissue=Cell, genotype=Genotype,
                treatment=Other, condition, group, control=ControlSample,
                paired_end, read_length) %>%
  dplyr::mutate(lab=stringr::str_to_title(lab)) %>%
  unique()

# Save 
usethis::use_data(rmapSamps, overwrite = TRUE, compress = "xz")
