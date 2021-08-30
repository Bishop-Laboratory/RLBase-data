#' Makes the "GeneExpSamples" table for RMapDB
#'
if (sys.nframe() == 0L) {
  
  require(rlang)
  require(magrittr)
  library(dplyr)
  source("analyses/Prepare-RMapDB-Tables/GeneRLoopOverlap_GenFeatRLoopOverlap/gene_rl_overlap__gf_rl_overlap.R")
  source("../RSeq/src/scripts/utils.R")
  load("../RSeq/src/data/available_genomes.rda")
  
  # TODO: Needs to become part of the snake pipe
  RLOOP_ORIG_MANIFEST <- "analyses/R_loop_map_accessions_11132020.xlsx"
  RMAPSAMPSRAW <- "analyses/rmap_full_11_25_with_study.csv"
  RMAPSAMPS <- "analyses/Prepare-RMapDB-Tables/rmap_samples.csv.xz"
  PRI <- "analyses/Prepare-RMapDB-Tables/GeneExpSamples_and_exp2RMap/rlorig_pri.rda"
  GENE_EXP_SAMPLES <- "analyses/Prepare-RMapDB-Tables/gene_exp_samples.csv"
  GENE_EXP_2_RMAP_SAMPLES <- "analyses/Prepare-RMapDB-Tables/gene_exp_to_rmap_samples.csv"
  
  # Get the public-run Info
  rlorig <- readxl::read_excel(RLOOP_ORIG_MANIFEST) %>%
    dplyr::select(-Run, -SpikeIn, -AddInfo, -Issues, -Reference) %>%
    dplyr::distinct(GSM, .keep_all = TRUE) %>%
    dplyr::filter(Species == "HS") %>%
    dplyr::filter(! grepl(x = GSM, pattern = "_|\\."))
  
  rmapraw <- read_csv(RMAPSAMPSRAW) %>%
    dplyr::select(- tidyselect:::where(is.numeric)) %>%
    dplyr::select(- c(mode_group:rlfs_method))
  
  if (! file.exists(PRI)) {
    rlorig_pri <- get_public_run_info(rlorig$GSM)
    save(rlorig_pri, file = PRI)
  } else {
    load(PRI)
  }
  
  # Add the study ID
  pairs <- c("SRP256283" = "SRP256284", "SRP058310"="SRP058311",
             "SRP214116" = "SRP214118", "SRP214116" = "SRP214119",
             "SRP193695" = "GSE130242")
  
  # Big data wrangle to get the mapping of exp samples to rmap samples
  geneexp_rmap_raw <- full_join(rmapraw,
            dplyr::select(rlorig_pri, accessions_original,
                          condition, sra_experiment),
            by = c("GSM" = "accessions_original",
                   "study" = "condition")) %>%
    select(GSM, study, sra_experiment) %>%
    unique() %>%
    mutate(study = unlist(purrr::map(study, function(x) if(x %in% pairs) names(pairs)[pairs == x] else x))) %>%
    right_join(rlorig, by = "GSM") %>%
    select(-GSM, SRX=sra_experiment) %>%
    distinct(SRX, .keep_all = TRUE) %>%
    group_by(study) %>% 
    mutate(keeper = "RNASEQ" %in% Type) %>%
    ungroup() %>%
    filter(keeper) %>%
    mutate(cond = paste0(study, "_", Cell, "_", Other)) %>%
    select(SRX, Type, cond) %>%
    arrange(cond) %>%
    filter(Type %in% c("DRIP", "DRIPc",
                       "Input", "RNASEQ",
                       "sDRIP", "ssDRIP")) %>%
    group_by(cond) %>%
    mutate(keeper2 = "RNASEQ" %in% Type & length(Type) > 1) %>%
    filter(keeper2) %>%
    select(-contains("keeper")) %>% 
    mutate(Type = ifelse(Type == "RNASEQ", "gene_exp_sample_id", "rmap_sample_id")) %>%
    filter("rmap_sample_id" %in% Type & length(Type) > 1) %>%
    ungroup() %>%
    mutate(row = row_number()) %>%
    pivot_wider( id_cols = c(row, cond), names_from = Type, values_from = SRX) %>%
    select(-row, condition = cond)
  
  
  # Join RMap and GeneExp
  genexp_2_rmap <- full_join(
    select(geneexp_rmap_raw, condition, rmap_sample_id) %>% na.omit(),
    select(geneexp_rmap_raw, condition, gene_exp_sample_id) %>% na.omit()
  )
  
  # Write out the gene exp samples
  genexp_2_rmap %>%
    select(gene_exp_sample_id, condition) %>%
    unique() %>%
    write_csv(GENE_EXP_SAMPLES)
  
  # Write out the gene_exp_to_rmap_samples
  genexp_2_rmap %>%
    select(-condition) %>%
    unique() %>%
    write_csv(GENE_EXP_2_RMAP_SAMPLES)
  
  # Compress
  message("Compress")
  system(paste0("xz -f ", GENE_EXP_SAMPLES))  
  system(paste0("xz -f ", GENE_EXP_2_RMAP_SAMPLES))  
  
}
