#' Make the RMap Samples db table for RMapDB
#' 
#' @param input_table CSV containing RMapDB data
#' 
#' @importFrom magrittr %>%
#' @import rlang
makeRMapSamples <- function(input_table) {
  readr::read_csv(input_table) %>%
    dplyr::mutate(
      is_input = FALSE,
      is_rnh_like = case_when(
        Condition %in% c(
          "ACTD", "IgG", "RNH", "WKKD"
        ) ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    dplyr::select(
      id = SRX,
      sample_name,
      genome,
      lab = Group,
      study_id = study,
      tissue = Cell,
      genotype = Genotype,
      treatment = Other,
      condition = Condition,
      control,
      mode,
      paired_end,
      is_input,
      is_rnh_like
    ) 
}


if (interactive()) {
  RMAPDB_CSV <- "analyses/rmap_full_11_25_with_study.csv"
  RMAP_SAMPLES_CSV <- "analyses/Prepare-RMapDB-Tables/rmap_samples.csv"
  
  readr::write_csv(makeRMapSamples(RMAPDB_CSV), RMAP_SAMPLES_CSV)
  system(paste0("xz ", RMAP_SAMPLES_CSV))
  
}

