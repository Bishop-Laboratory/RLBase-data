#' Make the SampleQualChar db table for RMapDB
#' 
#' @param input_table CSV containing RMapDB data
#' 
#' @importFrom magrittr %>%
#' @import rlang
makeSampleQualChar <- function(input_table) {
  readr::read_csv(input_table) %>%
    dplyr::select(
      SRX,
      dplyr::contains(
        "_"
      ) & ! 
        contains(
          c("ip_type", "paired_end", "sample_name",
          "clean_name", "strand_specific", 
          "corr_method", "rlfs_method",
          "mode_group", "Shearing_method")
        )
    ) %>%
    tidyr::pivot_longer(cols = ! all_of("SRX")) %>%
    dplyr::select(
      id = SRX,
      char_type = name,
      value
    ) 
}


if (interactive()) {
  RMAPDB_CSV <- "analyses/rmap_full_11_25_with_study.csv"
  RMAP_SAMPLES_CSV <- "analyses/Prepare-RMapDB-Tables/sample_quality_characteristics.csv"
  
  readr::write_csv(makeSampleQualChar(RMAPDB_CSV), RMAP_SAMPLES_CSV)
  system(paste0("xz ", RMAP_SAMPLES_CSV))
  
}

