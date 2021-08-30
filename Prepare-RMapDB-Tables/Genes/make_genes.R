#' Makes the "Genes" table for RMapDB
#'
#' @importFrom magrittr %>%
#' @import rlang
makeGenes <- function() {
  EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86 %>%
    AnnotationDbi::select(
      keys = AnnotationDbi::keys(.), 
      columns = c("SYMBOL", "GENEBIOTYPE", 
                  "SEQNAME", "GENESEQSTART", 
                  "GENESEQEND", "SEQSTRAND")
    ) %>% 
    dplyr::mutate(
      location = paste0(
        "chr", SEQNAME, ":", GENESEQSTART, 
        "-", GENESEQEND, ":", ifelse(SEQSTRAND == 1, "+", "-")
      )
    ) %>%
    dplyr::select(
      id = GENEID,
      symbol = SYMBOL,
      biotype = GENEBIOTYPE,
      location
    )
} 

if (interactive()) {
  
  GENES <- "analyses/Prepare-RMapDB-Tables/genes.csv"
  
  require(rlang)
  require(magrittr)
  
  readr::write_csv(makeGenes(), file = GENES)
  system(paste0("xz ", GENES))  # xzip
  
}
