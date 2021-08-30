#' Makes the "AliasToGene" table for RMapDB
#'
#' @importFrom magrittr %>%
#' @import rlang
makeAliasToGene <- function() {
  
  # Get ensembl data
  dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()
  sqlQuery <- 'SELECT * FROM alias, genes WHERE alias._id == genes._id;'
  sqlQuery2 <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
  
  ## Query
  dat1 <- DBI::dbGetQuery(dbCon, sqlQuery)[,c(-1)] %>%
    dplyr::rename(entrez_id=gene_id)
  dat2 <- DBI::dbGetQuery(dbCon, sqlQuery2)[,c(-1)] %>%
    dplyr::select(`_id`, description=gene_name)
  
  ## Wrangle
  aliasSymbol <- left_join(dat1, dat2,  by = "_id") %>%
    dplyr::select(-`_id`) %>%
    unique()
  
  # Get ensembl data
  Ens <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86 %>%
    AnnotationDbi::select(.,
                          keys = AnnotationDbi::keys(.,
                                                     keytype = "GENEID"),
                          columns = c("ENTREZID")) %>%
    dplyr::rename(
      gene_id=GENEID, entrez_id=ENTREZID
    ) %>%
    mutate(entrez_id = as.character(entrez_id))
  
  # Ensembl data with entrez data and finalize dataset
  full_join(
    x = Ens, 
    y = aliasSymbol,
    by = c("entrez_id")
  ) %>%
    dplyr::select(alias = alias_symbol,
                  gene = gene_id) %>%
    dplyr::filter(! is.na(gene) & ! stringr::str_detect(gene, pattern = "LRG")) %>%
    unique()
}

if (interactive()) {
  
  ALIAS_TO_GENE <- "analyses/Prepare-RMapDB-Tables/alias_to_gene.csv"
  
  require(rlang)
  require(magrittr)
  
  readr::write_csv(makeAliasToGene(), file = ALIAS_TO_GENE)
  system(paste0("xz -f ", ALIAS_TO_GENE))  # xzip
}



