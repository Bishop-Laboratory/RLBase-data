#' Build Genome Masks
#' Helper Function which builds GRanges of masked genomes for RegioneR
#' @return A named list of GRanges object with the masked chrom sizes. 
buildGenomeMasks <- function() {
  
  # Find genomes which have masked genomes available
  available_genomes <- BSgenome::available.genomes()
  masked_genomes <- available_genomes[grep(available_genomes, 
                                           pattern = ".+\\.masked")]
  
  # Download any masked genomes which are not already available
  to_install <- masked_genomes[! masked_genomes %in% installed.packages()]
  if (length(to_install) > 0) {
    BiocManager::install(to_install, ask = FALSE)
  }
  
  # Get the masked genomes
  genomes <- gsub(masked_genomes, replacement = "\\1",
                  pattern = "^BSgenome\\.[a-zA-Z]+\\.UCSC\\.([a-zA-Z0-9]+)\\.masked$")
  maskLst <- lapply(genomes, getMask2)
  names(maskLst) <- genomes
  
  return(maskLst)
  
}

#' Build Mask for Genome as GRanges
#' Adapted from `regioneR::getMask()` in the 
#' @return A GRanges object with the masked chrom sizes. 
getMask2 <- function(genome) {
  gn <- regioneR::characterToBSGenome(genome)
  chrs <- as.character(GenomicRanges::seqnames(GRanges(GenomeInfoDb::seqinfo(gn))))
  bsgenome <- gn
  chr.masks <- sapply(chrs, function(chr) {
    mm <- Biostrings::masks(bsgenome[[chr]])
    if(is.null(mm)) {
      return(NULL)
    } else {
      mm <- Biostrings::collapse(mm)[[1]]
      return(mm)
    }})
  
  chr.masks <- sapply(chrs, function(chr) {
    if(is.null(chr.masks[[chr]])) {
      return(NULL)
    } else {
      return(GenomicRanges::GRanges(seqnames=S4Vectors::Rle(rep(chr, length(chr.masks[[chr]]))), ranges=chr.masks[[chr]]))
    }
  })
  #Combine the mask for each chromosome into a single mask
  mask <- GenomicRanges::GRanges(seqinfo = seqinfo(bsgenome))
  for(chr in chrs) {
    if(!is.null(chr.masks[[chr]])) {
      suppressWarnings(mask <- c(mask, chr.masks[[chr]]))
    }
  }
  
  return(mask)
}

# Main execution of data building steps
main <- function() {
  # Get the list of genome masks if not already available
  if (! file.exists("data-raw/maskLst.rda")) {
    maskLst <- RSeqR:::buildGenomeMasks()
    save(maskLst, file = "misc/maskLst.rda", compress = "xz")
  } else {
    load("data-raw/maskLst.rda")
  }
  
  # Shrink size of mask genome list
  # 1. Only keep masks which have RLFS
  maskLst <- maskLst[sapply(names(maskLst), RSeqR:::checkRLFSAnno)]
  # 2. Remove large masks
  maskLst <- maskLst[sapply(maskLst, function(x) {
    object.size(x) < 500000
  })]
  names(maskLst) <- paste0(names(maskLst), ".masked")
  # 3. Remove ranges of size 1
  genomeMasks <- lapply(names(maskLst), function(genome) {
    newGr <- maskLst[[genome]]
    newGr[GenomicRanges::width(newGr) > 1]
  })
  names(genomeMasks) <- names(maskLst)
  
  # Save the data with xz compression
  usethis::use_data(genomeMasks, compress = "xz", overwrite = TRUE)
  
  
}




