#' builds the conda env through reticulate
#' @param envName name for new conda env
#' @param ... other arguments passed to reticulate::conda_create()
#' @return File path to created environment
#' @importFrom dplyr %>%
#' @importFrom rlang .data
buildCondaEnv <- function(envName, ...) {
  
  
  # Create env and install
  reticulate::conda_create(envname = envName, ...)
  
  # Get Env Path
  envPath <- reticulate::conda_list() %>%
    dplyr::filter(name == envName) %>%
    dplyr::mutate(python = gsub(python, pattern = "/bin/python", replacement = "")) %>%
    dplyr::pull(python)
  
  return(envPath)
}


#' Get list of effective genome sizes
#' @param genome the UCSC genome name to check
#' @param envPath the path to conda env containing "khmer"
#' @return A data frame containing the effective genome sizes.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
getEffectiveGenomeSizes <- function(genome, 
                                    envPath,
                                    lengths = c(36, 50, 75, 100, 125, 
                                                150, 200, 250, 300)) {
  
  fasta_file <- paste0("ftp://hgdownload.soe.ucsc.edu/goldenPath/",
                       genome, "/bigZips/", genome, ".fa.gz")
  outFile <- file.path(tempdir(), paste0(genome, ".fa.gz"))
  
  # Download genome fasta if not exists
  if (! file.exists(outFile)) {
    download.file(fasta_file, destfile = outFile)
  } 
  
  # Find eff genome size for all lengths
  sizes <- lapply(lengths, getEffGenSize,
                  faFile=outFile,
                  envPath=envPath)
  names(sizes) <- paste0("eff_genome_size_", lengths, "bp")
  sizedf <- data.frame(sizes) %>%
    dplyr::mutate(UCSC_orgID = genome)
  
  return(sizedf)
}


#' Get effective genome size
#' @param len The read length for which to obtain genome sizes
#' @param faFile The genome fasta file to analyze
#' @param envPath The path to the conda env with khmer installed
#' @param condaPath The path to the conda binary. Default = reticulate::conda_binary()
#' @return A data frame containing the effective genome sizes.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
getEffGenSize <- function(len,
                          faFile,
                          envPath) {
  
  # Get back to binary
  condaPath <- reticulate::conda_binary()
  
  # This line sources the conda shell script, activates the khmerenv, and runs unique-kmers.py
  cmd <- paste0(". ", gsub(condaPath, pattern = "/bin/conda", 
                           replacement = "/etc/profile.d/conda.sh"), 
                " && conda activate ", 
                envPath, " && unique-kmers.py -k ", len, " -R ",
                faFile, "_", len, ".txt ", faFile)
  # If already run, then not needed to run again
  if (! file.exists(paste0(faFile, "_", len, ".txt"))) {
    system(cmd)
  } else {
    print(paste0("Already ran length: ", len))
  }
  # Get the part of the unique-kmers.py output with the eff genome size
  lines <- readLines(paste0(faFile, "_", len, ".txt"))
  lines <- lines[grep(x = lines,"number of unique k-mers: ")]
  size <- gsub(lines, pattern = "number of unique k-mers: \t([0-9]+)$",
               replacement = "\\1")
  return(as.numeric(size))
}


#' Build Available Genomes
#' Helper Function which builds the UCSC available genomes dataset
#' @param test if TRUE, will only build info set for first 2 genomes in UCSC.
#' @param ... arguments passed to `getEffectiveGenomeSizes()`
#' @return A data frame of avaialable genomes and corresponding info. 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
buildAvailableGenomes <- function(test=FALSE, ...) {
  
  envName <- buildCondaEnv(packages = "khmer", 
                           envName = "khmerEnv",
                           channel = c("bioconda", "conda-forge"))
  
  # from http://rstudio-pubs-static.s3.amazonaws.com/562103_092b7264f392482e827440cf1521363c.html
  api_genome <- restfulr::RestUri("http://api.genome.ucsc.edu/list")
  response <- restfulr::read(api_genome$ucscGenomes)[['ucscGenomes']]
  
  # Wrangle the data from UCSC Genome API
  available_genomes <- do.call(rbind.data.frame, response) %>%
    tibble::rownames_to_column(var = "UCSC_orgID") %>%
    dplyr::group_by(UCSC_orgID) %>%
    dplyr::mutate(genes_available = checkGenes(UCSC_orgID),
                  homer_anno_available = checkHomer(UCSC_orgID),
                  year = as.numeric(gsub(x = description, 
                                         replacement = "\\1",
                                         pattern = ".+ ([0-9]+) \\(.+")))
  # If testing, only use small genomes
  smallgen <- c("eboVir3", "wuhCor1")
  if (test) {
    available_genomes <- dplyr::filter(available_genomes, .data$UCSC_orgID %in% smallgen)
  }
  
  # Get the effective genome sizes
  eff_gen_sizes <- available_genomes %>%
    dplyr::pull(UCSC_orgID) %>%
    lapply(function(genome) {
      getEffectiveGenomeSizes(genome = genome, envPath = envName, ...)
    }) %>%
    dplyr::bind_rows()
  
  # Bind final table and return
  return(
    dplyr::left_join(available_genomes, eff_gen_sizes, by = "UCSC_orgID")
  )
}

# Main execution of the data building steps
main <- function() {
  available_genomes <- buildAvailableGenomes()
  usethis::use_data(available_genomes, overwrite = TRUE)
}

