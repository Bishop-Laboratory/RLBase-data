#' Get RLFS bed files
#' 
#' Uses QmRLFS-finder.py to calculate RLFS and construct bed files
#' Assumes your working directory is "RLBase-data"
get_rlfs <- function(cores=1, doRLFS=TRUE) {
  
  load("misc-data/available_genomes.rda")
  helpers_dir <- "scripts/"
  
  script <- file.path(helpers_dir, "QmRLFS-finder", "QmRLFS-finder.py")
  outdir <- file.path("misc-data", "rlfs")
  dir.create(outdir, showWarnings = FALSE)
  
  genomes <- available_genomes$UCSC_orgID[available_genomes$genes_available]
  
  parallel::mclapply(genomes, function(genome_now) {
    message(genome_now)
    if (! file.exists(file.path(outdir, paste0(genome_now, ".rlfs.out.table.txt"))) && 
        ! file.exists(file.path(outdir, paste0(genome_now, ".fa")))) {
      print("Getting genome!")
      if (! file.exists(file.path(outdir, paste0(genome_now, ".fa.gz")))) {
        download.file(paste0("http://hgdownload.soe.ucsc.edu/goldenPath/", genome_now, "/bigZips/", genome_now, ".fa.gz"),
                      destfile = file.path(outdir, paste0(genome_now, ".fa.gz")), method = "curl")
      }
      R.utils::gunzip(file.path(outdir, paste0(genome_now, ".fa.gz")))
    }
    if (doRLFS) {
      cmd <- paste0("python ", script, " -i ", file.path(outdir, paste0(genome_now, ".fa")),
                    " -o ", file.path(outdir, paste0(genome_now, ".rlfs")))
      if (! file.exists(file.path(outdir, paste0(genome_now, ".rlfs.out.table.txt")))) {
        print("Calculating RLFS")
        system(cmd)
      }
      
      # Convert to bed6 and liftOver
      cmd <- paste0("awk 'FNR > 1 {print $3\"\\t\"$21}' ", file.path(outdir, paste0(genome_now, ".rlfs.out.table.txt")),
                    "| awk '{gsub(\":|-\",\"\\t\", $1); print $1\"\\t\"\".\"\"\\t\"\".\"\"\\t\"$2}'",
                    " | bedtools sort -i stdin | mergeBed -i stdin -s -c 4,5,6 -o distinct > ", 
                    file.path(outdir, paste0(genome_now, ".rlfs.bed")))
      if (file.exists(file.path(outdir, paste0(genome_now, ".rlfs.out.table.txt"))) &&
          ! file.exists(file.path(outdir, paste0(genome_now, ".rlfs.bed")))) {
        print("Convert to bed!")
        system(cmd)
      }
    }
  }, mc.cores = cores)
  
  # Clean up
  if (doRLFS) {
    bed_files <- list.files(outdir, pattern = "\\.bed", full.names = F)
    dir.create(paste0(outdir, '-beds'), showWarnings = FALSE)
    sapply(bed_files, function(filenow) {
      old_file <- file.path(outdir, filenow)
      new_file <- file.path(paste0(outdir, '-beds'), filenow)
      file.copy(old_file, new_file, overwrite = TRUE)
      return(new_file)
    })
  }
  
}

# Run the function
args = commandArgs(trailingOnly=TRUE)
file_paths <- get_rlfs(cores=args[1], doRLFS=as.logical(args[2]))