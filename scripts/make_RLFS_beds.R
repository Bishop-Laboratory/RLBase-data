# Get RLFSs
get_rlfs <- function() {
  
    
  AVAILABLE_GENOMES <- ""
  QM_RLFS_FINDER <- "scripts/"
  
  script <- file.path(helpers_dir, "external", "QmRLFS-finder.py")
  outdir <- file.path(helpers_dir, "data", "rlfs")
  dir.create(outdir, showWarnings = FALSE)
  
  genomes <- available_genomes$UCSC_orgID[available_genomes$homer_anno_available]
  
  lapply(genomes, function(genome_now) {
    message(genome_now)
    if (! file.exists(file.path(outdir, paste0(genome_now, ".rlfs.out.table.txt"))) && 
        ! file.exists(file.path(outdir, paste0(genome_now, ".fa")))) {
      print("Getting genome!")
      download.file(paste0("ftp://hgdownload.soe.ucsc.edu/goldenPath/", genome_now, "/bigZips/", genome_now, ".fa.gz"),
                    destfile = file.path(outdir, paste0(genome_now, ".fa.gz")))
      R.utils::gunzip(file.path(outdir, paste0(genome_now, ".fa.gz")))
    }
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
  })
  
  # Clean up
  bed_files <- list.files(outdir, pattern = "\\.bed", full.names = F)
  lapply(bed_files, function(filenow) {
    old_file <- file.path(outdir, filenow)
    new_file <- file.path(paste0(outdir, '-beds'), filenow)
    file.copy(old_file, new_file, overwrite = TRUE)
  })
  
}
