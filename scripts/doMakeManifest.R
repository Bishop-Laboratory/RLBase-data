suppressPackageStartupMessages(library(tidyverse))

# Get args
args <- commandArgs(trailingOnly = TRUE)

# Location of the catalog excel file
CATALOG=args[1]

# Make the output dir
dir.create("rmap-data/", showWarnings = FALSE)

# Get the catalog, condition map, and modes
catalog <- readxl::read_excel(CATALOG, sheet = "catalog")
condition_map <- readxl::read_excel(CATALOG, sheet = "condition_map")
modes <- readxl::read_excel(CATALOG, sheet = "modes")

# Modes that will not be built currently due to limitations in the pipeline
unbuildable <- modes %>%
  filter(bisulfite_seq) %>%
  pull(mode)

# Group samples
catalogGrouped <- catalog %>%
  mutate(group = case_when(
    mode %in% modes$mode ~ "rmap",
    mode %in% c("RNA-Seq") ~ "exp",
    TRUE ~ "other"
  ))

# Get the manifest for RMap
catalogGrouped %>%
  filter(group == "rmap",
         ! mode %in%  unbuildable) %>%
  unique() %>%
  write_csv(args[2])

# Echo out
message("\nDone -- file written: ", args[2], "\n")
