suppressPackageStartupMessages(library(tidyverse))


# Get args
args <- commandArgs(trailingOnly = TRUE)

# Location of the catalog excel file
if (! interactive()) {
  CATALOG <- args[1]
  MANIFEST <- args[2]
} else {
  CATALOG <- "rlbase-data/rlbase_catalog.xlsx"
  MANIFEST <- "rlbase-data/rlbase_manifest.csv"
}
message(CATALOG)


# Make the output dir
dir.create("rlbase-data/", showWarnings = FALSE)

# Get the catalog, condition map, and modes
catalog <- readxl::read_excel(CATALOG, sheet = "catalog")
condition_map <- readxl::read_excel(CATALOG, sheet = "condition_map")
modes <- readxl::read_excel(CATALOG, sheet = "modes")

# Group samples
catalogGrouped <- catalog %>%
  mutate(group = case_when(
    mode %in% modes$mode ~ "rl",
    mode %in% c("RNA-Seq") ~ "exp",
    TRUE ~ "other"
  ))

# Get the manifest for RMap
toBuild <- catalogGrouped %>%
  left_join(select(modes, -PMID), by = "mode") %>%
  filter(group != "other",
         ! bisulfite_seq | mode == "RNA-Seq") %>%
  unique() 

# Write csv
toBuild %>%
  write_csv(MANIFEST)


# Echo out
message("\nDone -- file written: ", args[2], "\n")
Sys.sleep(1)
toBuild %>% group_by(group) %>% tally() %>% print()
