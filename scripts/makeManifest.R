suppressPackageStartupMessages(library(tidyverse))


# Get args
args <- commandArgs(trailingOnly = TRUE)

# Location of the catalog excel file
if (! interactive()) {
  CATALOG <- args[1]
} else {
  CATALOG <- "rlbase-data/rlbase_catalog.xlsx"
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
  write_csv(args[2])

# Echo out
message("\nDone -- file written: ", args[2], "\n")
Sys.sleep(1)
toBuild %>% group_by(group) %>% tally() %>% print()
