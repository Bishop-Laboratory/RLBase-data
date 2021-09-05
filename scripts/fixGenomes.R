suppressPackageStartupMessages(library(tidyverse))

# Get args
args <- commandArgs(trailingOnly = TRUE)

# Location of the catalog excel file
CATALOG=args[1]
message(CATALOG)
CONFIG=args[2]
message(CONFIG)

# get the config
config <- read_tsv(CONFIG)

# Get the catalog, condition map, and modes
fixGenome <- readxl::read_excel(CATALOG, sheet = "fixgenome")

# Change config genome
fixed <- config %>%
  left_join(
    fixGenome, by = c("experiment_original" = "experiment")
  ) %>%
  mutate(
    genome = case_when(
      is.na(genome.y) ~ genome.x,
      TRUE ~ genome.y
    )
  ) %>%
  select(-genome.x, -genome.y) 

# Save to file
write_tsv(fixed, CONFIG)

# Echo out
message("\nDone -- file written: ", args[2], "\n")
Sys.sleep(1)
fixed %>% group_by(genome) %>% tally() %>% print()
