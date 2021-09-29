suppressPackageStartupMessages(library(tidyverse))
library(pbapply)
pbo <- pboptions(type="txt") 

# Helper function for extracting condition type from condition map
getlabel <- function(mode, condition, condition_map) {
  modeLgl <- map_lgl(condition_map$Mode, function(x) {
    grepl(mode, pattern = x)
  }) & condition == condition_map$Condition
  res <- condition_map[modeLgl,] %>%
    pivot_longer(cols = c("POS", "NEG", "NULL")) %>%
    filter(value) %>%
    pull(name)
  if (mode == "RNA-Seq") {
    res <- NA
  }
  return(res)
}

# Get args
args <- commandArgs(trailingOnly = TRUE)

# Location of the catalog excel file
if (! interactive()) {
  CATALOG=args[1]
  message(CATALOG)
  CONFIG=args[2]
  message(CONFIG)
} else {
  CATALOG = "rlbase-data/rlbase_catalog.xlsx"
  CONFIG = "rlbase-data/rlpipes-out/config.tsv"
}


# get the config
config <- read_tsv(CONFIG)

# Get the catalog, condition map, and modes
fixGenome <- readxl::read_excel(CATALOG, sheet = "fixgenome")
condition_map <- readxl::read_excel(CATALOG, sheet = "condition_map")

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
  dplyr::select(-genome.x, -genome.y)

# Get the condition type added in -- write the final manifest
condTbl <- fixed %>% 
  filter(group == "rl")
condTbl$conds <- pbsapply(
  seq(condTbl$experiment),
  function(i) {
    mode <- condTbl$mode[i]
    condition <- condTbl$condition[i]
    getlabel(
      mode, 
      condition, 
      condition_map
    )
  }
)
fixed <- condTbl %>%
  dplyr::select(experiment, label=conds) %>%
  right_join(fixed, by = "experiment") %>%
  mutate(condition = case_when(
    condition == "ActD" ~ "Input",
    TRUE ~ condition
  )) %>%
  unique()

# Remove duplicates
if ("label.x" %in% colnames(fixed)) {
  fixed <- dplyr::select(fixed, -label.x, -label.y) %>%
    unique()
}
if ("condType" %in% colnames(fixed)) {
  fixed <- dplyr::select(fixed, -condType) %>% unique()
}
fixed <- fixed %>% 
  mutate(
  label = ifelse(label == "NULL", "NEG", label)
)

# Save to file
write_tsv(fixed, CONFIG)

# Echo out
message("\nDone -- file written: ", args[2], "\n")
Sys.sleep(1)
fixed %>% group_by(genome) %>% tally() %>% print()
