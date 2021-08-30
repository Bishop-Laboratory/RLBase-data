# Script showing generation of RMapDB
# This script assumes RSeqCLI is installed

# 0. Libraries
library(tidyverse)

# 1. rmapSamps needed to be wrangled for RSeqCLI
dir.create("data-raw/RMapDB_generate/", showWarnings = FALSE)
RSeqR::rmapSamps %>%
  select(experiment = id,
         control, 
         mode,
         genome) %>%
  filter(! mode %in% c("bisDRIP", "SMRF")) %>%
  write_csv("data-raw/RMapDB_generate/samples.csv")

# 2. Get the genome sizes prepped

# 2. Build config.json (with activated conda env if needed)
system("RSeqCLI build data-raw/RMapDB_generate/ data-raw/RMapDB_generate/samples.csv")

# 3. Check the JSON
configs <- jsonlite::read_json("data-raw/RMapDB_generate/config.json", simplifyVector = TRUE)
as_tibble(configs) %>%
  unique() %>%  View()

