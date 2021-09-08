# Taken from the finalizeRLRegions.R script -- for plotting the 
# overlay of an idividual sample on the qval/signal dist of RLRegions
rlregions %>%
  mutate(isol = grepl(samples, pattern = "SRX3084736", perl = TRUE),
         isol = ifelse(isol, "overlap", "non-overlapping"),
         isol = factor(isol, levels=c("overlap", "non-overlapping"))) %>%
  arrange(desc(isol)) %>%
  ggplot(aes(x = medSignalVal, y = medQVal, color = isol)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10()
# Same but for plotting the sample qVal vs RLRegions avg qVal
rlregionsOl %>%
  filter(sample == "SRX3084736") %>%
  ggplot(aes(x = qVal, y = medQVal, color = is_rlfs)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10()

