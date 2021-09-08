


# Taken from the finalizeRLRegions.R script -- for plotting the 
# overlay of an idividual sample on the qval/signal dist of RLRegions
select(isct, rlregion=rlregion.x, sample=.source, 
       is_rlfs=is_rlfs.x, confidence_level=confidence_level.x, 
       signalVal=signalVal.y, pVal=pVal.y, qVal=qVal.y,
       contains("med"), contains("avg")) %>%
  group_by(rlregion) %>%
  mutate(
    samples = paste0(sample, collapse = "\n"),
    nSamples = length(samples)
  ) %>%
  distinct(rlregion, .keep_all = TRUE) %>%
  (isol = sum(sample == "SRX3084736")) %>%
  mutate(isol = ifelse(isol > 0, "overlap", "non-overlapping"),
         isol = factor(isol, levels=c("overlap", "non-overlapping"))) %>%
  arrange(desc(isol)) %>%
  ggplot(aes(x = medSignalVal, y = medQVal, color = isol)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10()
# Same but for plotting the sample qVal vs RLRegions avg qVal
select(isct, rlregion=rlregion.x, sample=.source, 
       is_rlfs=is_rlfs.x, confidence_level=confidence_level.x, 
       signalVal=signalVal.y, pVal=pVal.y, qVal=qVal.y,
       contains("med"), contains("avg")) %>%
  filter(sample == "SRX3084736") %>%
  ggplot(aes(x = qVal, y = medQVal, color = is_rlfs)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10()




