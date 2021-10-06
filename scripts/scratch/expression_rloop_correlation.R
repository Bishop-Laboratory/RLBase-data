#' Analyze the expression-r-loop correlation
rlsamples <- RLHub::rlbase_samples()
rlfsres <- RLHub::rlfs_res()
rlbaseRes <- RLHub::feat_enrich_samples()
gss <- RLHub::gs_signal()
rlregions <- RLHub::rlregions_meta()
geneexp <- RLHub::gene_exp()
vst <- geneexp@assays@data$vst
## Convert vst gene to rl ##
rlregion_to_gene <- rlregions %>%
  dplyr::select(rlregion, geneIDs) %>%
  dplyr::filter(! is.na(geneIDs)) %>%
  mutate(gene = map(geneIDs, function(x) {unlist(strsplit(x, split = ","))})) %>%
  dplyr::select(-geneIDs) %>%
  unnest(gene)
# TODO: Should this be a sum?
vstexp <- vst %>%
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -"gene") %>%
  left_join(rlregion_to_gene, by = "gene") %>%
  group_by(name, rlregion) %>%
  summarise(value = mean(value)) %>%
  dplyr::rename(exp = value)
condmap <- rlsamples %>%
  dplyr::select(rlsample, expsamples, exp_matchCond) %>%
  mutate(expsamples = map(expsamples, function(x) {unlist(strsplit(x, split = ","))})) %>%
  unnest(expsamples)
exp2cond <- unique(dplyr::select(condmap, -rlsample))
rl2cond <- unique(dplyr::select(condmap, -expsamples))
vstexp2 <- inner_join(exp2cond, vstexp,  by = c("expsamples" = "name")) %>%
  group_by(rlregion, exp_matchCond) %>%
  summarise(exp = mean(exp))
# Combine with R-loop signal
rlcounts <- RLHub::rlregions_counts()
vstrl <- rlcounts@assays@data$vst %>%
  as.data.frame() %>% 
  rownames_to_column(var = "rlregion") %>%
  pivot_longer(cols = -"rlregion") %>%
  dplyr::rename(rl = value)
# Match up to shared conds
vstrl2 <- inner_join(rl2cond, vstrl,  by = c("rlsample" = "name")) %>%
  group_by(rlregion, exp_matchCond) %>%
  summarise(rl = mean(rl))
# Combined expression and r-loop
vst_join <- inner_join(ungroup(vstrl2), ungroup(vstexp2), by = c("rlregion", "exp_matchCond"))
# Get corr
corr_estimate <- vst_join %>%
  group_by(rlregion) %>%
  group_split() %>%
  pbapply::pblapply(function(x) {
    ct <- suppressWarnings(cor.test(x$rl, x$exp, method="spearman"))
    tibble(
      rlregion = x$rlregion[1],
      pval = ct$`p.value`,
      estimate=ct$estimate
    )
  }) %>% bind_rows()
corr_estimate <- corr_estimate %>%
  mutate(corrPAdj = p.adjust(pval)) %>%
  dplyr::rename(corrPVal = pval,
                corrR = estimate)
corr_estimate_vst <- corr_estimate
rlregions <- rlregions %>% dplyr::select(
  -corrPAdj,
  -corrPVal,
  -corrR
)
rlregion_join <- left_join(rlregions, corr_estimate_vst, by = "rlregion")