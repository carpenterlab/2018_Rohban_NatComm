rm(list = ls())

library(dplyr)
source("moa_evaluations.R")

cr.melt.mean <- readRDS("cr_mean.rds")
cr.melt.cov <- readRDS("cr_cov.rds")

cr.mean <- cr.melt.mean %>%
  select(Var1, Var2, value) %>%
  reshape2::acast("Var1 ~ Var2") 

k.snf <- 7

cr.cov <- cr.melt.cov %>%
  select(Var1, Var2, value) %>%
  reshape2::acast("Var1 ~ Var2")

af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = 0.5)
af.2 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = 0.5)
af.snf <- SNFtool::SNF(list(af.1, af.2), K = k.snf, t = 10)
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.mean.snf <- af.snf

af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = 0.5)
af.2 <- SNFtool::affinityMatrix(Diff = 1 - cr.cov, K = k.snf, sigma = 0.5)
af.snf <- SNFtool::SNF(list(af.1, af.2), K = k.snf, t = 10)
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.mix <- af.snf

metadata <- cr.melt.mean %>%
  select(Var1, Metadata_moa.x) %>%
  unique() %>%
  mutate(Metadata_broad_sample = Var1, Metadata_moa = Metadata_moa.x) %>%
  select(-Var1, -Metadata_moa.x)

#cr.mean <- cr.mean.snf

k <- 10

d.mean <- cmpd_classification(sm = cr.mean, metadata = metadata, k0 = k)
d.mix <- cmpd_classification(sm = cr.mix, metadata = metadata, k0 = k)

d.mean.sel <- d.mean %>%
  filter(pass)

d.mix.sel <- d.mix %>%
  filter(pass)

d.mean.sel %>% NROW() %>% print
d.mix.sel %>% NROW() %>% print

diff.cmpds <- setdiff(d.mix.sel$Var1, d.mean.sel$Var1)
d.mix %>%
  filter(Var1 %in% diff.cmpds) %>%
  select(-Metadata_moa.x, -pass) %>%
  left_join(d.mean, by = "Var1") %>%
  select(-pass) %>%
  rename(p.val.mix = p.val.x, p.val.mean = p.val.y, Compound = Var1) %>%
  arrange(p.val.mix) %>%
  knitr::kable() %>%
  print()

diff.cmpds <- setdiff(d.mean.sel$Var1, d.mix.sel$Var1)
d.mean %>%
  filter(Var1 %in% diff.cmpds) %>%
  select(-Metadata_moa.x, -pass) %>%
  left_join(d.mix, by = "Var1") %>%
  select(-pass) %>%
  rename(p.val.mix = p.val.y, p.val.mean = p.val.x, Compound = Var1) %>%
  arrange(p.val.mean) %>%
  knitr::kable() %>%
  print()
