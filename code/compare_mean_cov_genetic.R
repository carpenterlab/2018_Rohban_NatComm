rm(list = ls())

library(dplyr)
source("moa_evaluations.R")

cr.melt.mean <- readRDS("cr_mean.rds")
cr.melt.cov <- readRDS("cr_cov.rds")

cr.mean <- cr.melt.mean %>%
  select(Var1, Var2, value) %>%
  reshape2::acast("Var1 ~ Var2") 

k.snf <- round(NROW(cr.mean)^0.5 / 4)

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

ppi <- readRDS("../input/pr.pr.binary.rds")
ppi.extra <- readRDS("../input/pr.pr.interaction.rds")

metadata.2 <- readr::read_csv("../input/metadata_TA_2.csv")
metadata.1 <- readr::read_csv("../input/metadata_TA.csv")

metadata.2 <- metadata.2 %>%
  semi_join(metadata.1, by = "Metadata_broad_sample")

protein.interaction.comparison(sm = cr.mean, metadata = metadata.2, interaction.data = ppi, extra.interaction.data = ppi.extra, quant = 0.95, suffix = "mean")
protein.interaction.comparison(sm = cr.mix, metadata = metadata.2, interaction.data = ppi, extra.interaction.data = ppi.extra, quant = 0.95, suffix = "mix")

#print(d.mean)
#print(d.mix)