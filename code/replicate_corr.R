rm(list = ls())
library(dplyr)
library(stringr)
library(foreach)
source("rep.corr.func.R")

fls <- list.files("../output")

get.profiles <- function(profile.type) {
  profiles.nrm <- foreach (fl = fls, .combine = rbind) %do% {
    if (profile.type == "cov") {
      x <- readr::read_csv(paste0("../output/", fl))  
    } else {
      pl <- str_split(fl, "_")[[1]][1]
      x <- readr::read_csv(paste0("../input/", pl, "_normalized_variable_selected.csv"))  
    }
    x
  }
  
  profiles.nrm
}

Pf.mean <- get.profiles("mean")
Pf.cov <- get.profiles("cov")

feat.mean <- Pf.mean %>% colnames()
feat.mean <- feat.mean[which(!str_detect(feat.mean, "Metadata_"))]

feat.cov <- Pf.cov %>% colnames()
feat.cov <- feat.cov[which(!str_detect(feat.cov, "Metadata_"))]

thr.mean <- non.rep.cor(list(data = Pf.mean), grp.var = "Metadata_broad_sample", feat.mean, quant = 0.95)
thr.cov <- non.rep.cor(list(data = Pf.cov), grp.var = "Metadata_broad_sample", feat.cov, quant = 0.95)

rep.mean <- rep.cor(list(data = Pf.mean), grp.var = "Metadata_broad_sample", feat.mean)
rep.cov <- rep.cor(list(data = Pf.cov), grp.var = "Metadata_broad_sample", feat.cov)

rep.mean %>%
  ungroup() %>%
  filter(!is.na(Metadata_broad_sample)) %>%
  mutate(tot = n()) %>%
  filter(cr > thr.mean) %>%
  mutate(frac.strong = n()/tot) %>%
  select(frac.strong) %>%
  slice(1) %>%
  print

rep.cov %>%
  ungroup() %>%
  filter(!is.na(Metadata_broad_sample)) %>%
  mutate(tot = n()) %>%
  filter(cr > thr.cov) %>%
  mutate(frac.strong = n()/tot) %>%
  select(frac.strong) %>%
  slice(1) %>%
  print
