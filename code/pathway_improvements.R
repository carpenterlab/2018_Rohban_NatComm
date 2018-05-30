rm(list = ls())

library(dplyr)
library(ggplot2)
library(stringr)
library(corrplot)
library(igraph)

load("workspace.RData")
set.seed(24)

thr <- 0.995

x <- sm.median.mad.cov %>% 
  filter(value > quantile(value, thr) & 
           same.moa) %>% 
  mutate(Metadata_moa.x = ifelse(str_length(Metadata_moa.x) < str_length(Metadata_moa.y), Metadata_moa.x, Metadata_moa.y)) %>%
  group_by(Metadata_moa.x) %>%
  tally() %>%
  arrange(-n) %>% 
  rename(Metadata_moa = Metadata_moa.x, n.median.mad.cov = n)

y <- sm.median.mad.2 %>% 
  filter(value > quantile(value, thr) & 
           same.moa) %>% 
  mutate(Metadata_moa.x = ifelse(str_length(Metadata_moa.x) < str_length(Metadata_moa.y), Metadata_moa.x, Metadata_moa.y)) %>%
  group_by(Metadata_moa.x) %>%
  tally() %>%
  arrange(-n) %>% 
  rename(Metadata_moa = Metadata_moa.x, n.median.mad.conc = n)

z <- sm.median.mad.2 %>% 
  filter(same.moa) %>% 
  mutate(Metadata_moa.x = ifelse(str_length(Metadata_moa.x) < str_length(Metadata_moa.y), Metadata_moa.x, Metadata_moa.y)) %>%
  group_by(Metadata_moa.x) %>%
  tally() %>%
  arrange(-n) %>% 
  rename(Metadata_moa = Metadata_moa.x, n.tot = n)

x %>%
  full_join(y, by = "Metadata_moa") %>%
  full_join(z, by = "Metadata_moa") %>%
  mutate(frac.median.mad.cov = n.median.mad.cov/n.tot, 
         frac.median.mad = n.median.mad.conc/n.tot) %>%
  select(-n.median.mad.cov, -n.median.mad.conc) %>%
  mutate(frac.median.mad.cov = ifelse(is.na(frac.median.mad.cov), 0, frac.median.mad.cov)) %>%
  mutate(frac.median.mad = ifelse(is.na(frac.median.mad), 0, frac.median.mad)) %>%
  arrange(-frac.median.mad.cov + frac.median.mad) %>%
  select(Metadata_moa, frac.median.mad.cov, frac.median.mad, n.tot) %>%
  mutate(frac.median.mad.cov = round(100 * frac.median.mad.cov)) %>%
  mutate(frac.median.mad = round(100 * frac.median.mad)) %>%
  filter(frac.median.mad != 0 | frac.median.mad.cov != 0) %>%
  htmlTable::htmlTable()
