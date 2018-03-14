# per moa
rm(list = ls())

library(dplyr)
library(ggplot2)

source("moa_evaluations.R")

load("workspace.RData")

mean.res <- moa_recall(sm = cr.mean, metadata = metadata, n.cores = 3, N = 100)
median.mad.cov.res <- moa_recall(sm = cr.median.mad.cov, metadata = metadata, n.cores = 3, N = 100)
mix.res <- moa_recall(sm = cr.mix, metadata = metadata, n.cores = 3, N = 100)
median.mad.2.res <- moa_recall(sm = cr.median.mad.2, metadata = metadata, n.cores = 3, N = 100)
median.mad.res <- moa_recall(sm = cr.median.mad, metadata = metadata, n.cores = 3, N = 100)

mean.n <- mean.res %>% filter(!is.na(p.value)) %>% mutate(p.value = p.adjust(p.value)) %>% filter(p.value < 0.05) %>% NROW()
median.mad.cov.n <- median.mad.cov.res %>% filter(!is.na(p.value)) %>% mutate(p.value = p.adjust(p.value)) %>% filter(p.value < 0.05) %>% NROW()
mix.n <- mix.res %>% filter(!is.na(p.value)) %>% mutate(p.value = p.adjust(p.value)) %>% filter(p.value < 0.05) %>% NROW()
median.mad.2.n <- median.mad.2.res %>% filter(!is.na(p.value)) %>% mutate(p.value = p.adjust(p.value)) %>% filter(p.value < 0.05) %>% NROW()
median.mad.n <- median.mad.res %>% filter(!is.na(p.value)) %>% mutate(p.value = p.adjust(p.value)) %>% filter(p.value < 0.05) %>% NROW()

D <- data.frame(method = c("median", "median+mad (concatenated)", "median+median (SNF)", "median+mad (SNF)", "median+mad+cov. (SNF)"), no.moas = c(mean.n, median.mad.2.n, mix.n, median.mad.n, median.mad.cov.n))

D <- D %>% mutate(method = factor(method, levels = rev(c("median", "median+mad (concatenated)", "median+median (SNF)", "median+mad (SNF)", "median+mad+cov. (SNF)"))))

ggplot(D, aes(x = method, y = no.moas, fill = method, order = method)) + geom_bar(stat = "identity") + scale_x_discrete(breaks = NULL)

save.image("per_moa.RData")
