rm(list = ls())
library(dplyr)
library(ggplot2)

load("workspace.RData")

source("moa_evaluations.R")

metadata.t <- metadata %>%
  group_by(Metadata_broad_sample) %>%
  slice(1) %>%
  ungroup()

cmpd_classification <- Vectorize(cmpd_classification, "k0")
k <- seq(from = 8, to = 15, by = 1)

d.mean <- cmpd_classification(cr.mean, metadata.t, k, not.same.batch = not.same.batch) 
d.mix <- cmpd_classification(cr.mix, metadata.t, k, not.same.batch = not.same.batch) 
d.median.mad <- cmpd_classification(cr.median.mad, metadata.t, k, not.same.batch = not.same.batch) 
d.median.mad.2 <- cmpd_classification(cr.median.mad.2, metadata.t, k, not.same.batch = not.same.batch) 
d.median.mad.cov <- cmpd_classification(cr.median.mad.cov, metadata.t, k, not.same.batch = not.same.batch) 

l.mean <- lapply(d.mean[5, ], function(x) sum(x)) 
l.mix <- lapply(d.mix[5, ], function(x) sum(x))
l.median.mad <- lapply(d.median.mad[5, ], function(x) sum(x)) 
l.median.mad.2 <- lapply(d.median.mad.2[5, ], function(x) sum(x)) 
l.median.mad.cov <- lapply(d.median.mad.cov[5, ], function(x) sum(x)) 

D <- data.frame(method = "median", k = k, tp = (unlist(l.mean)))
D <- rbind(D, 
           data.frame(method = "median+median (SNF)", k = k, tp = (unlist(l.mix))))
D <- rbind(D, 
           data.frame(method = "median+mad (SNF)", k = k, tp = (unlist(l.median.mad))))
D <- rbind(D, 
           data.frame(method = "median+mad (concatenated)", k = k, tp = (unlist(l.median.mad.2))))
D <- rbind(D, 
           data.frame(method = "median+mad+cov. (SNF)", k = k, tp = (unlist(l.median.mad.cov))))

lvls <- c("median+mad+cov. (SNF)", "median+mad (SNF)", "median+median (SNF)", "median+mad (concatenated)", "median")
D <- D %>% mutate(method = factor(method, levels = lvls))

g <- ggplot(D, aes(x = k, y = tp, color = method, order = as.character(method))) + 
  geom_point() + 
  geom_line() + 
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_continuous(breaks = k, minor_breaks = k) +
  ylab("No. of treatments") + 
  ggtitle("No. of treatments with a \n relevant match in their k-NNs") +
  theme_bw() +
  theme(axis.text = element_text(size=20), text = element_text(size=15)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.title=element_blank())
print(g) 
ggsave("enr_classification.png", g)
