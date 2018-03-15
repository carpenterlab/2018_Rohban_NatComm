rm(list = ls())
library(ggplot2)
library(dplyr)

load("../../CDRP_2nd_moment/code/per_moa.RData")
g <- g + theme_bw() + theme(axis.text = element_text(size=20), text = element_text(size=15)) 
ggsave("per_moa_cdrp.png", g)

load("../../BBBC022_2nd_moment/code/per_moa.RData")
g <- g + theme_bw() + theme(axis.text = element_text(size=20), text = element_text(size=15))
ggsave("per_moa_bbbc022.png", g)

load("../../TA_ORF_2nd_moment/code/per_moa.RData")
g <- g + theme_bw() + theme(axis.text = element_text(size=20), text = element_text(size=15))
ggsave("per_moa_taorf.png", g)