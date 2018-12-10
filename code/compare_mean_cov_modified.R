#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = T)

extends <- methods::extends

'compare_mean_cov_modified
Usage:
compare_mean_cov_modified -p <perturbation_type>
Options:
-h --help                                         Show this screen.
-p <perturbation_type> --pert=<perturbation_type> Either chemical or genetic' -> doc

opts <- docopt::docopt(doc)

p <- as.character(opts[["pert"]])

library(dplyr)
library(ggplot2)

source("moa_evaluations.R")

enrichment.based.classification <- FALSE
k.snf <- 7     # neighborhood size in SNF
t <- 10
k <- 1:10      # k top hits are used for classification
genetic <- (p == "genetic")
not.same.batch <- (!genetic)
snf.med.mad <- T

if (genetic) {
  not.same.batch <- F
}

cr.melt.mean <- readRDS("cr_median.rds")
cr.melt.cov <- readRDS("cr_cov.rds")
cr.melt.mad <- readRDS("cr_mad.rds")  
cr.melt.median.mad.2 <- readRDS("cr_median+mad.rds")
FA.cr.melt.mean <- readRDS("FA_cr_mean.rds")
FA.cr.melt.median <- readRDS("FA_cr_median.rds")
FA.cr.melt.mad <- readRDS("FA_cr_mad.rds")
pc.cr.melt.mean <- readRDS("pc_cr_mean.rds")
pc.cr.melt.median <- readRDS("pc_cr_median.rds")
pc.cr.melt.mad <- readRDS("pc_cr_mad.rds")
#dp.cr.melt.mean <- readRDS("dp_cr_mean.rds")

sigma.mean <- 0.5 
sigma.cov <- 0.5 
sigma.mad <- 0.5 

cr.mean <- cr.melt.mean %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2") 

cr.median.mad.2 <- cr.melt.median.mad.2 %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2") 

cr.mad <- cr.melt.mad %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2") 

cr.cov <- cr.melt.cov %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2")
FA.cr.mean <- FA.cr.melt.mean %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2")

FA.cr.median <- FA.cr.melt.median %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2")

FA.cr.mad <- FA.cr.melt.mad %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2")
pc.cr.mean <- pc.cr.melt.mean %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2") 
pc.cr.median <- pc.cr.melt.median %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2") 
pc.cr.mad <- pc.cr.melt.mad %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2")
# dp.cr.mean <- dp.cr.melt.mean %>%
#   select(Var1, Var2, value) %>%
#   group_by(Var1, Var2) %>%
#   summarise(value = max(value)) %>%
#   reshape2::acast("Var1 ~ Var2")

cr.mean <- sim_normalize(cr.mean) 
cr.cov <- sim_normalize(cr.cov) 
cr.median.mad.2 <- sim_normalize(cr.median.mad.2) 
cr.mad <- sim_normalize(cr.mad)
FA.cr.mean <- sim_normalize(FA.cr.mean)
FA.cr.median <- sim_normalize(FA.cr.median)
FA.cr.mad <- sim_normalize(FA.cr.mad)
pc.cr.mean <- sim_normalize(pc.cr.mean) 
pc.cr.median <- sim_normalize(pc.cr.median) 
pc.cr.mad <- sim_normalize(pc.cr.mad)
#dp.cr.mean <- sim_normalize(dp.cr.mean)

d <- apply(cr.mean, 1, function(x) !(sum(is.na(x)) >= (NROW(cr.mean) -1 )))
cr.mean <- cr.mean[d, d]

d <- apply(cr.cov, 1, function(x) !(sum(is.na(x)) >= (NROW(cr.cov) -1 )))
cr.cov <- cr.cov[d, d]

d <- apply(cr.mad, 1, function(x) !(sum(is.na(x)) >= (NROW(cr.mad) -1 )))
cr.mad <- cr.mad[d, d]

d <- apply(FA.cr.mean, 1, function(x) !(sum(is.na(x)) >= (NROW(FA.cr.mean) -1 )))
FA.cr.mean <- FA.cr.mean[d, d]

d <- apply(FA.cr.median, 1, function(x) !(sum(is.na(x)) >= (NROW(FA.cr.median) -1 )))
FA.cr.median <- FA.cr.median[d, d]

d <- apply(FA.cr.mad, 1, function(x) !(sum(is.na(x)) >= (NROW(FA.cr.mad) -1 )))
FA.cr.median <- FA.cr.mad[d, d]

d <- apply(pc.cr.mean, 1, function(x) !(sum(is.na(x)) >= (NROW(pc.cr.mean) -1 )))
pc.cr.mean <- pc.cr.mean[d, d]
d <- apply(pc.cr.median, 1, function(x) !(sum(is.na(x)) >= (NROW(pc.cr.median) -1 )))
pc.cr.median <- pc.cr.median[d, d]
d <- apply(pc.cr.mad, 1, function(x) !(sum(is.na(x)) >= (NROW(pc.cr.mad) -1 )))
pc.cr.mad <- pc.cr.mad[d, d]

#d <- apply(dp.cr.mean, 1, function(x) !(sum(is.na(x)) >= (NROW(dp.cr.mean) -1 )))
#dp.cr.mean <- dp.cr.mean[d, d]



#cm.rn <- setdiff(intersect(intersect(rownames(cr.mad), rownames(cr.mean)), intersect(rownames(cr.cov), rownames(cr.median.mad.2))), NA)
cm.rn <- intersect(intersect(rownames(cr.mad), rownames(cr.mean)), intersect(rownames(cr.cov), rownames(cr.median.mad.2))) 
FA.cm.rn <- intersect(intersect(rownames(FA.cr.mean), rownames(FA.cr.median)), intersect(rownames(FA.cr.mean), rownames(FA.cr.mad)))
pc.cm.rn <- intersect(intersect(rownames(pc.cr.mean), rownames(pc.cr.median)), intersect(rownames(pc.cr.mean), rownames(pc.cr.mad)))
int <- intersect(cm.rn, FA.cm.rn)
int1 <- intersect(int, pc.cm.rn)
#int2 <- intersect(int1, rownames(dp.cr.mean))
cm.rn <-setdiff(int1, NA)

cr.mean <- cr.mean[cm.rn, cm.rn]
cr.cov <- cr.cov[cm.rn, cm.rn]
cr.mad <- cr.mad[cm.rn, cm.rn]
cr.median.mad.2 <- cr.median.mad.2[cm.rn, cm.rn]
FA.cr.mean <- FA.cr.mean[cm.rn, cm.rn]
FA.cr.median <- FA.cr.median[cm.rn, cm.rn]
FA.cr.mad <- FA.cr.mad[cm.rn, cm.rn]
pc.cr.mean <- pc.cr.mean[cm.rn, cm.rn] 
pc.cr.median <- pc.cr.median[cm.rn, cm.rn]
pc.cr.mad <- pc.cr.mad[cm.rn, cm.rn]  
#dp.cr.mean <- dp.cr.mean[cm.rn, cm.rn] 

cr.mean[is.na(cr.mean)] <- 0
cr.median.mad.2[is.na(cr.median.mad.2)] <- 0
cr.mad[is.na(cr.mad)] <- 0
cr.cov[is.na(cr.cov)] <- 0
FA.cr.mean[is.na(FA.cr.mean)] <- 0
FA.cr.median[is.na(FA.cr.median)] <- 0
FA.cr.mad[is.na(FA.cr.mad)] <- 0
pc.cr.mean[is.na(pc.cr.mean)] <- 0
pc.cr.median[is.na(pc.cr.median)] <- 0
pc.cr.mad[is.na(pc.cr.mad)] <- 0
#dp.cr.mean[is.na(dp.cr.mean)] <- 0

af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = sigma.mean)
af.2 <- SNFtool::affinityMatrix(Diff = 1 - cr.mad, K = k.snf, sigma = sigma.mad)
af.snf <- SNFtool::SNF(list(af.1, af.2), K = k.snf, t = t)
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.median.mad <- af.snf

af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = sigma.mean)
af.2 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = sigma.mean)
af.snf <- SNFtool::SNF(list(af.1, af.2), K = k.snf, t = t)
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.mix <- af.snf

af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = sigma.mean)
af.2 <- SNFtool::affinityMatrix(Diff = 1 - cr.mad, K = k.snf, sigma = sigma.mad)
af.3 <- SNFtool::affinityMatrix(Diff = 1 - cr.cov, K = k.snf, sigma = sigma.cov)
af.snf <- SNFtool::SNF(list(af.1, af.2, af.3), K = k.snf, t = round(3/2 * t))
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.median.mad.cov <- af.snf

metadata <- cr.melt.mean %>%
  select(Var1, Metadata_moa.x, Metadata_Plate_Map_Name.x) %>%
  unique() %>%
  mutate(Metadata_broad_sample = Var1, Metadata_moa = Metadata_moa.x, Metadata_Plate_Map_Name = Metadata_Plate_Map_Name.x) %>%
  select(-Var1, -Metadata_moa.x, -Metadata_Plate_Map_Name.x)

metadata <- metadata %>% mutate(Metadata_moa = str_to_lower(Metadata_moa))

#cmpd_classification <- Vectorize(cmpd_classification, "k0")
cmpd_knn_classification <- Vectorize(cmpd_knn_classification, "k0")

if (enrichment.based.classification) {
  d.mean <- cmpd_classification(sm = cr.mean, metadata = metadata, k0 = k)
  d.mix <- cmpd_classification(sm = cr.mix, metadata = metadata, k0 = k)
  
  if (length(k) == 1) {
    d.mean <- (d.mean %>% t) %>% apply(., 2, function(x) {y <- data.frame(as.vector(x)); names(y) <- names(x); y}) 
    d.mix <- (d.mix %>% t) %>% apply(., 2, function(x) {y <- data.frame(as.vector(x)); names(y) <- names(x); y}) 
    
    d.mean <- d.mean %>% as.data.frame()
    d.mix <- d.mix %>% as.data.frame()
    
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
  } else {
    l.mean <- lapply(d.mean["pass", ], function(x) sum(x)) 
    l.mix <- lapply(d.mix["pass", ], function(x) sum(x))
    
    D <- data.frame(method = "mean", k = k, tp = (unlist(l.mean)))
    D <- rbind(D, 
               data.frame(method = "mean+cov.", k = k, tp = (unlist(l.mix))))
    g <- ggplot(D, aes(x = k, y = tp, color = method)) + 
      geom_point() + 
      geom_line() + 
      scale_y_continuous(limits = c(0, NA)) +
      scale_x_continuous(breaks = k, minor_breaks = k) +
      ylab("No. of treatment with a \n same MOA/Pathway treatment in their k-NNs")
    print(g)    
    ggsave("classification_comparison.png", g, width = 7, height = 5)
  }
} else {
  d.mean <- cmpd_knn_classification(cr.mean, metadata, k, not.same.batch = not.same.batch) 
  d.mix <- cmpd_knn_classification(cr.mix, metadata, k, not.same.batch = not.same.batch) 
  d.median.mad <- cmpd_knn_classification(cr.median.mad, metadata, k, not.same.batch = not.same.batch) 
  d.median.mad.2 <- cmpd_knn_classification(cr.median.mad.2, metadata, k, not.same.batch = not.same.batch) 
  d.median.mad.cov <- cmpd_knn_classification(cr.median.mad.cov, metadata, k, not.same.batch = not.same.batch) 
  
  if (length(k) == 1) {
    d.mean %>% NROW %>% print
    d.mix %>% NROW %>% print
    d.diff <- d.mean %>%
      full_join(d.mix, by = "Var1") %>%
      mutate(Metadata_moa = ifelse(is.na(Metadata_moa.x.x), Metadata_moa.x.y, Metadata_moa.x.x)) %>%
      select(-Metadata_moa.x.x, -Metadata_moa.x.y) %>%
      mutate(pass.mean = ifelse(is.na(pass.x), "No", "Yes"),
             pass.mix = ifelse(is.na(pass.y), "No", "Yes")) %>%
      select(-pass.x, -pass.y) %>% 
      arrange(pass.mix) 
    
    d.diff %>%  
      filter(pass.mix == "Yes" & pass.mean == "No") %>%
      arrange(Metadata_moa) %>%
      htmlTable::htmlTable() %>%
      print
    
    d.diff %>%  
      filter(pass.mix == "No" & pass.mean == "Yes") %>%
      arrange(Metadata_moa) %>%
      htmlTable::htmlTable() %>%
      print
  } else {
    l.mean <- lapply(d.mean[3, ], function(x) sum(x)) 
    l.mix <- lapply(d.mix[3, ], function(x) sum(x))
    l.median.mad <- lapply(d.median.mad[3, ], function(x) sum(x)) 
    l.median.mad.2 <- lapply(d.median.mad.2[3, ], function(x) sum(x)) 
    l.median.mad.cov <- lapply(d.median.mad.cov[3, ], function(x) sum(x)) 
    
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
    ggsave("classification_comparison.png", g, width = 7, height = 5)
  }
}

top.prec <- c(seq(from = 0.98, to = 0.997, by = 0.002))
enrichment_top_conn <- Vectorize(enrichment_top_conn, vectorize.args = "top.perc")
sm.mean <- perpare_sm(sm = cr.mean, metadata = metadata)
sm.mix <- perpare_sm(sm = cr.mix, metadata = metadata)
sm.median.mad <- perpare_sm(sm = cr.median.mad, metadata = metadata)
sm.median.mad.2 <- perpare_sm(sm = cr.median.mad.2, metadata = metadata)
sm.median.mad.cov <- perpare_sm(sm = cr.median.mad.cov, metadata = metadata)
FA.sm.mean <- perpare_sm(sm = FA.cr.mean, metadata = metadata)
FA.sm.median <- perpare_sm(sm = FA.cr.median, metadata = metadata)
FA.sm.mad <- perpare_sm(sm = FA.cr.mad, metadata = metadata)
pc.sm.mean <- perpare_sm(sm = pc.cr.mean, metadata = metadata)
pc.sm.median <- perpare_sm(sm = pc.cr.median, metadata = metadata)
pc.sm.mad <- perpare_sm(sm = pc.cr.mad, metadata = metadata)
#dp.sm.mean <- perpare_sm(sm = dp.cr.mean, metadata = metadata)



saveRDS(sm.median.mad, "sm_median_mad.rds")
saveRDS(sm.median.mad.cov, "sm_median_mad_cov.rds")
saveRDS(sm.mix, "sm_median_median.rds")

mean.res <- enrichment_top_conn(sm = sm.mean, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
mix.res <- enrichment_top_conn(sm = sm.mix, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
median.mad.res <- enrichment_top_conn(sm = sm.median.mad, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
median.mad.2.res <- enrichment_top_conn(sm = sm.median.mad.2, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
median.mad.cov.res <- enrichment_top_conn(sm = sm.median.mad.cov, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
FA.mean.res <- enrichment_top_conn(sm = FA.sm.mean, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
FA.median.res <- enrichment_top_conn(sm = FA.sm.median, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
FA.mad.res <- enrichment_top_conn(sm = FA.sm.mad, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
pc.mean.res <- enrichment_top_conn(sm = pc.sm.mean, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
pc.median.res <- enrichment_top_conn(sm = pc.sm.median, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
pc.mad.res <- enrichment_top_conn(sm = pc.sm.mad, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
#dp.mean.res <- enrichment_top_conn(sm = dp.sm.mean, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)

mean.res <- mean.res[3,] %>% unlist %>% unname()
mix.res <- mix.res[3,] %>% unlist %>% unname()
median.mad.res <- median.mad.res[3,] %>% unlist %>% unname()
median.mad.2.res <- median.mad.2.res[3,] %>% unlist %>% unname()
median.mad.cov.res <- median.mad.cov.res[3,] %>% unlist %>% unname()
FA.mean.res <- FA.mean.res[3,] %>% unlist %>% unname()
FA.median.res <- FA.median.res[3,] %>% unlist %>% unname()
FA.mad.res <- FA.mad.res[3,] %>% unlist %>% unname()
pc.mean.res <- pc.mean.res[3,] %>% unlist %>% unname()
pc.median.res <- pc.median.res[3,] %>% unlist %>% unname()
pc.mad.res <- pc.mad.res[3,] %>% unlist %>% unname()
#dp.mean.res <- dp.mean.res[3,] %>% unlist %>% unname()


D1 <- data.frame(top.prec = top.prec * 100, odds.ratio = mean.res, method = "median")
D2 <- data.frame(top.prec = top.prec * 100, odds.ratio = median.mad.res, method = "median+mad (SNF)")
D3 <- data.frame(top.prec = top.prec * 100, odds.ratio = mix.res, method = "median+median (SNF)")
D4 <- data.frame(top.prec = top.prec * 100, odds.ratio = median.mad.2.res, method = "median+mad (concatenated)")
D5 <- data.frame(top.prec = top.prec * 100, odds.ratio = median.mad.cov.res, method = "median+mad+cov. (SNF)")
D6 <- data.frame(top.prec = top.prec * 100, odds.ratio = FA.mean.res, method = "FA.mean")
D7 <- data.frame(top.prec = top.prec * 100, odds.ratio = FA.median.res, method = "FA.median")
D8 <- data.frame(top.prec = top.prec * 100, odds.ratio = FA.mad.res, method = "FA.mad")
D9 <- data.frame(top.prec = top.prec * 100, odds.ratio = pc.mean.res, method = "pc.mean")
D10 <- data.frame(top.prec = top.prec * 100, odds.ratio = pc.median.res, method = "pc.median")
D11 <- data.frame(top.prec = top.prec * 100, odds.ratio = pc.mad.res, method = "pc.mad")
#D12 <- data.frame(top.prec = top.prec * 100, odds.ratio = dp.mean.res, method = "dp.mean")

D <- rbind(D1, D2)
D <- rbind(D, D3)
D <- rbind(D, D4)
D <- rbind(D, D5)
D <- rbind(D, D6)
D <- rbind(D, D7)
D <- rbind(D, D8)
D <- rbind(D, D9)
D <- rbind(D, D10)
D <- rbind(D, D11)
#D <- rbind(D, D12)

#lvls <- sort(unique(as.character(D$method)))
#lvls <- c("median+mad+cov. (SNF)", "median+mad (SNF)", "median+median (SNF)", "median+mad (concatenated)", "median")
lvls <- c("median+mad+cov. (SNF)", "median+mad (SNF)", "median+median (SNF)", "median+mad (concatenated)",
          "median", "FA.mean", "FA.median", "FA.mad", 
          "pc.mean", "pc.median", "pc.mad")
D <- D %>% mutate(method = factor(method, levels = lvls))
D <- D %>% mutate(top.prec = 100 - top.prec)

g <- ggplot(D, aes(x = top.prec, y = odds.ratio, color = method, order = method)) + 
  geom_point() + 
  geom_line() + 
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_continuous(breaks = 100 - rev(top.prec[seq(from = 1, to = length(top.prec), by = 2)] * 100), minor_breaks = 100 - rev(top.prec * 100)) +
  ylab("Folds of enrichment") + 
  xlab("p") +
  ggtitle("Folds of enrichment for top p% connections \n to have same MOAs/Pathways") +
  theme_bw() +
  theme(axis.text = element_text(size=20), text = element_text(size=15)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.title=element_blank())
print(g) 
ggsave("global_comparison.png", g, width = 7, height = 5)
#save.image("workspace.RData")
