rm(list = ls())

library(dplyr)
library(ggplot2)

source("moa_evaluations.R")

enrichment.based.classification <- FALSE
k.snf <- 7     # neighborhood size in SNF
k <- 1:10      # k top hits are used for classification

cr.melt.mean <- readRDS("cr_mean.rds")
cr.melt.cov <- readRDS("cr_cov.rds")

cr.mean <- cr.melt.mean %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2") 

cr.cov <- cr.melt.cov %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2")

d <- apply(cr.mean, 1, function(x) !(sum(is.na(x)) >= (NROW(cr.mean) -1 )))
cr.mean <- cr.mean[d, d]

d <- apply(cr.cov, 1, function(x) !(sum(is.na(x)) >= (NROW(cr.cov) -1 )))
cr.cov <- cr.cov[d, d]

cm.rn <- setdiff(intersect(rownames(cr.mean), rownames(cr.cov)), NA)
cr.mean <- cr.mean[cm.rn, cm.rn]
cr.cov <- cr.cov[cm.rn, cm.rn]

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

cmpd_classification <- Vectorize(cmpd_classification, "k0")
cmpd_knn_classification <- Vectorize(cmpd_knn_classification, "k0")

if (enrichment.based.classification) {
  d.mean <- cmpd_classification(sm = cr.mean, metadata = metadata, k0 = k)
  d.mix <- cmpd_classification(sm = cr.mix, metadata = metadata, k0 = k)
  
  if (length(k) == 1) {
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
      ylab("No. of compounds with a \n same MOA compound in their k-NNs")
    print(g)    
    ggsave("classification_comparison.png", g, width = 7, height = 5)
  }
} else {
  d.mean <- cmpd_knn_classification(cr.mean, metadata, k) 
  d.mix <- cmpd_knn_classification(cr.mix, metadata, k) 
  
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
      ylab("No. of compounds with a \n same MOA compound in their k-NNs")
    print(g) 
    ggsave("classification_comparison.png", g, width = 7, height = 5)
  }
}
