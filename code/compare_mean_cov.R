rm(list = ls())

library(dplyr)
library(ggplot2)

source("moa_evaluations.R")

enrichment.based.classification <- FALSE
k.snf <- 7     # neighborhood size in SNF
k <- 1:10      # k top hits are used for classification
not.same.batch <- T

cr.melt.mean <- readRDS("cr_mean.rds")
cr.melt.cov <- readRDS("cr_cov.rds")
cr.melt.median.mad <- readRDS("cr_median+mad.rds")

cr.mean <- cr.melt.mean %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2") 

cr.median.mad <- cr.melt.median.mad %>%
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

cm.rn <- setdiff(intersect(intersect(rownames(cr.median.mad), rownames(cr.mean)), rownames(cr.cov)), NA)
cr.mean <- cr.mean[cm.rn, cm.rn]
cr.cov <- cr.cov[cm.rn, cm.rn]
cr.median.mad <- cr.median.mad[cm.rn, cm.rn]

cr.mean[is.na(cr.mean)] <- 0
cr.median.mad[is.na(cr.median.mad)] <- 0
cr.cov[is.na(cr.cov)] <- 0

af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = 0.5)
af.2 <- SNFtool::affinityMatrix(Diff = 1 - cr.cov, K = k.snf, sigma = 0.5)
af.snf <- SNFtool::SNF(list(af.1, af.2), K = k.snf, t = 10)
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.mix <- af.snf

metadata <- cr.melt.mean %>%
  select(Var1, Metadata_moa.x, Metadata_Plate_Map_Name.x) %>%
  unique() %>%
  mutate(Metadata_broad_sample = Var1, Metadata_moa = Metadata_moa.x, Metadata_Plate_Map_Name = Metadata_Plate_Map_Name.x) %>%
  select(-Var1, -Metadata_moa.x, -Metadata_Plate_Map_Name.x)

cmpd_classification <- Vectorize(cmpd_classification, "k0")
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
      ylab("No. of compounds with a \n same MOA compound in their k-NNs")
    print(g)    
    ggsave("classification_comparison.png", g, width = 7, height = 5)
  }
} else {
  d.mean <- cmpd_knn_classification(cr.mean, metadata, k, not.same.batch = not.same.batch) 
  d.mix <- cmpd_knn_classification(cr.mix, metadata, k, not.same.batch = not.same.batch) 
  d.median.mad <- cmpd_knn_classification(cr.median.mad, metadata, k, not.same.batch = not.same.batch) 
  
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
    
    D <- data.frame(method = "mean", k = k, tp = (unlist(l.mean)))
    D <- rbind(D, 
               data.frame(method = "mean+cov.", k = k, tp = (unlist(l.mix))))
    D <- rbind(D, 
               data.frame(method = "median+mad", k = k, tp = (unlist(l.median.mad))))
    D <- D %>% mutate(method = factor(method, levels = sort(unique(as.character(D$method)))))
    
    g <- ggplot(D, aes(x = k, y = tp, color = method, order = as.character(method))) + 
      geom_point() + 
      geom_line() + 
      scale_y_continuous(limits = c(0, NA)) +
      scale_x_continuous(breaks = k, minor_breaks = k) +
      ylab("No. of compounds with a \n same MOA compound in their k-NNs")
    print(g) 
    ggsave("classification_comparison.png", g, width = 7, height = 5)
  }
}

top.prec <- seq(from = 0.98, to = 0.999, by = 0.002)
enrichment_top_conn <- Vectorize(enrichment_top_conn, vectorize.args = "top.perc")
mean.res <- enrichment_top_conn(sm = cr.mean, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
mix.res <- enrichment_top_conn(sm = cr.mix, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
median.mad.res <- enrichment_top_conn(sm = cr.median.mad, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)

mean.res <- mean.res["estimate",] %>% unlist %>% unname()
mix.res <- mix.res["estimate",] %>% unlist %>% unname()
median.mad.res <- median.mad.res["estimate",] %>% unlist %>% unname()

D1 <- data.frame(top.prec = top.prec * 100, odds.ratio = mean.res, method = "mean")
D2 <- data.frame(top.prec = top.prec * 100, odds.ratio = median.mad.res, method = "median+mad")
D3 <- data.frame(top.prec = top.prec * 100, odds.ratio = mix.res, method = "mean+cov.")

D <- rbind(D1, D2)
D <- rbind(D, D3)
D <- D %>% mutate(method = factor(method, levels = sort(unique(as.character(D$method)))))
D <- D %>% mutate(top.prec = 100 - top.prec)

g <- ggplot(D, aes(x = top.prec, y = odds.ratio, color = method, order = method)) + 
  geom_point() + 
  geom_line() + 
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_continuous(breaks = 100 - rev(top.prec[seq(from = 1, to = length(top.prec), by = 2)] * 100), minor_breaks = 100 - rev(top.prec * 100)) +
  ylab("No. of folds of enrichment \n for top p% conn. to have same MOA") + 
  xlab("p")
print(g) 
ggsave("global_comparison.png", g, width = 7, height = 5)
