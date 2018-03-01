rm(list = ls())

library(dplyr)
library(ggplot2)

source("moa_evaluations.R")

enrichment.based.classification <- FALSE
k.snf <- 7     # neighborhood size in SNF
t <- 10
k <- 1:10      # k top hits are used for classification
not.same.batch <- F

cr.melt.mean <- readRDS("cr_median.rds")
cr.melt.cov <- readRDS("cr_cov.rds")
cr.melt.mad <- readRDS("cr_mad.rds")  
cr.melt.median.mad.2 <- readRDS("cr_median+mad.rds") 

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

cr.mean <- sim_normalize(cr.mean) 
cr.cov <- sim_normalize(cr.cov) 
cr.median.mad.2 <- sim_normalize(cr.median.mad.2) 
cr.mad <- sim_normalize(cr.mad) 

metad <- readr::read_csv("../input/metadata_TA_2.csv")
gene.s <- function(x) {
  str_split(x, "_")[[1]][1]
}
allele.s <- function(x) {
  str_split(x, "_")[[1]][2]
}
gene.s <- Vectorize(gene.s)
allele.s <- Vectorize(allele.s)

subs <- metad %>%
  mutate(gene = gene.s(Metadata_Treatment),
         allele = allele.s(Metadata_Treatment))  %>%
  filter(str_detect(allele, "WT")) %>%
  group_by(gene) %>%
  slice(1) %>%
  ungroup() %>%
  select(Metadata_broad_sample) %>%
  as.matrix() %>%
  as.vector()

cr.mean <- cr.mean[subs, subs]

d <- apply(cr.mean, 1, function(x) !(sum(is.na(x)) >= (NROW(cr.mean) -1 )))
cr.mean <- cr.mean[d, d]

d <- apply(cr.cov, 1, function(x) !(sum(is.na(x)) >= (NROW(cr.cov) -1 )))
cr.cov <- cr.cov[d, d]

d <- apply(cr.mad, 1, function(x) !(sum(is.na(x)) >= (NROW(cr.mad) -1 )))
cr.mad <- cr.mad[d, d]

cm.rn <- setdiff(intersect(intersect(rownames(cr.mad), rownames(cr.mean)), intersect(rownames(cr.cov), rownames(cr.median.mad.2))), NA)
cr.mean <- cr.mean[cm.rn, cm.rn]
cr.cov <- cr.cov[cm.rn, cm.rn]
cr.mad <- cr.mad[cm.rn, cm.rn]
cr.median.mad.2 <- cr.median.mad.2[cm.rn, cm.rn]

cr.mean[is.na(cr.mean)] <- 0
cr.median.mad.2[is.na(cr.median.mad.2)] <- 0
cr.mad[is.na(cr.mad)] <- 0
cr.cov[is.na(cr.cov)] <- 0

af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = sigma.mean)
af.2 <- SNFtool::affinityMatrix(Diff = 1 - cr.mad, K = k.snf, sigma = sigma.mad)
af.snf <- SNFtool::SNF(list(af.1, af.2), K = k.snf, t = t)
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.median.mad <- af.snf

af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = sigma.mean)
af.2 <- SNFtool::affinityMatrix(Diff = 1 - cr.cov, K = k.snf, sigma = sigma.cov)
af.snf <- SNFtool::SNF(list(af.1, af.2), K = k.snf, t = t)
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.mix <- af.snf

af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = sigma.mean)
af.2 <- SNFtool::affinityMatrix(Diff = 1 - cr.mad, K = k.snf, sigma = sigma.mad)
af.3 <- SNFtool::affinityMatrix(Diff = 1 - cr.cov, K = k.snf, sigma = sigma.cov)
af.snf <- SNFtool::SNF(list(af.1, af.2, af.3), K = k.snf, t = t)
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.median.mad.cov <- af.snf

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
      ylab("No. of genes with a \n same Pathway gene in their k-NNs")
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
    #D <- rbind(D, 
    #           data.frame(method = "median+cov. (SNF)", k = k, tp = (unlist(l.mix))))
    D <- rbind(D, 
               data.frame(method = "median+mad (SNF)", k = k, tp = (unlist(l.median.mad))))
    D <- rbind(D, 
               data.frame(method = "median+mad (concatenated)", k = k, tp = (unlist(l.median.mad.2))))
    D <- rbind(D, 
               data.frame(method = "median+mad+cov. (SNF)", k = k, tp = (unlist(l.median.mad.cov))))
    
    lvls <- c("median+mad+cov. (SNF)", "median+mad (SNF)", "median+mad (concatenated)", "median")
    D <- D %>% mutate(method = factor(method, levels = lvls))
    
    g <- ggplot(D, aes(x = k, y = tp, color = method, order = as.character(method))) + 
      geom_point() + 
      geom_line() + 
      scale_y_continuous(limits = c(0, NA)) +
      scale_x_continuous(breaks = k, minor_breaks = k) +
      ylab("No. of genes with a \n same Pathway gene in their k-NNs")
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

mean.res <- enrichment_top_conn(sm = sm.mean, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
mix.res <- enrichment_top_conn(sm = sm.mix, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
median.mad.res <- enrichment_top_conn(sm = sm.median.mad, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
median.mad.2.res <- enrichment_top_conn(sm = sm.median.mad.2, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)
median.mad.cov.res <- enrichment_top_conn(sm = sm.median.mad.cov, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)


mean.res <- mean.res[3,] %>% unlist %>% unname()
mix.res <- mix.res[3,] %>% unlist %>% unname()
median.mad.res <- median.mad.res[3,] %>% unlist %>% unname()
median.mad.2.res <- median.mad.2.res[3,] %>% unlist %>% unname()
median.mad.cov.res <- median.mad.cov.res[3,] %>% unlist %>% unname()

D1 <- data.frame(top.prec = top.prec * 100, odds.ratio = mean.res, method = "median")
D2 <- data.frame(top.prec = top.prec * 100, odds.ratio = median.mad.res, method = "median+mad (SNF)")
D3 <- data.frame(top.prec = top.prec * 100, odds.ratio = mix.res, method = "median+cov. (SNF)")
D4 <- data.frame(top.prec = top.prec * 100, odds.ratio = median.mad.2.res, method = "median+mad (concatenated)")
D5 <- data.frame(top.prec = top.prec * 100, odds.ratio = median.mad.cov.res, method = "median+mad+cov. (SNF)")

D <- rbind(D1, D2)
#D <- rbind(D, D3)
D <- rbind(D, D4)
D <- rbind(D, D5)

#lvls <- sort(unique(as.character(D$method)))
lvls <- c("median+mad+cov. (SNF)", "median+mad (SNF)", "median+mad (concatenated)", "median")
D <- D %>% mutate(method = factor(method, levels = lvls))
D <- D %>% mutate(top.prec = 100 - top.prec)

g <- ggplot(D, aes(x = top.prec, y = odds.ratio, color = method, order = method)) + 
  geom_point() + 
  geom_line() + 
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_continuous(breaks = 100 - rev(top.prec[seq(from = 1, to = length(top.prec), by = 2)] * 100), minor_breaks = 100 - rev(top.prec * 100)) +
  ylab("No. of folds of enrichment \n for top p% conn. to have same Pathway") + 
  xlab("p")
print(g) 
ggsave("global_comparison.png", g, width = 7, height = 5)