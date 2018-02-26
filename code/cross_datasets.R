rm(list = ls())

library(dplyr)
library(stringr)
library(foreach)
library(readr)
library(magrittr)
library(ggplot2)
source("moa_evaluations.R")

mean.na <- function(x) {mean(x, na.rm = T)}

k.snf <- 15
feat.list <- readr::read_csv("../input/feature_list.txt", col_names = F)
feat.list <- feat.list %>% unlist() %>% unname()

metadata.cdrp <- readr::read_csv("../input/metadata_CDRP.csv")
metadata.repurp <- NULL

read.and.summarize <- function(profile.type, path, feat.list, metadata.df) {
  fls <- list.files(paste0(path, "/output"))
  feat.list.s <- feat.list
  
  profiles.nrm <- foreach (fl = fls, .combine = rbind) %do% {
    if (profile.type == "cov") {
      x <- readr::read_csv(paste0(path, "/output/", fl))  
    } else if (profile.type == "mean") {
      pl <- str_split(fl, "_")[[1]][1]
      x <- readr::read_csv(paste0(path, "/input/", pl, "_normalized.csv"))  
      x <- x %>% select(matches("Metadata_"), one_of(feat.list))
    } else if (profile.type == "median") {
      pl <- str_split(fl, "_")[[1]][1]
      init <- list.dirs(paste0(path, "/backend"), recursive = F)
      fl.name <- paste0(init, "/", pl, "/", pl, "_normalized_median_mad.csv")
      
      if (file.exists(fl.name)) {
        x <- readr::read_csv(fl.name)    
        if (!is.null(feat.list)) {
          x <- x %>%
            select(matches("Metadata_"), one_of(paste0(feat.list, "_median")))
          feat.list.s <- paste0(feat.list, "_median")
        }
      } else {
        x <- NULL
        warning(paste0("Plate ", pl, " is missing."))
      }
    } else if (str_detect(profile.type, "\\+")) {
      p1 <- str_split(profile.type, "\\+")[[1]][1]  
      p2 <- str_split(profile.type, "\\+")[[1]][2]  
      
      if (!is.null(feat.list)) {
        p1 <- str_split(profile.type, "\\+")[[1]][1]  
        p2 <- str_split(profile.type, "\\+")[[1]][2]  
        feat.list.s <- c(paste0(feat.list, "_", p1), paste0(feat.list, "_", p2))
      }
      
      pl <- str_split(fl, "_")[[1]][1]
      init <- list.dirs(paste0(path, "/backend"), recursive = F)
      fl.name <- paste0(init, "/", pl, "/", pl, "_normalized_", p1, "_", p2, ".csv")
      if (file.exists(fl.name)) {
        x <- readr::read_csv(fl.name)    
      } else {
        x <- NULL
        warning(paste0("Plate ", pl, " is missing."))
      }
    }
    
    x
  }
  
  variable.names <- colnames(profiles.nrm)
  variable.names <- variable.names[which(!str_detect(variable.names, "Metadata_"))]

  # in some special cases of mean profiles (e.g. CDRP), there seems to be Inf, and NA value.
  # this is to treat those cases
  if (profile.type != "cov") {
    meta.cols <- setdiff(colnames(profiles.nrm), variable.names)
    
    if (is.null(feat.list)) {
      ids <- apply(profiles.nrm[,variable.names], 2, function(x) !any(is.na(x) | is.nan(x) | is.infinite(x) | sd(x) > 10)) %>% which
      variable.names <- variable.names[ids]
    } else {
      variable.names <- feat.list.s
    }
    
    profiles.nrm <- profiles.nrm %>% select(one_of(c(meta.cols, variable.names)))
  }
  
  print(length(fls))
  
  if (!"Metadata_mmoles_per_liter" %in% colnames(profiles.nrm)) {
    profiles.nrm <- profiles.nrm %>%
      mutate(Metadata_mmoles_per_liter = 10)
  }
  
  prf <- profiles.nrm %>% 
    group_by(Metadata_broad_sample, Metadata_mmoles_per_liter) %>%
    summarise_at(.vars = variable.names, .funs = "mean.na")
  
  prf %<>%
    arrange(abs(Metadata_mmoles_per_liter - 10)) %>%
    group_by(Metadata_broad_sample) %>%
    slice(1) %>%
    ungroup()
  
  if (is.null(metadata.df)) {
    prf %<>% 
      left_join(profiles.nrm %>% select(Metadata_broad_sample, Metadata_moa) %>% unique, by = "Metadata_broad_sample")
  } else {
    prf %<>% 
      left_join(metadata.df %>% select(Metadata_broad_sample, Metadata_moa) %>% unique, by = "Metadata_broad_sample")
  }
  
  profiles.nrm <- prf
  feats <- colnames(prf)
  feats <- feats[which(!str_detect(feats, "Metadata_"))]
  
  return(list(data = profiles.nrm, feats = feats))
}

metadata.cdrp <- metadata.cdrp %>%
  mutate(Metadata_broad_sample = str_sub(Metadata_broad_sample, 1, 13))

match.brds <- function(cr, k) {
  tot <- cr %>% 
    reshape2::melt() %>%
    mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
    group_by(Var1) %>%
    summarise(pass = any(Var1 == Var2)) %>%
    ungroup() %>%
    select(pass) %>%
    as.matrix() %>%
    as.vector() %>%
    sum
  
  ps <- cr %>% 
    reshape2::melt() %>%
    mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
    arrange(-value) %>%
    group_by(Var1) %>%
    slice(1:k) %>%
    summarise(pass = any(Var1 == Var2)) %>%
    ungroup() %>%
    select(pass) %>%
    as.matrix() %>%
    as.vector() %>%
    sum
  
  return(ps/tot)
}

Pf.1.mean <- read.and.summarize("median", "../../CDRP_2nd_moment", feat.list, metadata.cdrp)
Pf.2.mean <- read.and.summarize("median", "../../repurp_2nd_moment", feat.list, metadata.repurp)

Pf.1.cov <- read.and.summarize("cov", "../../CDRP_2nd_moment", feat.list, metadata.cdrp)
Pf.2.cov <- read.and.summarize("cov", "../../repurp_2nd_moment", feat.list, metadata.repurp)

Pf.1.median.mad <- read.and.summarize("median+mad", "../../CDRP_2nd_moment", feat.list, metadata.cdrp)
Pf.2.median.mad <- read.and.summarize("median+mad", "../../repurp_2nd_moment", feat.list, metadata.repurp)

Pf.1.mean$data %<>% 
  mutate(Metadata_Treatment = str_sub(Metadata_broad_sample, 1, 13))

Pf.2.mean$data %<>% 
  mutate(Metadata_Treatment = str_sub(Metadata_broad_sample, 1, 13))

Pf.1.cov$data %<>% 
  mutate(Metadata_Treatment = str_sub(Metadata_broad_sample, 1, 13))

Pf.2.cov$data %<>% 
  mutate(Metadata_Treatment = str_sub(Metadata_broad_sample, 1, 13))

Pf.1.median.mad$data %<>% 
  mutate(Metadata_Treatment = str_sub(Metadata_broad_sample, 1, 13))

Pf.2.median.mad$data %<>% 
  mutate(Metadata_Treatment = str_sub(Metadata_broad_sample, 1, 13))

compute.corr <- function(Pf.1, Pf.2) {
  cr <- cor(Pf.1$data[, Pf.1$feats] %>% t, Pf.2$data[, Pf.2$feats] %>% t)
  rownames(cr) <- Pf.1$data$Metadata_Treatment
  colnames(cr) <- Pf.2$data$Metadata_Treatment
  
  return(cr)
}


v1 <- Pf.1.mean$data %>%
  select(Metadata_Treatment) %>%
  unique %>%
  as.matrix() %>%
  as.vector()
v2 <- Pf.1.cov$data %>%
  select(Metadata_Treatment) %>%
  unique %>%
  as.matrix() %>%
  as.vector()
v3 <- Pf.1.median.mad$data %>%
  select(Metadata_Treatment) %>%
  unique %>%
  as.matrix() %>%
  as.vector()

v <- intersect(intersect(v1, v2), v3)

w1 <- Pf.2.mean$data %>%
  select(Metadata_Treatment) %>%
  unique %>%
  as.matrix() %>%
  as.vector()
w2 <- Pf.2.cov$data %>%
  select(Metadata_Treatment) %>%
  unique %>%
  as.matrix() %>%
  as.vector()
w3 <- Pf.2.median.mad$data %>%
  select(Metadata_Treatment) %>%
  unique %>%
  as.matrix() %>%
  as.vector()

w <- intersect(intersect(w1, w2), w3)

Pf.1.mean$data %<>% filter(Metadata_Treatment %in% v)
Pf.1.cov$data %<>% filter(Metadata_Treatment %in% v)
Pf.1.median.mad$data %<>% filter(Metadata_Treatment %in% v)

Pf.2.mean$data %<>% filter(Metadata_Treatment %in% w)
Pf.2.cov$data %<>% filter(Metadata_Treatment %in% w)
Pf.2.median.mad$data %<>% filter(Metadata_Treatment %in% w)

cr.1.mean <- compute.corr(Pf.1.mean, Pf.1.mean)
cr.1.cov <- compute.corr(Pf.1.cov, Pf.1.cov)

cr.2.mean <- compute.corr(Pf.2.mean, Pf.2.mean)
cr.2.cov <- compute.corr(Pf.2.cov, Pf.2.cov)

cr.1.2.median.mad <- compute.corr(Pf.1.median.mad, Pf.2.median.mad)
cr.1.2.mean <- compute.corr(Pf.1.mean, Pf.2.mean)
cr.1.2.cov <- compute.corr(Pf.1.cov, Pf.2.cov)

cr.2.1.mean <- compute.corr(Pf.2.mean, Pf.1.mean)
cr.2.1.cov <- compute.corr(Pf.2.cov, Pf.1.cov)

cr.mean <- rbind(cbind(cr.1.mean, cr.1.2.mean), cbind(cr.2.mean, cr.2.1.mean))
cr.cov <- rbind(cbind(cr.1.cov, cr.1.2.cov), cbind(cr.2.cov, cr.2.1.cov))

m11 <- matrix(NA, NROW(Pf.1.mean$data), NROW(Pf.1.mean$data))
m22 <- matrix(NA, NROW(Pf.2.mean$data), NROW(Pf.2.mean$data))
rownames(m11) <- rownames(cr.1.2.mean)
colnames(m11) <- rownames(cr.1.2.mean)
rownames(m22) <- colnames(cr.1.2.mean)
colnames(m22) <- colnames(cr.1.2.mean)

c1 <- rbind(cbind(m11, cr.1.2.mean), cbind(t(cr.1.2.mean), m22))

m11 <- matrix(NA, NROW(Pf.1.cov$data), NROW(Pf.1.cov$data))
m22 <- matrix(NA, NROW(Pf.2.cov$data), NROW(Pf.2.cov$data))
rownames(m11) <- rownames(cr.1.2.cov)
colnames(m11) <- rownames(cr.1.2.cov)
rownames(m22) <- colnames(cr.1.2.cov)
colnames(m22) <- colnames(cr.1.2.cov)

c2 <- rbind(cbind(m11, cr.1.2.cov), cbind(t(cr.1.2.cov), m22))

c1[(is.na(c1) | is.nan(c1) | is.infinite(c1))] <- -2
c2[(is.na(c2) | is.nan(c2) | is.infinite(c2))] <- -2

af.1 <- SNFtool::affinityMatrix(Diff = 1 - c1, K = k.snf, sigma = 0.5)
af.2 <- SNFtool::affinityMatrix(Diff = 1 - c2, K = k.snf, sigma = 0.5)
af.snf <- SNFtool::SNF(list(af.1, af.2), K = k.snf, t = 10)
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.mix <- af.snf

 match.brds <- Vectorize(match.brds, vectorize.args = "k")

rng <- 1:30 

cr.1 <- cr.mean[1:NROW(Pf.1.mean$data), (NROW(Pf.1.mean$data) + 1):NCOL(cr.mean)]
cr.2 <- cr.cov[1:NROW(Pf.1.mean$data), (NROW(Pf.1.mean$data) + 1):NCOL(cr.cov)]
cr.3 <- cr.mix[1:NROW(Pf.1.mean$data), (NROW(Pf.1.mean$data) + 1):NCOL(cr.mix)]
cr.4 <- cr.1.2.median.mad

v1 <- match.brds(cr.mean[1:NROW(Pf.1.mean$data), (NROW(Pf.1.mean$data) + 1):NCOL(cr.mean)], rng) %>% unlist
v2 <- match.brds(cr.3, rng) %>% unlist
v3 <- match.brds(cr.4, rng) %>% unlist

D1 <- data.frame(method = "mean", k = rng, y = v1)
D2 <- data.frame(method = "mean+cov.", k = rng, y = v2)
D3 <- data.frame(method = "median+mad", k = rng, y = v3)

D <- rbind(D1, D2)
D <- rbind(D, D3)
D <- D %>% mutate(method = factor(method, levels = sort(unique(as.character(D$method)))))

g <- ggplot2::ggplot(D, aes(x = k, y = y, color = method, order = method)) + geom_line() + geom_point() + ylab("how many compounds find themselves \n within their kNN in the other batch")
g

metadata.repurp <- Pf.2.cov$data %>%
  select(matches("Metadata_")) %>%
  unique %>%
  mutate(Metadata_broad_sample = str_sub(Metadata_broad_sample, 1, 13))

top.precs <- seq(from = 0.98, to = 0.999, by = 0.002)
a1 <- lapply(top.precs, function(x) enrichment_top_conn_cross(sm = cr.1, metadata1 = metadata.cdrp, metadata2 = metadata.repurp, top.perc = x)$test$estimate) %>% unlist
a3 <- lapply(top.precs, function(x) enrichment_top_conn_cross(sm = cr.3, metadata1 = metadata.cdrp, metadata2 = metadata.repurp, top.perc = x)$test$estimate) %>% unlist
a4 <- lapply(top.precs, function(x) enrichment_top_conn_cross(sm = cr.4, metadata1 = metadata.cdrp, metadata2 = metadata.repurp, top.perc = x)$test$estimate) %>% unlist


D1 <- data.frame(p = top.precs, folds = a1, method = "median")
D2 <- data.frame(p = top.precs, folds = a3, method = "median+cov.")
D3 <- data.frame(p = top.precs, folds = a4, method = "median+mad")
D <- rbind(D1, D2)
D <- rbind(D, D3)
D <- D %>% mutate(method = factor(method, levels = sort(unique(as.character(D$method)))))
D <- D %>% mutate(p = 1 - p)

g <- ggplot(D, aes(x = p * 100, y = folds, color = method)) + 
  geom_point() + 
  geom_line() + 
  scale_y_continuous(limits = c(0, NA)) + 
  scale_x_continuous(breaks = 100 - rev(top.precs[seq(from = 1, to = length(top.precs), by = 2)] * 100), minor_breaks = 100 - rev(top.precs * 100)) +
  ylab("No. of folds of enrichment \n for top p% conn. to have same MOA") + 
  xlab("p")
print(g) 
ggsave("global_comparison_across.png", g, width = 7, height = 5)


g
