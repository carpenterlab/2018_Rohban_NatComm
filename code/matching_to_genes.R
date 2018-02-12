rm(list = ls())

library(dplyr)
library(stringr)
library(foreach)
library(readr)
library(magrittr)
source("rep.corr.func.R")

k.snf <- 7
feat.list <- readr::read_csv("../input/feature_list.txt", col_names = F)
feat.list <- feat.list %>% unlist() %>% unname()

get.gene <- function(x) {str_split(x, "_")[[1]][1]}
get.allele <- function(x) {str_split(x, "_")[[1]][2]}

get.gene <- Vectorize(get.gene)
get.allele <- Vectorize(get.allele)

metadata.cdrp <- readr::read_csv("../input/metadata_CDRP.csv")
metadata.ta <- readr::read_csv("../input/metadata_TA.csv")
metadata.ta %<>% 
  mutate(Metadata_Name = `Gene Allele Name`) %>%
  mutate(Metadata_Gene = get.gene(Metadata_Name), Metadata_Allele = get.allele(Metadata_Name)) %>%
  filter(str_detect(Metadata_Allele, "WT")) 

read.and.summarize <- function(profile.type, path, feat.list, metadata.df) {
  fls <- list.files(paste0(path, "/output"))
  
  profiles.nrm <- foreach (fl = fls, .combine = rbind) %do% {
    if (profile.type == "cov") {
      x <- readr::read_csv(paste0(path, "/output/", fl))  
    } else {
      pl <- str_split(fl, "_")[[1]][1]
      x <- readr::read_csv(paste0(path, "/input/", pl, "_normalized.csv"))  
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
      variable.names <- feat.list
    }
    
    profiles.nrm <- profiles.nrm %>% select(one_of(c(meta.cols, variable.names)))
  }
  
  print(length(fls))
  
  if (!"Metadata_mmoles_per_liter" %in% colnames(profiles.nrm)) {
    profiles.nrm <- profiles.nrm %>%
      mutate(Metadata_mmoles_per_liter = 10)
  }
  
  thr <- non.rep.cor(list(data = profiles.nrm), c("Metadata_broad_sample", "Metadata_mmoles_per_liter"), variable.names, 0.95)
  u <- rep.cor(list(data = profiles.nrm), c("Metadata_broad_sample", "Metadata_mmoles_per_liter"), variable.names)
  
  prf <- profiles.nrm %>% 
    group_by(Metadata_broad_sample, Metadata_mmoles_per_liter) %>%
    summarise_at(.vars = variable.names, .funs = "mean")
  
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
  
  return(list(profile = list(data = profiles.nrm, feats = feats), thr = thr, rep.cor = u))
}

get.weak.trt <- function(Pf) {
  u <- Pf$rep.cor %>%
    filter(cr < Pf$thr) %>%
    select(Metadata_broad_sample) %>%
    as.matrix() %>%
    as.vector()
  
  u <- setdiff(u, NA)
  return(u)
}

Pf.1.mean <- read.and.summarize("mean", "../../CDRP_2nd_moment", feat.list, metadata.cdrp)
Pf.2.mean <- read.and.summarize("mean", "../../TA_ORF_2nd_moment", feat.list, metadata.ta)

Pf.1.cov <- read.and.summarize("cov", "../../CDRP_2nd_moment", feat.list, metadata.cdrp)
Pf.2.cov <- read.and.summarize("cov", "../../TA_ORF_2nd_moment", feat.list, metadata.ta)

compute.corr <- function(Pf.1, Pf.2) {
  cr <- cor(Pf.1$profile$data[, Pf.1$profile$feats] %>% t, Pf.2$profile$data[, Pf.2$profile$feats] %>% t)
  rownames(cr) <- Pf.1$profile$data$Metadata_broad_sample
  colnames(cr) <- Pf.2$profile$data$Metadata_broad_sample
  
  w.1 <- get.weak.trt(Pf.1)
  w.2 <- get.weak.trt(Pf.2)
  
  #cr[setdiff(w.1, NA), ] <- NA
  #cr[, setdiff(w.2, NA)] <- NA
  
  return(cr)
}

cr.1.2.mean <- compute.corr(Pf.1.mean, Pf.2.mean)
cr.1.2.cov <- compute.corr(Pf.1.cov, Pf.2.cov)

m11 <- matrix(NA, NROW(Pf.1.mean$profile$data), NROW(Pf.1.mean$profile$data))
m22 <- matrix(NA, NROW(Pf.2.mean$profile$data), NROW(Pf.2.mean$profile$data))
rownames(m11) <- rownames(cr.1.2.mean)
colnames(m11) <- rownames(cr.1.2.mean)
rownames(m22) <- colnames(cr.1.2.mean)
colnames(m22) <- colnames(cr.1.2.mean)

c1 <- rbind(cbind(m11, cr.1.2.mean), cbind(t(cr.1.2.mean), m22))

m11 <- matrix(NA, NROW(Pf.1.cov$profile$data), NROW(Pf.1.cov$profile$data))
m22 <- matrix(NA, NROW(Pf.2.cov$profile$data), NROW(Pf.2.cov$profile$data))
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


cr.1 <- cr.1.2.mean 
cr.2 <- cr.1.2.cov 
cr.3 <- cr.mix[1:NROW(cr.1), (NROW(cr.1) + 1):NCOL(cr.mix)]

  #pmax(cr.1/(sd(cr.1, na.rm = T)), cr.2/(sd(cr.2, na.rm = T)), na.rm = T)

gene_enr <- function(cr, metadata.1, metadata.2, quant) {
  is.in.list.one <- function(gene.list, gene) {
    gene %in% str_split(gene.list, ", ")[[1]]
  }
  
  is.in.list <- function(gene.list, gene) {
    lapply(gene, function(x) (x %in% str_split(gene.list, ", ")[[1]])) %>% unlist %>% any
  }
  is.in.list.one <- Vectorize(is.in.list.one)
  is.in.list.multi <- Vectorize(is.in.list, vectorize.args = "gene.list")
  
  cr.annot <- cr %>% 
    reshape2::melt() %>%
    filter(!is.na(value)) %>%
    filter(Var1 %in% cmpd) %>%
    left_join(metadata.1, by = c("Var1" = "Metadata_broad_sample")) %>%
    left_join(metadata.2, by = c("Var2" = "Metadata_broad_sample")) %>%
    filter(!is.na(Metadata_target) & !is.na(Metadata_Gene)) %>%
    select(Var1, Var2, Metadata_Gene, Metadata_target, value) %>%
    mutate(consistent = is.in.list.one(Metadata_target, Metadata_Gene))
  
  res <- cr.annot %>%
    filter(consistent) %>%
    filter(value > quantile(value, quant, na.rm = T)) 
  
  V <- rbind(c(cr.annot %>% filter(consistent & value > quantile(value, quant, na.rm = T)) %>% NROW, 
               cr.annot %>% filter(!consistent & value > quantile(value, quant, na.rm = T)) %>% NROW),
             c(cr.annot %>% filter(consistent & value < quantile(value, quant, na.rm = T)) %>% NROW, 
               cr.annot %>% filter(!consistent & value < quantile(value, quant, na.rm = T)) %>% NROW))
  
  list(V = V,
    test = fisher.test(V, alternative = "greater"),
       res.top = cr.annot %>%
         filter(value > quantile(value, quant)), 
        res.valid = cr.annot %>%
         filter(consistent & value > quantile(value, quant)),
       res = res,
    no.cmpds = length(cmpd))
}

w.1.mean <- get.weak.trt(Pf.1.mean)
w.1.cov <- get.weak.trt(Pf.1.cov)
w.2.mean <- get.weak.trt(Pf.2.mean)
w.2.cov <- get.weak.trt(Pf.2.cov)

cr.1[w.1.mean, ] <- NA
cr.1[, w.2.mean] <- NA

cr.2[w.1.cov, ] <- NA
cr.2[, w.2.cov] <- NA

cr.3[intersect(w.1.mean, w.1.cov), ] <- NA
cr.3[, intersect(w.2.mean, w.2.cov)] <- NA

a1 <- gene_enr(cr.1, metadata.cdrp, metadata.ta, 0.95)
a2 <- gene_enr(cr.2, metadata.cdrp, metadata.ta, 0.95)
a3 <- gene_enr(cr.3, metadata.cdrp, metadata.ta, 0.95)

a1$test
a2$test
a3$test

a1$res.valid %>% htmlTable::htmlTable()
a2$res.valid %>% htmlTable::htmlTable()
a3$res.valid %>% htmlTable::htmlTable()
