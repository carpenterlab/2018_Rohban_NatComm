#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = T)

extends <- methods::extends

'evaluate
Usage:
method -m <method_name> -e <metadata_file>
Options:
-h --help                                         Show this screen.
-m <method_name> --method=<method_name>           Profiling method name, which could be mean or cov.
-e <metadata_file> --metadata=<metadata_file>     Path to a csv file containing the association between the Metadata_broad_sample and Metadata_moa. This could be skipped if it is present in the profiles.
' -> doc

opts <- docopt::docopt(doc)

p <- opts[["method"]]
meta.file <- opts[["metadata"]]

library(dplyr)
library(foreach)
library(stringr)
library(readr)
library(magrittr)
library(SNFtool)
library(ggplot2)

if (!is.null(meta.file)) {
  metadata.df <- readr::read_csv(meta.file)  
} else {
  metadata.df <- NULL
}


type.eval <- "global" # global or classification or lift
profile.type <- p
print(p)
mix1 <- "mean"
mix2 <- "cov"
quant <- 0.99

read.and.summarize <- function(profile.type) {
  fls <- list.files("../output")
  
  profiles.nrm <- foreach (fl = fls, .combine = rbind) %do% {
    if (profile.type == "cov") {
      x <- readr::read_csv(paste0("../output/", fl))  
    } else {
      pl <- str_split(fl, "_")[[1]][1]
      x <- readr::read_csv(paste0("../input/", pl, "_normalized_variable_selected.csv"))  
    }
    x
  }
  
  variable.names <- colnames(profiles.nrm)
  variable.names <- variable.names[which(!str_detect(variable.names, "Metadata_"))]

  # in some special cases of mean profiles (e.g. CDRP), there seems to be Inf, and NA value.
  # this is to treat those cases
  if (profile.type != "cov") {
    meta.cols <- setdiff(colnames(profiles.nrm), variable.names)
    ids <- apply(profiles.nrm[,variable.names], 2, function(x) !any(is.na(x) | is.nan(x) | is.infinite(x) | sd(x) > 10)) %>% which
    variable.names <- variable.names[ids]
    profiles.nrm <- profiles.nrm %>% select(one_of(c(meta.cols, variable.names)))
  }
  
  print(length(fls))
  
  if (!"Metadata_mmoles_per_liter" %in% colnames(profiles.nrm)) {
    profiles.nrm <- profiles.nrm %>%
      mutate(Metadata_mmoles_per_liter = 10)
  }
  
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
  
  return(list(data = profiles.nrm, feats = feats))
}

evaluate.moa <- function(cr, profiles.meta, quant = 0.95, type.eval = "global", k = 1) {
  cr.melt <- cr %>% reshape2::melt()
  
  cr.melt <- cr.melt %>% left_join(profiles.meta, by = c("Var1" = "Metadata_broad_sample")) %>% left_join(profiles.meta, by = c("Var2" = "Metadata_broad_sample")) 
  
  match.moas <- function(moa1, moa2) {
    if (is.na(moa1) | is.na(moa2)) {
      return(FALSE)
    }
    x <- str_split(moa1, "\\|")[[1]]
    y <- str_split(moa2, "\\|")[[1]]
    return(any(x %in% y) | any(y %in% x))
  }
  match.moas <- Vectorize(match.moas)
  saveRDS(cr.melt, paste0("cr_",  profile.type, ifelse(profile.type == "mix", paste0(mix1, mix2), ""), ".rds"))
  
  if (type.eval == "global") {
    cr.melt <- cr.melt %>% filter(Var1 < Var2)
    
    thr <- quantile(cr.melt$value, quant, na.rm = T)
    
    v11 <- cr.melt %>% filter(value > thr & match.moas(Metadata_moa.x, Metadata_moa.y)) %>% NROW()
    v12 <- cr.melt %>% filter(value > thr & !match.moas(Metadata_moa.x, Metadata_moa.y)) %>% NROW()
    v21 <- cr.melt %>% filter(value <= thr & match.moas(Metadata_moa.x, Metadata_moa.y)) %>% NROW()
    v22 <- cr.melt %>% filter(value <= thr & !match.moas(Metadata_moa.x, Metadata_moa.y)) %>% NROW()
    
    V <- rbind(c(v11, v12), c(v21, v22))
    
    return(fisher.test(V, alternative = "greater")$estimate)
  } else if (type.eval == "classification") {
    res <- cr.melt %>%
      filter(Var1 != Var2) %>% 
      arrange(-value) %>%
      group_by(Var1, Metadata_moa.x) %>%
      slice(1:k) %>%
      summarise(good = any(match.moas(Metadata_moa.x, Metadata_moa.y))) %>%
      ungroup() %>%
      filter(good) 
    return(res)
  } else if (type.eval == "lift") {
    same.moa <- match.moas
    
    u <- cr.melt %>%
      filter(as.character(Var1) < as.character(Var2) &
               Var1 != "DMSO" & Var2 != "DMSO" & !is.na(Var1) & !is.na(Var2)) %>%
      arrange(-value) %>%
      mutate(same.moa = same.moa(Metadata_moa.x, Metadata_moa.y))
      
      x <- u$same.moa %>% as.numeric() %>% cumsum()
      x <- x/x[length(x)]
      y <- (!u$same.moa) %>% as.numeric() %>% cumsum()
      y <- y/y[length(y)]
      
      u <- data.frame(x - y)
      
      colnames(u) <- profile.type
      u <- cbind(data.frame(n = 1:NROW(u)), u)
      
      return(u)
  }
}

if (profile.type != "mix") {
  Pf <- read.and.summarize(profile.type = profile.type)
  profiles.nrm <- Pf$data
  feats <- Pf$feats
  if (!is.null(metadata.df)) {
    profiles.meta <- metadata.df %>% select("Metadata_broad_sample", "Metadata_moa") %>% unique
  } else {
    profiles.meta <- profiles.nrm %>% select("Metadata_broad_sample", "Metadata_moa") %>% unique
  }
  
  cr <- cor(profiles.nrm[, feats] %>% t)
  rownames(cr) <- profiles.nrm$Metadata_broad_sample
  colnames(cr) <- profiles.nrm$Metadata_broad_sample
  
  res <- evaluate.moa(cr = cr, profiles.meta = profiles.meta, quant = quant, type.eval = type.eval)
  if (type.eval == "classification") {
    print(NROW(res))
    saveRDS(res, paste0("res_",  profile.type, "_classification.rds"))
  } else {
    print(res)  
  }
  
  if (type.eval == "lift") {
    saveRDS(res, paste0("lift_",  profile.type, ".rds"))
  }
} else {
  Pf.1 <- read.and.summarize(profile.type = mix1)
  
  if (mix2 != mix1) {
    Pf.2 <- read.and.summarize(profile.type = mix2)
  } else {
    Pf.2 <- Pf.1
  }
    
  profiles.nrm.1 <- Pf.1$data
  feats.1 <- Pf.1$feats
  if (!is.null(metadata.df)) {
    profiles.meta <- metadata.df %>% select("Metadata_broad_sample", "Metadata_moa") %>% unique
  } else {
    profiles.meta <- profiles.nrm %>% select("Metadata_broad_sample", "Metadata_moa") %>% unique
  }
  
  cr.1 <- cor(profiles.nrm.1[, feats.1] %>% t)
  rownames(cr.1) <- profiles.nrm.1$Metadata_broad_sample
  colnames(cr.1) <- profiles.nrm.1$Metadata_broad_sample
  
  profiles.nrm.2 <- Pf.2$data
  feats.2 <- Pf.2$feats
  
  cr.2 <- cor(profiles.nrm.2[, feats.2] %>% t)
  rownames(cr.2) <- profiles.nrm.2$Metadata_broad_sample
  colnames(cr.2) <- profiles.nrm.2$Metadata_broad_sample
  
  cr.1 <- cr.1[rownames(cr.2), colnames(cr.2)]
  
  af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.1, K = 7, sigma = 0.5)
  af.2 <- SNFtool::affinityMatrix(Diff = 1 - cr.2, K = 7, sigma = 0.5)
  
  
  af.snf <- SNFtool::SNF(list(af.1, af.2), K = 7, t = 10)
  rownames(af.snf) <- rownames(af.1)
  colnames(af.snf) <- colnames(af.1)
  
  print(evaluate.moa(cr = af.snf, profiles.meta = profiles.meta, quant = quant, type.eval = type.eval))
}
