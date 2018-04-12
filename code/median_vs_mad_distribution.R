rm(list = ls())
library(dplyr)
library(foreach)
library(stringr)
library(readr)
library(magrittr)
library(SNFtool)
library(ggplot2)

plate.list <- read.csv("../input/processed_plates_CDRP_bio.txt", header = F) %>% as.matrix() %>% as.vector()
feat.list <- read.csv("../input/feature_list.txt", header = F) %>% as.matrix() %>% as.vector()
mean.na <- function(x) {mean(x, na.rm = T)}
metadata.df <- readr::read_csv("../input/metadata_CDRP.csv")  

read.and.summarize <- function(profile.type) {
  fls <- paste0(plate.list, "_covariance.csv")
  feat.list.s <- feat.list
  
  profiles.nrm <- foreach (fl = fls, .combine = rbind) %do% {
    if (profile.type == "cov") {
      if (file.exists(paste0("../output/", fl))) {
        x <- readr::read_csv(paste0("../output/", fl))    
      } else {
        x <- NULL
      }
    } else if (profile.type == "mean") {
      pl <- str_split(fl, "_")[[1]][1]
      x <- readr::read_csv(paste0("../input/", pl, "_normalized.csv"))  
      if (!is.null(feat.list)) {
        x <- x %>%
          select(matches("Metadata_"), one_of(feat.list))
      }
    } else if (profile.type == "median") {
      pl <- str_split(fl, "_")[[1]][1]
      init <- list.dirs("../backend", recursive = F)
      fl.name <- paste0(init, "/", pl, "/", pl, "_normalized_median_mad.csv")
      
      if (file.exists(fl.name)) {
        x <- readr::read_csv(fl.name)    
      } else {
        x <- NULL
        warning(paste0("Plate ", pl, " is missing."))
      }
      
      if (!is.null(feat.list) & ! is.null(x)) {
        x <- x %>%
          select(matches("Metadata_"), one_of(paste0(feat.list, "_median")))
        feat.list.s <- paste0(feat.list, "_median")
      }
    } else if (profile.type == "mad") {
      pl <- str_split(fl, "_")[[1]][1]
      init <- list.dirs("../backend", recursive = F)
      fl.name <- paste0(init, "/", pl, "/", pl, "_normalized_median_mad.csv")
      
      if (file.exists(fl.name)) {
        x <- readr::read_csv(fl.name)    
      } else {
        x <- NULL
        warning(paste0("Plate ", pl, " is missing."))
      }
      
      if (!is.null(feat.list) & ! is.null(x)) {
        x <- x %>%
          select(matches("Metadata_"), one_of(paste0(feat.list, "_mad")))
        feat.list.s <- paste0(feat.list, "_mad")
      }
    }
    else if (str_detect(profile.type, "\\+")) {
      p1 <- str_split(profile.type, "\\+")[[1]][1]  
      p2 <- str_split(profile.type, "\\+")[[1]][2]  
      
      if (!is.null(feat.list)) {
        p1 <- str_split(profile.type, "\\+")[[1]][1]  
        p2 <- str_split(profile.type, "\\+")[[1]][2]  
        feat.list.s <- c(paste0(feat.list, "_", p1), paste0(feat.list, "_", p2))
      }
      
      pl <- str_split(fl, "_")[[1]][1]
      init <- list.dirs("../backend", recursive = F)
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
  
  
  prf <- profiles.nrm
  feats <- colnames(prf)
  feats <- feats[which(!str_detect(feats, "Metadata_"))]
  
  return(list(data = profiles.nrm, feats = feats))
}

prf.median <- read.and.summarize("median")
prf.mad <- read.and.summarize("mad")

plt <- function(x, nm) {
  x <- x[x > quantile(x, 0.01) & x < quantile(x, 0.99)]
  x <- data.frame(col = x)
  g <- ggplot(x, aes(x = col)) + geom_histogram(aes(y = ..density..), bins = 30) + geom_density() + xlab(nm) + theme(axis.text = element_text(size=20), text = element_text(size=25)) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  return(g)
}

g1 <- plt(prf.median$data$Nuclei_Texture_InfoMeas2_DNA_5_0_median, "Nuclei_Texture_InfoMeas2_DNA_5_0_median")
ggsave("median_hist.png", g1)

g2 <- plt(prf.mad$data$Nuclei_Texture_InfoMeas2_DNA_5_0_mad, "Nuclei_Texture_InfoMeas2_DNA_5_0_mad")
ggsave("mad_hist.png", g2)
