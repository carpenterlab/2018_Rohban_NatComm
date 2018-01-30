library(dplyr)
library(magrittr)
library(cytominer)
library(foreach)
library(stringr)
library(readr)
library(doParallel)

source("generate_component_matrix.R")

profile.plate <- function(pl, n.components = 3000, rand.density = 0.1, cores = 4) {
    
  doParallel::registerDoParallel(cores = cores)
  
  if (!file.exists(paste0("../input/", pl, ".sqlite"))) {
    system(command = paste0("aws s3 cp 's3://imaging-platform/projects/2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad/workspace/backend/2016_04_01_a549_48hr_batch1/", 
                            pl, 
                            "/", 
                            pl, 
                            ".sqlite' ../input/",
                            pl,
                            ".sqlite"))
  }
  
  if (!file.exists(paste0("../input/", pl, "_normalized_variable_selected.csv"))) {
    system(command = paste0("aws s3 cp 's3://imaging-platform/projects/2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad/workspace/backend/2016_04_01_a549_48hr_batch1/", 
                            pl, 
                            "/", 
                            pl, 
                            "_normalized_variable_selected.csv' ../input/",
                            pl,
                            "_normalized_variable_selected.csv"))
  }
  
  prf <- readr::read_csv(paste0("../input/",
                          pl,
                         "_normalized_variable_selected.csv"))
  
  sites.all <- unique(prf$Metadata_Well)
  
  variables <- colnames(prf)
  variables <- variables[which(!str_detect(variables, "Metadata_"))]
  prf.metadata <- setdiff(colnames(prf), variables)
  
  sqlite_file <- paste0("../input/", pl, ".sqlite")
  db <- DBI::dbConnect(RSQLite::SQLite(), sqlite_file)
  RSQLite::initExtension(db)
  
  image_object_join_columns <- c("TableNumber", "ImageNumber")
  strata <- c("Image_Metadata_Plate", "Image_Metadata_Well")
  
  image <- dplyr::tbl(src = db, "image") %>%
    dplyr::select(c(image_object_join_columns, strata))
  cells <- dplyr::tbl(src = db, "cells")
  cytoplasm <- dplyr::tbl(src = db, "cytoplasm")
  nuclei <- dplyr::tbl(src = db, "nuclei")
  
  image.coll <- image %>%
    select(Image_Metadata_Well, ImageNumber, TableNumber) %>%
    dplyr::collect()
  
  t <- proc.time()
  profiles <- foreach (sites = sites.all, .combine = rbind) %do% {
    
    saveRDS(sites, paste0("../tmp/", sites, ".rds"))  
    
    image.sub <- image.coll %>% 
      filter(Image_Metadata_Well == sites) 
    
    append_operation_tag <- function(s) stringr::str_c(s, operation, sep = "_")
    
    dt.sub <- foreach (i = 1:NROW(image.sub), .combine = rbind) %dopar% {
      cells.sub <- cells %>% 
        filter(ImageNumber == image.sub$ImageNumber[i] & TableNumber == image.sub$TableNumber[i]) %>%
        dplyr::collect()
      
      cytoplasm.sub <- cytoplasm %>% 
        filter(ImageNumber == image.sub$ImageNumber[i] & TableNumber == image.sub$TableNumber[i]) %>%
        dplyr::collect()
      
      nuclei.sub <- nuclei %>% 
        filter(ImageNumber == image.sub$ImageNumber[i] & TableNumber == image.sub$TableNumber[i]) %>%
        dplyr::collect()
      
      dt <- cells.sub %>%
        dplyr::inner_join(cytoplasm.sub, by = "ObjectNumber") %>%
        dplyr::inner_join(nuclei.sub, by = "ObjectNumber") %>%
        dplyr::inner_join(image.sub, by = image_object_join_columns) 
        
      dt
    }
    
    all.variables <- colnames(dt.sub)
    all.variables <- all.variables[which(str_detect(all.variables, "Cells") | str_detect(all.variables, "Cytoplasm") | str_detect(all.variables, "Nuclei"))]
    metadata <- setdiff(colnames(dt.sub), all.variables)
    
    dt.sub <- dt.sub %>%
      select(one_of(c(metadata, variables)))
    
    dt.sub[, variables] <- apply(dt.sub[, variables], 2, function(x) as.numeric(x))
    
    profile <- cytominer::covariance(population = dt.sub, variables = variables)
    profile <- cbind(profile, data.frame(Metadata_Plate = pl, Metadata_Well = sites))
    profile
  }
  t2 <- proc.time()
  
  cov.variables <- colnames(profiles)
  cov.variables <- cov.variables[which(str_detect(cov.variables, "Cells_") | str_detect(cov.variables, "Cytoplasm_") | str_detect(cov.variables, "Nuclei_"))]
  cov.metadata <- setdiff(colnames(profiles), cov.variables)
  
  dmso.ids <- prf %>% 
    filter(Metadata_broad_sample == "DMSO") %>%
    select(Metadata_Plate, Metadata_Well)
  
  samples.nrm <- profiles %>%
    semi_join(dmso.ids, by = c("Metadata_Plate", "Metadata_Well"))
  
  mn <- apply(samples.nrm %>% select(one_of(cov.variables)), 2, function(x) mean(x, na.rm = T))
  sdv <- apply(samples.nrm %>% select(one_of(cov.variables)), 2, function(x) sd(x, na.rm = T))
  
  dt <- profiles[, cov.variables]
  
  dt.nrm <- scale(dt, center = mn, scale = sdv)
  
  profiles.nrm <- cbind(dt.nrm, profiles[, cov.metadata])
  
  profiles.nrm.meta <- profiles.nrm[, cov.metadata] %>% 
    left_join(prf %>% select(one_of(prf.metadata)), by = c("Metadata_Plate", "Metadata_Well"))
  
  profiles.nrm <- cbind(profiles.nrm[, cov.variables], profiles.nrm.meta)
  cov.metadata <- setdiff(colnames(profiles.nrm), cov.variables)
  
  if (!file.exists("../input/random_projection.rds")) {
    rand.proj <- generate_component_matrix(n_features = length(cov.variables), n_components = n.components, density = rand.density)
    saveRDS(rand.proj, "../input/random_projection.rds")
  } else {
    rand.proj <- readRDS("../input/random_projection.rds")
  }
  
  dt <- as.matrix(profiles.nrm[, cov.variables])
  dt[is.na(dt)] <- 0
  
  profiles.nrm.red <- dt %*% as.matrix(rand.proj)
  profiles.nrm <- cbind(as.data.frame(profiles.nrm.red), profiles.nrm[, cov.metadata])
  
  readr::write_csv(profiles.nrm, paste0("../output/", pl, "_covariance.csv"))
  
  system(paste0("rm ../input/", pl, ".sqlite"))
}
