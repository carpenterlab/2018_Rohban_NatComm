## plot median vs covariance profile distributions
rm(list = ls())

library(dplyr)
library(magrittr)
library(cytominer)
library(foreach)
library(stringr)
library(readr)
library(doParallel)

pl <- "24278"
nrm.column <- "Metadata_broad_sample"
nrm.value <- "DMSO"

project.name <- "CDRPBIO-BBBC036-Bray"
batch.name <- "CDRP"

dir.create("../tmp_sql")

system(command = paste0("aws s3 cp 's3://cellpainting-datasets/",
                        project.name, 
                        "/workspace/backend/",
                        batch.name,
                        "/", 
                        pl, 
                        "/", 
                        pl, 
                        ".sqlite' ../tmp_sql/",
                        pl,
                        ".sqlite"))

system(command = paste0("aws s3 cp 's3://cellpainting-datasets/",
                        project.name, 
                        "/workspace/backend/",
                        batch.name,
                        "/", 
                        pl, 
                        "/", 
                        pl, 
                        "_normalized.csv' ../tmp_sql/",
                        pl,
                        "_normalized.csv"))

feat.list <- readr::read_csv("../input/feature_list.txt", col_names = F) %>% as.matrix() %>% as.vector()

prf <- readr::read_csv(paste0("../tmp_sql/",
                              pl,
                              "_normalized.csv"))

sites.all <- unique(prf$Metadata_Well)

variables <- colnames(prf)
variables <- variables[which(!str_detect(variables, "Metadata_"))]
prf.metadata <- setdiff(colnames(prf), variables)
if (!is.null(feat.list)) {
  variables <- feat.list
}

sqlite_file <- paste0("../tmp_sql/", pl, ".sqlite")
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

saveRDS(profiles, paste0("../tmp_sql/cov_profiles_", pl, ".rds"))