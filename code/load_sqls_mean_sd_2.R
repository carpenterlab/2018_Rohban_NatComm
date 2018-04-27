library(dplyr)
library(magrittr)
library(cytominer)
library(foreach)
library(stringr)
library(readr)
library(doParallel)

profile.plate.traditional.2 <- function(pl, project.name, batch.name, operation, nrm.column, nrm.value, cores = 2, feat.list = NULL) {
  
  p1 <- str_split(operation, "\\+")[[1]][1]
  p2 <- str_split(operation, "\\+")[[1]][2]
  
  if (file.exists(paste0("../backend/", batch.name, "/", pl, "/", pl, "_normalized_", p1, "_", p2, ".csv"))) {
    return(NULL)
  }
  
  doParallel::registerDoParallel(cores = cores)
  
  sql.path <- paste0("../backend/", batch.name, "/", pl, "/", pl, ".sqlite")
  out.path <- paste0("../backend/", batch.name, "/", pl)
  out.file <- paste0(out.path, "/", pl, ".csv")
  
  if (!file.exists(sql.path)) {
    system(command = paste0("aws s3 cp 's3://cellpainting-datasets/",
                            project.name, 
                            "/workspace/backend/",
                            batch.name,
                            "/", 
                            pl, 
                            "/", 
                            pl, 
                            ".sqlite' ",
                            sql.path))
  }
  
  if (!file.exists(paste0("../input/", pl, "_normalized.csv"))) {
    system(command = paste0("aws s3 cp 's3://cellpainting-datasets/",
                            project.name, 
                            "/workspace/backend/",
                            batch.name,
                            "/", 
                            pl, 
                            "/", 
                            pl, 
                            "_normalized.csv' ../input/",
                            pl,
                            "_normalized.csv"))
  }
  
  prf <- readr::read_csv(paste0("../input/",
                                pl,
                                "_normalized.csv"))
  
  sites.all <- unique(prf$Metadata_Well)
  
  variables <- colnames(prf)
  variables <- variables[which(!str_detect(variables, "Metadata_"))]
  prf.metadata <- setdiff(colnames(prf), variables)
  if (!is.null(feat.list)) {
    variables <- feat.list
  }
  
  sqlite_file <- sql.path
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
  
  if (!dir.exists("../tmp_trad/")) {
    dir.create("../tmp_trad/")
  }
  
  profiles <- foreach (sites = sites.all, .combine = rbind) %do% {
    
    saveRDS(sites, paste0("../tmp_trad/", sites, ".rds"))  
    
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
    
    profile <- cytominer::aggregate(population = dt.sub, strata = c("Image_Metadata_Well"), variables = variables, operation = operation)
    profile <- cbind(profile, data.frame(Metadata_Plate = pl, Metadata_Well = sites))
    profile
  }
  t2 <- proc.time()
  
  cov.variables <- colnames(profiles)
  cov.variables <- cov.variables[which(str_detect(cov.variables, "Cells_") | str_detect(cov.variables, "Cytoplasm_") | str_detect(cov.variables, "Nuclei_"))]
  
  cov.metadata <- setdiff(colnames(profiles), cov.variables)
  
  dmso.ids <- prf %>% 
    filter((!!rlang::sym(nrm.column)) == nrm.value) %>%
    select(Metadata_Plate, Metadata_Well)
  
  samples.nrm <- profiles %>% 
    mutate(Metadata_Plate = as.character(Metadata_Plate)) %>%
    semi_join(dmso.ids %>%
                mutate(Metadata_Plate = as.character(Metadata_Plate)), 
              by = c("Metadata_Plate", "Metadata_Well"))
  
  mn <- apply(samples.nrm %>% select(one_of(cov.variables)), 2, function(x) mean(x, na.rm = T))
  sdv <- apply(samples.nrm %>% select(one_of(cov.variables)), 2, function(x) sd(x, na.rm = T))
  
  dt <- profiles[, cov.variables]
  
  dt.nrm <- scale(dt, center = mn, scale = sdv)
  
  profiles.nrm <- cbind(dt.nrm, profiles[, cov.metadata])
  
  profiles.nrm.meta <- profiles.nrm[, cov.metadata] %>% 
    mutate(Metadata_Plate = as.character(Metadata_Plate)) %>%
    left_join(prf %>% select(one_of(prf.metadata)) %>%
                mutate(Metadata_Plate = as.character(Metadata_Plate)), by = c("Metadata_Plate", "Metadata_Well"))
  
  profiles.nrm <- cbind(profiles.nrm[, cov.variables], profiles.nrm.meta)
  
  p1 <- str_split(operation, "\\+")[[1]][1]
  p2 <- str_split(operation, "\\+")[[1]][2]
  
  readr::write_csv(profiles.nrm, paste0("../backend/", batch.name, "/", pl, "/", pl, "_normalized_", p1, "_", p2, ".csv"))
  
  system(paste0(paste0("rm ", sql.path)))
}
