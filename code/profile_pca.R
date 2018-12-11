#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = T)

extends <- methods::extends

'profile_pca
Usage:
profile_pca -n <project_name> -b <batch_name> -p <plate_list_path> -f <feat_list_path> -l <norm_column> -v <norm_value>

Options:
-h --help                                         Show this screen.
-n <project_name> --name=<project_name>           Project name on s3.
-b <batch_name> --batch=<batch_name>              Batch name. 
-p <plate_list_path> --plate=<plate_list_path>    Path of the plate list.
-f <feat_list_file> --feats=<feat_list_file>      Path to the file containing the list of features. 
-l <norm_column> --col=<norm_column>              Column name to be used to select samples for normalization.
-v <norm_value> --value=<norm_value>              Value of the mentioned column which indicates the sample.
' -> doc

opts <- docopt::docopt(doc)

plate.list.path <- opts[["plate"]] 
batch.name <- opts[["batch"]] 
feat.list.path <- opts[["feats"]] 
project.name <- opts[["name"]]
col.name <- opts[["col"]]
col.val <- opts[["value"]]

library(dplyr)
library(magrittr)
library(cytominer)
library(foreach)
library(stringr)
library(readr)
library(iterators)
library(doParallel)
library(reshape2)
library(psych)
library(ggplot2)
library(readbulk)

feat.list <- read.table(feat.list.path, header = F) %>% as.matrix() %>% as.vector() %>% unlist() %>% unname()
plate.list <- read.table(plate.list.path, header = F) %>% as.matrix() %>% as.vector() %>% unlist()
variables <- feat.list

dir.create("../PCA/", recursive = T)

# Combining all the metadata information from normalized csv of all plates
f.path <- NULL
for (p in 1:length(plate.list)) {
  f.path[p]<- paste0("../input/", plate.list[p], "_normalized.csv")
}

# reading sqlite
read_sql<- function(sql.path) {
  db <- DBI::dbConnect(RSQLite::SQLite(), sql.path)
  RSQLite::initExtension(db)
  
  image <- RSQLite::dbReadTable(conn = db, "Image")
  cells <- RSQLite::dbReadTable(conn = db, "Cells")
  nuclei <- RSQLite::dbReadTable(conn = db, "Nuclei")
  cytoplasm <- RSQLite::dbReadTable(conn = db, "Cytoplasm")
  
  dt <- cells %>%
    left_join(cytoplasm, by = c("TableNumber", "ImageNumber", "ObjectNumber")) %>%
    left_join(nuclei, by = c("TableNumber", "ImageNumber", "ObjectNumber")) %>%
    left_join(image, by = c("TableNumber", "ImageNumber"))
  
  return(dt)
  
}

dmso.all <- data.frame(readr::read_csv("../FA/dmso/dmso.all.csv"),  stringsAsFactors =F)

dmeta <- colnames(dmso.all)[which(str_detect(colnames(dmso.all), "Metadata"))]
#variables <- setdiff(colnames(dmso.all), dmeta)

# Standardization step
mn <- apply(dmso.all %>% select(one_of(variables)) %>% as.matrix(), 2, function(x) mean(x, na.rm = T))
sdv <- apply(dmso.all %>% select(one_of(variables)) %>% as.matrix(), 2, function(x) sd(x, na.rm = T))

# doing the Principle component analysis

pca <- prcomp(dmso.all[, variables], center=mn, scale=sdv)
pc <- (pca$rotation[,1:50])

profile <- NULL
for (i in 1:length(plate.list)) {
  prf <- read.csv(f.path[i], stringsAsFactors = FALSE)

  # Extracting metadata
  meta <- colnames(prf)[which(str_detect(colnames(prf), "Metadata_"))]
  pmeta <- prf %>% select(one_of(meta)) %>% dplyr::collect()
  
  # reading sqlite file
  pl <- plate.list[i]
  sql.path <- paste0("../backend/", batch.name, "/", pl, "/", pl, ".sqlite")
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
  sql_data <- as.data.frame(read_sql(sql.path))
  
  image.col <- sql_data %>%
    select(Image_Metadata_Well, Image_Metadata_Plate) %>%
    colnames()
  sql_data <- sql_data %>%
    select(image.col, variables) %>%
    dplyr::collect()
  
  # removing NAs
  sql_data[is.na(sql_data)] <- 0
  sql_data <- merge(sql_data, pmeta,
                    by.x = c("Image_Metadata_Plate","Image_Metadata_Well"),
                    by.y = c("Metadata_Plate", "Metadata_Well"))
  
  metadata <- colnames(sql_data)[which(str_detect(colnames(sql_data), "Metadata_"))]
  metadata <- sql_data %>% select(one_of(metadata)) %>% dplyr::collect()
  
  dmso <- sql_data %>%
    filter((!!rlang::sym(col.name)) == col.val) %>%
    dplyr::collect()
  
  #system(paste0(paste0("rm ", sql.path)))
  
  # calculating mean and standard_deviation of control condition
  mn <- apply(dmso %>% select(one_of(variables)), 2, function(x) mean(x, na.rm = T))
  sdv <- apply(dmso %>% select(one_of(variables)), 2, function(x) sd(x, na.rm = T))
  
  sql_data <- scale(sql_data[, variables], center = mn, scale = sdv)
  
  # prediction of PCs for test dataset
  prediction <- sql_data %*% pc
  
  # converting matrix into dataframe
  prediction <- as.data.frame(prediction)
  
  profiles.nrm  <- cbind(metadata, prediction)
  
  # selecting variables
  pc_variables <- colnames(profiles.nrm)
  pc_metavariables <- pc_variables[which(str_detect(pc_variables, "Metadata"))]
  pc_variables <- setdiff(colnames(profiles.nrm), pc_metavariables)
  
  operation <- c("mean", "median", "mad")
  
  for (l in 1:length(operation)) {
    profiles.nrm[, pc_variables] <- apply(profiles.nrm[, pc_variables],
                                          2, function(x) as.numeric(x))
    profile <- cytominer::aggregate(population = profiles.nrm,
                                    strata = c("Image_Metadata_Plate", "Image_Metadata_Well"),
                                    variables = pc_variables,
                                    operation = operation[l])
    
    
    profile <- cbind(profile, pmeta)
    readr::write_csv(profile, paste0("../PCA/", plate.list[i], "_PC_", operation[l], ".csv"))
    
  }
}