#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = T)

extends <- methods::extends

'sampe_dmso
Usage:
sampe_dmso -n <project_name> -b <batch_name> -p <plate_list_path> -f <feat_list_path> -l <norm_column> -v <norm_value>

Options:
-h --help                                         Show this screen.
-n <project_name> --name=<project_name>           Project name.
-b <batch_name> --batch=<batch_name>              Batch name. 
-p <plate_list_path> --plate=<plate_list_path>    Path of the plate list.
-f <feat_list_path> --feats=<feat_list_path>      Path to the file containing the list of features. 
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
library(xtable)
library(ggplot2)
library(readbulk)

feat.list <- read.table(feat.list.path, header = F) %>% as.matrix() %>% as.vector() %>% unlist() %>% unname()
plate.list <- read.table(plate.list.path, header = F) %>% as.matrix() %>% as.vector() %>% unlist()
variables <- feat.list

dir.create("../FA/dmso/", recursive = T)

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

# Combining all the metadata information from normalized csv of all plates
for (p in 1:length(plate.list)) {
  f.path<- paste0("../input/", plate.list[p], "_normalized.csv")
  
  pl <- plate.list[p]
  
  dmso <- NULL
  
  prf <- read.csv(f.path, stringsAsFactors = FALSE)
  
  # Extracting metadata
  meta <- colnames(prf)[which(str_detect(colnames(prf), "Metadata_"))]
  pmeta <- prf %>% select(one_of(meta)) %>% dplyr::collect()
  
  # reading sqlite file
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
  
  sql_data <- read_sql(sql.path)
  image.col <- sql_data %>% select(Image_Metadata_Well, Image_Metadata_Plate) %>% colnames()
  sql_data <- sql_data %>% select(image.col, variables) %>% dplyr::collect()
  
  # removing NAs
  sql_data[is.na(sql_data)] <- 0
  sql_data <- merge(sql_data, pmeta,
                    by.x = c("Image_Metadata_Plate","Image_Metadata_Well"),
                    by.y = c("Metadata_Plate", "Metadata_Well"))
  
  metadata <- colnames(sql_data)[which(str_detect(colnames(sql_data), "Metadata_"))]
  dmso <- sql_data %>%
    filter((!!rlang::sym(col.name)) == col.val) %>%
    dplyr::collect()
  
  set.seed(123)
  sample_size <- 455
  train_indx <- sample(seq_len(nrow(dmso)), size = sample_size)
  dmso <- dmso[train_indx, ]
  readr::write_csv(dmso, paste0("../FA/dmso/", plate.list[p], "_dmso", ".csv"))
  #system(paste0(paste0("rm ", sql.path)))
}

# combining all dmso cells collected from each replicate plate
path <- "../FA/dmso/"
dmso.all <- read_bulk(directory = path, extension = "_dmso.csv", stringsAsFactors=FALSE)
readr::write_csv(dmso.all, paste0("../FA/dmso/", "dmso.all", ".csv"))