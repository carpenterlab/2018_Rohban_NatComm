#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = T)

extends <- methods::extends

'profile
Usage:
profile -n <project_name> -b <batch_name> -p <plate_number> -d <no_component> -r <random_density> -c <cores>
Options:
-h --help                                         Show this screen.
-n <project_name> --name=<project_name>           Project name on s3.
-b <batch_name> --batch=<batch_name>              Batch name. 
-p <plate_number> --plate=<plate_number>          Plate number.
-d <no_component> --dim=<no_component>            Number of random projections.
-r <random_density> --rdensity=<random_density>   Density of non-zero values in sparse random projections.
-c <cores> --cores=<cores>                        Number of cores to parallelize. ' -> doc

opts <- docopt::docopt(doc)

source("load_sqls.R")

p <- as.numeric(opts[["dim"]])
rand.density <- as.numeric(opts[["rdensity"]])
cores <- as.numeric(opts[["cores"]])
pl <- opts[["plate"]]
proj.name <- opts[["name"]]
batch.name <- opts[["batch"]]

profile.plate(pl = pl, project.name = proj.name, batch.name = batch.name, n.components = p, rand.density = rand.density, cores = cores)  
