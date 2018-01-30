#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = T)

extends <- methods::extends

'profile
Usage:
profile -p <plate_number> -d <no_component> -r <random_density> -c <cores>
Options:
-h --help                                         Show this screen.
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

profile.plate(pl = pl, n.components = p, rand.density = rand.density, cores = cores)  
