library(cytotools)

profile.plate.traditional <- function(pl, project.name, batch.name, operation) {
  
  doParallel::registerDoParallel(cores = cores)
  
  sql.path <- paste0("../backend/", batch.name, "/", pl, "/", pl, ".sqlite")
  out.path <- paste0("../backend/", batch.name, "/", pl)
  out.file <- paste0(out.path, "/", pl, ".csv")
  
  if (!file.exists(sql.path)) {
    system(command = paste0("aws s3 cp 's3://imaging-platform/projects/",
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
  
  cytotools::aggregate(sqlite_file = sql.path, 
                       output_file = out.file,
                       operation = operation,
                       variables = "all") 
  
  system(paste0(paste0("rm ", sql.path)))
}
