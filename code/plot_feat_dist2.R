rm(list = ls())
gc()

library(dplyr)
library(ggplot2)

profiles <- readRDS("../tmp_sql/cov_profiles_24278.rds")

plt <- function(x, ft.name) {
  x <- x[x > quantile(x, 0.01) & x < quantile(x, 0.99)]
  x <- data.frame(col = x)
  g <- ggplot(x, aes(x = col)) + geom_histogram(aes(y = ..density..), bins = 30) + geom_density()
  
  return(g)
}

y <- readr::read_csv("../backend/CDRP/24278/24278_normalized_median_mad.csv")
x <- profiles$Nuclei_Texture_Variance_RNA_10_0__Nuclei_Texture_InverseDifferenceMoment_RNA_10_0
z <- y$Nuclei_Texture_InverseDifferenceMoment_RNA_10_0_median
w <- y$Nuclei_Texture_Variance_RNA_10_0_median

x <- x[which(!is.na(x))]
z <- z[which(!is.na(z))]
w <- w[which(!is.na(w))]

plt(x, "cov(Cells_AreaShape_Compactness, Cells_AreaShape_Area)")
plt(z, "median(Cells_Neighbors_AngleBetweenNeighbors_Adjacent)")
plt(w, "median(Cells_AreaShape_Area)")
