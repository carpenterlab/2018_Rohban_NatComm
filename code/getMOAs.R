library(dplyr)
library(stringr)
library(RCurl)
library(readr)

y <- readr::read_csv("cdrp_missing_brds.txt", col_names = F) %>% as.matrix() %>% as.vector()
D <- data.frame(Metadata_broad_sample = c(), Metadata_moa = c())

for (yi in y) {
  yj <- str_sub(yi, 1, 13)
  x <- RCurl::getURL(paste0("https://api.clue.io/api/perts?filter={%22where%22:{%22pert_id%22:%22", yj, "%22}}&user_key=c0768ab9f945a12687e03d459d577c3a"))
  if (x != "[]") {
    moa <- stringr::str_split(stringr::str_split(stringr::str_split(x, "moa")[[1]][2], "\\:\\[")[[1]][2], "\"")[[1]][2]
    D <- rbind(D, data.frame(Metadata_broad_sample = yi, Metadata_moa = moa))
  }
}

D2 <- D %>% filter(!is.na(Metadata_moa))

readr::write_csv(D2, "moa_cdrp_added.csv")
