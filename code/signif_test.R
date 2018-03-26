rm(list = ls())

library(dplyr)
library(ggplot2)
library(stringr)

load("workspace.RData")

signif.test <- function(sm1, sm2, top.perc, not.same.batch = F) {
  if (not.same.batch) {
    sm1 <- sm1 %>%
      filter((is.na(Metadata_Plate_Map_Name.x) & !is.na(Metadata_Plate_Map_Name.y))
             | (is.na(Metadata_Plate_Map_Name.y) & !is.na(Metadata_Plate_Map_Name.x))
             | (Metadata_Plate_Map_Name.x != Metadata_Plate_Map_Name.y))
    sm2 <- sm2 %>%
      filter((is.na(Metadata_Plate_Map_Name.x) & !is.na(Metadata_Plate_Map_Name.y))
             | (is.na(Metadata_Plate_Map_Name.y) & !is.na(Metadata_Plate_Map_Name.x))
             | (Metadata_Plate_Map_Name.x != Metadata_Plate_Map_Name.y))
  }
  
  thr1 <- quantile(sm1$value, top.perc, na.rm = T)
  thr2 <- quantile(sm2$value, top.perc, na.rm = T)
  
  v11 <- sm1 %>%
    filter(value > thr1 & same.moa) %>%
    NROW
  
  v12 <- sm1 %>%
    filter(value < thr1 & same.moa) %>%
    NROW
  
  v21 <- sm2 %>%
    filter(value > thr2 & same.moa) %>%
    NROW
  
  v22 <- sm2 %>%
    filter(value < thr2 & same.moa) %>%
    NROW
  
  V <- rbind(c(v11, v12), c(v21, v22))
  rownames(V) <- c("first data", "second data")
  colnames(V) <- c("similar - validated", "non-similar - validated")
  print(V)
  
  f1 <- (fisher.test(x = V, 
                     alternative = "greater"))

  
  v11 <- sm1 %>%
    filter(value > thr1 & same.moa) %>%
    NROW
  
  v12 <- sm1 %>%
    filter(value > thr1 & !same.moa) %>%
    NROW
  
  v21 <- sm2 %>%
    filter(value > thr2 & same.moa) %>%
    NROW
  
  v22 <- sm2 %>%
    filter(value > thr2 & !same.moa) %>%
    NROW
  
  V <- rbind(c(v11, v12), c(v21, v22))
  rownames(V) <- c("first data", "second data")
  colnames(V) <- c("similar - validated", "similar - non-validated")
  print(V)
  
  f2 <- (fisher.test(x = V, 
                     alternative = "greater"))
  return(list(f1, f2))
}

signif.test(sm.median.mad.cov, sm.median.mad.2, top.perc = 0.995, not.same.batch = not.same.batch)
